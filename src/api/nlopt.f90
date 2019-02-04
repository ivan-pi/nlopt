module nlopt_user_func_mod

    use, intrinsic :: iso_c_binding, only: c_double, c_int
    implicit none
    private

    public :: nlopt_user_func
    public :: nlopt_user_mfunc
    public :: nlopt_user_precond

    !
    ! these types are to be extended by the user
    !

    type, abstract :: nlopt_user_func
    contains
        procedure(nlopt_func_if), deferred, public :: eval
    end type

    type, abstract :: nlopt_user_mfunc
    contains
        procedure(nlopt_mfunc_if), deferred, public :: eval
    end type

    type, abstract :: nlopt_user_precond
    contains
        procedure(nlopt_precond_if), deferred, public :: eval
    end type


    abstract interface
        real(c_double) function nlopt_func_if(this,n,x,grad)
            import nlopt_user_func, c_double, c_int
            class(nlopt_user_func), intent(in) :: this
            integer(c_int), intent(in) :: n
            real(c_double), intent(in) :: x(n)
            real(c_double), intent(out), optional :: grad(n)
        end function
        subroutine nlopt_mfunc_if(this,m,result,n,x,grad)
            import nlopt_user_mfunc, c_double, c_int
            class(nlopt_user_mfunc), intent(in) :: this
            integer(c_int), intent(in) :: m
            real(c_double), intent(out) :: result(m)
            integer(c_int), intent(in) :: n
            real(c_double), intent(in) :: x(n)
            real(c_double), intent(out), optional :: grad(m,n)
        end subroutine
        subroutine nlopt_precond_if(this,n,x,v,vpre)
            import nlopt_user_precond, c_double, c_int
            class(nlopt_user_precond), intent(in) :: this
            integer(c_int), intent(in) :: n
            real(c_double), intent(in) :: x(n)
            real(c_double), intent(in) :: v(n)
            real(c_double), intent(out) :: vpre(n)
        end subroutine
    end interface

end module

module adaptor_mod

    use iso_c_binding
    use nlopt_user_func_mod
    implicit none
    private

    public :: adaptor, func_aux, mfunc_aux, precond_aux

    type :: adaptor
        class(nlopt_user_func),    pointer :: f   => null()
        class(nlopt_user_mfunc),   pointer :: mf  => null()
        class(nlopt_user_precond), pointer :: pre => null()
    contains
        final :: destroy
    end type

    type :: adaptor_list
        type(adaptor) :: me
        type(adaptor), pointer :: next
    end type

contains

    subroutine destroy(this)
        type(adaptor), intent(inout) :: this
        nullify(this%f)
        nullify(this%mf)
        nullify(this%pre)
    end subroutine

    real(c_double) function func_aux(n,x,grad,data) bind(c)
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(out), optional :: grad(n)
        type(c_ptr), value :: data
        
        type(adaptor), pointer :: fdata
        call c_f_pointer(data,fdata) ! fdata contains both callback and data
        
        func_aux = fdata%f%eval(n,x,grad)
    end function

    subroutine mfunc_aux(m,result,n,x,grad,data) bind(c)
        integer(c_int), intent(in), value :: m
        real(c_double), intent(out) :: result(m)
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(out), optional :: grad(m,n)
        type(c_ptr), value :: data
        
        type(adaptor), pointer :: fdata
        call c_f_pointer(data,fdata) ! fdata contains both callback and data

        call fdata%mf%eval(m,result,n,x,grad)
    end subroutine

    subroutine precond_aux(n,x,v,vpre,data) bind(c)
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(in) :: v(n)
        real(c_double), intent(out) :: vpre(n)
        type(c_ptr), value :: data
        
        type(adaptor), pointer :: fdata
        call c_f_pointer(data,fdata) ! fdata contains both callback and data

        call fdata%pre%eval(n,x,v,vpre)
    end subroutine

end module

module void_module

    use, intrinsic::iso_c_binding, only: c_int, c_double
    implicit none
    private

    public :: func_if


    type, abstract, public :: nlopt_void
    end type

    abstract interface
        real(c_double) function func_if(n,x,grad,data)
            import c_int, c_double, nlopt_void
            integer(c_int), intent(in), value :: n
            real(c_double), intent(in) :: x(n)
            real(c_double), intent(out), optional :: grad(n)
            class(nlopt_void) :: data
        end function
    end interface

end module

module nlopt

    use, intrinsic :: iso_c_binding
    use nlopt_enum, only: NLOPT_FAILURE, NLOPT_FORCED_STOP
    use nlopt_c_interface, func => nlopt_func, mfunc => nlopt_mfunc
    use nlopt_user_func_mod, only: nlopt_user_func, nlopt_user_mfunc, nlopt_user_precond
    use adaptor_mod, only: adaptor, func_aux, mfunc_aux, precond_aux
    use void_module, only: nlopt_void, func_if

    implicit none
    private

    public :: opt
    public :: nlopt_version, nlopt_version_major, &
              nlopt_version_minor, nlopt_version_bugfix
    public :: nlopt_srand, nlopt_srand_time, algorithm_name

    ! Expose abstract types in order for the user to extend them!
    public :: nlopt_user_func, nlopt_user_mfunc, nlopt_user_precond

    public :: nlopt_void

    type, private :: nlopt_void_handle
        class(nlopt_void), pointer :: data => null()
        procedure(func_if), pointer, nopass :: func => null()
        type(nlopt_void_handle), pointer :: next => null()
    ! contains
        ! final :: destroy_nlopt_void_handle
    end type

    type, private :: nlopt_void_handle_list
        integer :: num_nodes = 0
        type(nlopt_void_handle), pointer :: head => null()
        type(nlopt_void_handle), pointer :: tail => null()
    ! contains
        ! final :: destroy_nlopt_void_handle_list
    end type

    type :: opt
        private
        type(c_ptr) :: o = c_null_ptr ! the "nlopt_opt" object
        integer(c_int) :: last_result = NLOPT_FAILURE
        real(c_double) :: last_optf = huge(1._c_double)
        integer(c_int) :: forced_stop_reason = NLOPT_FORCED_STOP

        type(adaptor), pointer :: objective => null() ! keep handle to objective function on Fortran side

        type(nlopt_void_handle), pointer :: objective_handle => null()
        type(nlopt_void_handle_list), public :: eq_constraint_list
        type(nlopt_void_handle_list), public :: neq_constraint_list
        ! type(nlopt_void_handle_list), pointer :: eq_mconstraint_handles => null()
        ! type(nlopt_void_handle_list), pointer :: neq_mconstraint_handles => null()
    contains

        procedure, public :: optimize

        procedure, public :: last_optimize_result
        procedure, public :: last_optimum_value

        procedure, public :: get_algorithm
        procedure, public :: get_algorithm_name
        procedure, public :: get_dimension

        procedure, private :: set_min_objective_oo
        procedure, public :: set_min_objective_new
        generic, public :: set_min_objective => set_min_objective_oo, set_min_objective_new

        procedure, private :: set_max_objective_oo
        procedure, private :: set_max_objective_new
        generic, public :: set_max_objective => set_max_objective_oo, set_max_objective_new

        procedure, public :: remove_inequality_constraints
        procedure, private :: add_inequality_constraint_classic
        procedure, public :: add_inequality_constraint_new
        procedure, private :: add_inequality_mconstraint_classic
        procedure, private :: add_inequality_constraint_oo
        procedure, private :: add_inequality_mconstraint_oo
        generic, public :: add_inequality_constraint => add_inequality_constraint_classic, add_inequality_constraint_oo

        procedure, public :: remove_equality_constraints
        procedure, public :: add_equality_constraint_cptr
        procedure, public :: add_equality_mconstraint_cptr
        procedure, public :: add_equality_constraint_oo
        procedure, public :: add_equality_mconstraint_oo
        generic, public :: add_equality_constraint => add_equality_constraint_cptr, add_equality_constraint_oo

        procedure, private :: set_lower_bounds_array
        procedure, private :: set_lower_bounds_scalar
        generic, public :: set_lower_bounds => set_lower_bounds_array, set_lower_bounds_scalar
        procedure, public :: get_lower_bounds

        procedure, private :: set_upper_bounds_array
        procedure, private :: set_upper_bounds_scalar
        generic, public :: set_upper_bounds => set_upper_bounds_array, set_upper_bounds_scalar
        procedure, public :: get_upper_bounds
        
        ! stopping criteria
        procedure, public :: set_stopval
        procedure, public :: get_stopval
        procedure, public :: set_ftol_rel
        procedure, public :: get_ftol_rel
        procedure, public :: set_ftol_abs
        procedure, public :: get_ftol_abs
        procedure, public :: set_xtol_rel
        procedure, public :: get_xtol_rel
        procedure, public :: set_xtol_abs
        procedure, public :: get_xtol_abs

        procedure, public :: set_maxeval
        procedure, public :: get_maxeval

        procedure, public :: get_numevals

        procedure, public :: set_maxtime
        procedure, public :: get_maxtime

        procedure, public :: set_force_stop
        procedure, public :: get_force_stop
        procedure, public :: force_stop

        procedure, public :: get_errmsg

        !
        ! more algorithm-specific parameters
        !
        procedure, public :: set_local_optimizer

        procedure, public :: set_population
        procedure, public :: get_population

        procedure, public :: set_vector_storage
        procedure, public :: get_vector_storage

        procedure, public :: set_default_initial_step
        procedure, private :: set_initial_step_array
        procedure, private :: set_initial_step_scalar
        generic, public :: set_initial_step => set_initial_step_array, set_initial_step_scalar
        procedure, public :: get_initial_step

        ! Overload assignment
        procedure, private :: assign_opt
        generic, public :: assignment(=) => assign_opt

        ! Destructor
        final :: destroy
    end type opt

    ! Overwrite Constructors
    interface opt
        module procedure new_opt
        module procedure copy_opt
    end interface

    ! interface nlopt_void_handle
    !     module procedure void_handle_constructor
    ! end interface

contains

    ! function void_handle_constructor(f,data) result(node)
    !     procedure(func_if) :: f
    !     class(nlopt_void), intent(in), target :: data
    !     type(nlopt_void_handle), pointer :: node

    !     allocate(node)
    !     node%func => f
    !     node%data => data
    ! end function

    subroutine void_handle_list_append(this,node)
        type(nlopt_void_handle_list), intent(inout) :: this
        type(nlopt_void_handle), pointer :: node

        if (associated(this%tail)) then
            allocate(this%tail%next,source=node)
            this%tail => this%tail%next
        else
            allocate(this%head,source=node)
            this%tail => this%head
        end if
        this%num_nodes = this%num_nodes + 1
    end subroutine

    type(opt) function new_opt(a,n)
        integer(c_int), intent(in) :: a, n
        new_opt%o = nlopt_create(a,n)
        new_opt%last_result = NLOPT_FAILURE
        new_opt%last_optf = huge(new_opt%last_optf)
        new_opt%forced_stop_reason = NLOPT_FORCED_STOP
        if (.not. c_associated(new_opt%o)) stop "bad allocation"
    end function

    ! Copy constructor
    type(opt) function copy_opt(f)
        type(opt), intent(in) :: f
        copy_opt%o = nlopt_copy(f%o)
        copy_opt%last_result = f%last_result
        copy_opt%last_optf = f%last_optf
        copy_opt%forced_stop_reason = f%forced_stop_reason
        if (c_associated(f%o) .and. (.not. c_associated(copy_opt%o))) stop "bad allocation" 
    end function

    ! Assignment
    subroutine assign_opt(lhs,rhs)
        class(opt), intent(inout) :: lhs
        class(opt), intent(in) :: rhs
        if (c_associated(lhs%o,rhs%o)) return ! self-assignment
        call destroy(lhs) ! perform finalization
        lhs%o = nlopt_copy(rhs%o) ! nlopt_copy leads to a memory leak here?
        if (c_associated(rhs%o) .and. (.not. c_associated(lhs%o))) stop "bad allocation"
        lhs%last_result = rhs%last_result
        lhs%last_optf = rhs%last_optf
        lhs%forced_stop_reason = rhs%forced_stop_reason
        lhs%objective => rhs%objective ! transfer objective
        ! constraints
    end subroutine

    ! Finalizer/destructor
    subroutine destroy(this)
        type(opt) :: this

        ! Fortran handle to objective
        if (associated(this%objective)) then
            deallocate(this%objective)
            nullify(this%objective)
        end if

        call nlopt_destroy(this%o)
    end subroutine

    ! Do the optimization
    integer(c_int) function optimize(this,x,opt_f)
        class(opt), intent(inout) :: this
        real(c_double), intent(inout) :: x(nlopt_get_dimension(this%o))
        real(c_double), intent(inout) :: opt_f
        integer(c_int) :: ret
        
        this%forced_stop_reason = NLOPT_FORCED_STOP
        ret = nlopt_optimize(this%o,x,opt_f)
        this%last_result = ret
        this%last_optf = opt_f
        optimize = ret
    end function


    pure integer(c_int) function last_optimize_result(this)
        class(opt), intent(in) :: this
        last_optimize_result = this%last_result
    end function
    pure real(c_double) function last_optimum_value(this)
        class(opt), intent(in) :: this
        last_optimum_value = this%last_optf
    end function

    !
    ! Accessors
    !

    integer(c_int) function get_algorithm(this)
        class(opt), intent(in) :: this
        if (.not. c_associated(this%o)) stop "uninitialized type(nlopt_opt) instance"
        get_algorithm = nlopt_get_algorithm(this%o)
    end function
    function get_algorithm_name(this) result(name)
        class(opt), intent(in) :: this
        character(len=:,kind=c_char), allocatable :: name
        if (.not. c_associated(this%o)) stop "uninitialized type(nlopt_opt) instance"
        name = algorithm_name(this%get_algorithm())
    end function
    pure integer(c_int) function get_dimension(this)
        class(opt), intent(in) :: this
        ! if (.not. c_associated(this%o)) stop "uninitialized type(nlopt_opt) instance"
        get_dimension = nlopt_get_dimension(this%o)
    end function

    !
    ! Set the objective function
    !
    subroutine set_min_objective_new(this,f,f_data,ires)
        class(opt), intent(inout) :: this
        procedure(func_if) :: f
        class(nlopt_void), target :: f_data
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret

        type(c_ptr) :: c_handle
        type(c_funptr) :: c_fun_handle

        if (associated(this%objective_handle)) then
            nullify(this%objective_handle)
        end if
        allocate(this%objective_handle)
        this%objective_handle%func => f
        this%objective_handle%data => f_data

        c_handle = c_loc(this%objective_handle)
        c_fun_handle = c_funloc(nlopt_function_poly_c)

        ret = nlopt_set_min_objective(this%o,c_fun_handle,c_handle)
        if (present(ires)) ires = ret
    end subroutine
    subroutine set_max_objective_new(this,f,f_data,ires)
        class(opt), intent(inout) :: this
        procedure(func_if) :: f
        class(nlopt_void), target :: f_data
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret

        type(nlopt_void_handle), target :: p_handle
        type(c_ptr) :: c_handle
        type(c_funptr) :: c_fun_handle

        p_handle%data => f_data
        p_handle%func => f

        c_handle = c_loc(p_handle)
        c_fun_handle = c_funloc(nlopt_function_poly_c)

        ret = nlopt_set_max_objective(this%o,c_fun_handle,c_handle)
        if (present(ires)) ires = ret
    end subroutine
    real(c_double) function nlopt_function_poly_c(n,x,grad,func_data) bind(c)
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n)
        real(c_double), intent(out), optional :: grad(n)
        type(c_ptr), value :: func_data
        type(nlopt_void_handle), pointer :: p_handle
        call c_f_pointer(func_data,p_handle)
        nlopt_function_poly_c = p_handle%func(n,x,grad,p_handle%data)
    end function
    !
    ! Set the objective function (object-oriented way)
    !
    subroutine set_min_objective_oo(this,f,ires)
        class(opt), intent(inout) :: this
        class(nlopt_user_func), intent(in), target :: f
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret

        if (associated(this%objective)) nullify(this%objective)
        allocate(this%objective)
        this%objective%f => f

        ret = nlopt_set_min_objective(this%o,c_funloc(func_aux),c_loc(this%objective))
        if (present(ires)) ires = ret
    end subroutine
    subroutine set_max_objective_oo(this,f,ires)
        class(opt), intent(inout) :: this
        class(nlopt_user_func), intent(in), target :: f
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret

        if (associated(this%objective)) nullify(this%objective)
        allocate(this%objective)
        this%objective%f => f

        ret = nlopt_set_max_objective(this%o,c_funloc(func_aux),c_loc(this%objective))
        if (present(ires)) ires = ret
    end subroutine

    !
    ! Nonlinear constraints
    !
    subroutine remove_inequality_constraints(this,ires)
        class(opt), intent(inout) :: this
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_remove_inequality_constraints(this%o)
        if (present(ires)) ires = ret
    end subroutine
    subroutine add_inequality_constraint_classic(this,fc,fc_data,tol,ires)
        class(opt), intent(inout) :: this
        procedure(func) :: fc
        type(c_ptr), value :: fc_data
        real(c_double), optional :: tol
        integer(c_int), intent(out), optional :: ires
        real(c_double) :: tol_
        integer(c_int) :: ret

        tol_ = 0.0_c_double
        if (present(tol)) tol_ = tol

        ret = nlopt_add_inequality_constraint(this%o,c_funloc(fc),fc_data,tol_)
        if (present(ires)) ires = ret
    end subroutine
    subroutine add_inequality_constraint_new(this,fc,fc_data,tol,ires)
        class(opt), intent(inout) :: this
        procedure(func_if) :: fc
        class(nlopt_void), target :: fc_data
        real(c_double), optional :: tol
        integer(c_int), intent(out), optional :: ires
        real(c_double) :: tol_
        integer(c_int) :: ret

        type(nlopt_void_handle), pointer :: handle
        type(c_ptr) :: c_handle
        type(c_funptr) :: c_fun_handle

        tol_ = 0.0_c_double
        if (present(tol)) tol_ = tol

        allocate(handle)
        handle%func => fc
        handle%data => fc_data
        call void_handle_list_append(this%neq_constraint_list,handle)
        c_handle = c_loc(handle)
        c_fun_handle = c_funloc(nlopt_function_poly_c)

        ret = nlopt_add_inequality_constraint(this%o,c_fun_handle,c_handle,tol_)
        if (present(ires)) ires = ret
    end subroutine
    subroutine add_inequality_mconstraint_classic(this,m,fc,fc_data,tol,ires)
        class(opt), intent(inout) :: this
        integer(c_int), intent(in) :: m
        procedure(mfunc) :: fc
        type(c_ptr), value :: fc_data
        real(c_double), intent(in), optional :: tol(nlopt_get_dimension(this%o))
        integer(c_int), intent(out), optional :: ires
        real(c_double) :: tol_(nlopt_get_dimension(this%o))
        integer(c_int) :: ret

        tol_ = 0.0_c_double
        if (present(tol)) tol_ = tol
        
        ret = nlopt_add_inequality_mconstraint(this%o,m,c_funloc(fc),fc_data,tol_)
        if (present(ires)) ires = ret
    end subroutine
    subroutine add_inequality_constraint_oo(this,fc,tol,ires,cptr)
        class(opt), intent(inout) :: this
        class(nlopt_user_func), intent(in), target :: fc
        real(c_double), intent(in), optional :: tol
        integer(c_int), intent(out), optional :: ires
        real(c_double) :: tol_
        integer(c_int) :: ret
        type(adaptor), pointer :: fc_data => null()
        type(c_ptr), intent(out), optional :: cptr

        tol_ = 0.0_c_double
        if (present(tol)) tol_ = tol

        allocate(fc_data)
        fc_data%f => fc

        if (present(cptr)) cptr = c_loc(fc_data)

        ret = nlopt_add_inequality_constraint(this%o,c_funloc(func_aux),c_loc(fc_data),tol_)
        if (present(ires)) ires = ret
    end subroutine
    subroutine add_inequality_mconstraint_oo(this,m,fc,tol,ires)
        class(opt), intent(inout) :: this
        integer(c_int), intent(in) :: m
        class(nlopt_user_mfunc), intent(in), target :: fc
        real(c_double), intent(in), optional :: tol(nlopt_get_dimension(this%o))
        integer(c_int), intent(out), optional :: ires
        real(c_double) :: tol_(nlopt_get_dimension(this%o))
        integer(c_int) :: ret
        type(adaptor), pointer :: fc_data => null()
        tol_ = 0.0_c_double
        if (present(tol)) tol_ = tol
        
        allocate(fc_data)
        fc_data%mf => fc
        ret = nlopt_add_inequality_mconstraint(this%o,m,c_funloc(mfunc_aux),c_loc(fc_data),tol_)
        if (present(ires)) ires = ret
    end subroutine

    subroutine remove_equality_constraints(this,ires)
        class(opt), intent(inout) :: this
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_remove_inequality_constraints(this%o)
        if (present(ires)) ires = ret
    end subroutine
    subroutine add_equality_constraint_cptr(this,h,h_data,tol,ires)
        class(opt), intent(inout) :: this
        procedure(func) :: h
        type(c_ptr), value :: h_data
        real(c_double), optional :: tol
        integer(c_int), intent(out), optional :: ires
        real(c_double) :: tol_
        integer(c_int) :: ret

        tol_ = 0.0_c_double
        if (present(tol)) tol_ = tol

        ret = nlopt_add_equality_constraint(this%o,c_funloc(h),h_data,tol_)
        if (present(ires)) ires = ret
    end subroutine
    subroutine add_equality_mconstraint_cptr(this,m,h,h_data,tol,ires)
        class(opt), intent(inout) :: this
        integer(c_int), intent(in) :: m
        procedure(mfunc) :: h
        type(c_ptr), value :: h_data
        real(c_double), optional :: tol(nlopt_get_dimension(this%o))
        integer(c_int), intent(out), optional :: ires
        real(c_double) :: tol_(nlopt_get_dimension(this%o))
        integer(c_int) :: ret

        tol_ = 0
        if (present(tol)) tol_ = tol

        ret = nlopt_add_equality_mconstraint(this%o,m,c_funloc(h),h_data,tol_)
        if (present(ires)) ires = ret
    end subroutine
    subroutine add_equality_constraint_oo(this,h,tol)
        class(opt), intent(inout) :: this
        class(nlopt_user_func), intent(in), target :: h
        real(c_double), intent(in), optional :: tol
        real(c_double) :: tol_
        integer(c_int) :: ret
        type(adaptor), pointer :: h_data => null()

        tol_ = 0.0_c_double
        if (present(tol)) tol_ = tol

        allocate(h_data)
        h_data%f => h
        ret = nlopt_add_equality_constraint(this%o,c_funloc(func_aux),c_loc(h_data),tol_)
    end subroutine
    subroutine add_equality_mconstraint_oo(this,m,h,tol)
        class(opt), intent(inout) :: this
        integer(c_int), intent(in) :: m
        class(nlopt_user_mfunc), intent(in), target :: h
        real(c_double), intent(in), optional :: tol(nlopt_get_dimension(this%o))
        real(c_double) :: tol_(nlopt_get_dimension(this%o))
        integer(c_int) :: ret
        type(adaptor), pointer :: h_data => null()
        tol_ = 0.0_c_double
        if (present(tol)) tol_ = tol
        
        allocate(h_data)
        h_data%mf => h
        ret = nlopt_add_equality_mconstraint(this%o,m,c_funloc(mfunc_aux),c_loc(h_data),tol_)
    end subroutine

    subroutine set_lower_bounds_array(this,lb,ires)
        class(opt), intent(inout) :: this
        real(c_double), intent(in) :: lb(nlopt_get_dimension(this%o))
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_set_lower_bounds(this%o,lb)
        if (present(ires)) ires = ret
    end subroutine
    subroutine set_lower_bounds_scalar(this,lb,ires)
        class(opt), intent(inout) :: this
        real(c_double), intent(in) :: lb
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_set_lower_bounds1(this%o,lb)
        if (present(ires)) ires = ret
    end subroutine
    subroutine get_lower_bounds(this,lb,ires)
        class(opt), intent(in) :: this
        real(c_double), intent(out) :: lb(nlopt_get_dimension(this%o))
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_get_lower_bounds(this%o,lb)
        if (present(ires)) ires = ret
    end subroutine
    subroutine set_upper_bounds_array(this,ub,ires)
        class(opt), intent(inout) :: this
        real(c_double), intent(in) :: ub(nlopt_get_dimension(this%o))
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_set_upper_bounds(this%o,ub)
        if (present(ires)) ires = ret
    end subroutine
    subroutine set_upper_bounds_scalar(this,ub,ires)
        class(opt), intent(inout) :: this
        real(c_double), intent(in) :: ub
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_set_upper_bounds1(this%o,ub)
        if (present(ires)) ires = ret
    end subroutine
    subroutine get_upper_bounds(this,ub,ires)
        class(opt), intent(in) :: this
        real(c_double), intent(out) :: ub(nlopt_get_dimension(this%o))
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_get_upper_bounds(this%o,ub)
        if (present(ires)) ires = ret
    end subroutine


    real(c_double) function get_stopval(this)
        class(opt), intent(in) :: this
        get_stopval = nlopt_get_stopval(this%o)
    end function
    subroutine set_stopval(this,stopval,ires)
        class(opt), intent(inout) :: this
        real(c_double), intent(in) :: stopval
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_set_stopval(this%o,stopval)
        if (present(ires)) ires = ret
    end subroutine

    real(c_double) function get_ftol_rel(this)
        class(opt), intent(in) :: this
        get_ftol_rel = nlopt_get_ftol_rel(this%o)
    end function
    subroutine set_ftol_rel(this,tol,ires)
        class(opt), intent(inout) :: this
        real(c_double), intent(in) :: tol
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_set_ftol_rel(this%o,tol)
        if (present(ires)) ires = ret
    end subroutine

    real(c_double) function get_ftol_abs(this)
        class(opt), intent(in) :: this
        get_ftol_abs = nlopt_get_ftol_abs(this%o)
    end function
    subroutine set_ftol_abs(this,tol,ires)
        class(opt), intent(inout) :: this
        real(c_double), intent(in) :: tol
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_set_ftol_abs(this%o,tol)
        if (present(ires)) ires = ret
    end subroutine

    real(c_double) function get_xtol_rel(this)
        class(opt), intent(in) :: this
        get_xtol_rel = nlopt_get_xtol_rel(this%o)
    end function
    subroutine set_xtol_rel(this,tol,ires)
        class(opt), intent(inout) :: this
        real(c_double), intent(in) :: tol
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_set_xtol_rel(this%o,tol)
        if (present(ires)) ires = ret
    end subroutine
    subroutine set_xtol_abs(this,tol,ires)
        class(opt), intent(inout) :: this
        real(c_double), intent(in) :: tol(nlopt_get_dimension(this%o))
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_set_xtol_abs(this%o,tol)
        if (present(ires)) ires = ret
    end subroutine
    subroutine get_xtol_abs(this,tol,ires)
        class(opt), intent(in) :: this
        real(c_double), intent(out) :: tol(nlopt_get_dimension(this%o))
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_get_xtol_abs(this%o,tol)
        if (present(ires)) ires = ret
    end subroutine


    integer(c_int) function get_maxeval(this)
        class(opt), intent(in) :: this
        ! if (.not. associated(this%o)) 
        get_maxeval = nlopt_get_maxeval(this%o)
    end function
    subroutine set_maxeval(this,maxeval,ires)
        class(opt), intent(inout) :: this
        integer(c_int), intent(in) :: maxeval
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_set_maxeval(this%o,maxeval)
        if (present(ires)) ires = ret
    end subroutine

    integer(c_int) function get_numevals(this)
        class(opt), intent(in) :: this
        get_numevals = nlopt_get_numevals(this%o)
    end function


    real(c_double) function get_maxtime(this)
        class(opt), intent(in) :: this
        get_maxtime = nlopt_get_maxtime(this%o)
    end function
    subroutine set_maxtime(this,maxtime,ires)
        class(opt), intent(inout) :: this
        real(c_double), intent(in) :: maxtime
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_set_maxtime(this%o,maxtime)
        if (present(ires)) ires = ret
    end subroutine

    integer(c_int) function get_force_stop(this)
        class(opt), intent(in) :: this
        ! if (.not. associated(this%o)) 
        get_force_stop = nlopt_get_force_stop(this%o)
    end function
    subroutine set_force_stop(this,ival,ires)
        class(opt), intent(inout) :: this
        integer(c_int), intent(in) :: ival
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_set_force_stop(this%o,ival)
        if (present(ires)) ires = ret
    end subroutine
    subroutine force_stop(this)
        class(opt), intent(inout) :: this
        call this%set_force_stop(1_c_int)
    end subroutine


    function get_errmsg(this) result(errmsg)
        class(opt), intent(in) :: this
        character(len=:,kind=c_char), allocatable :: errmsg
        type(c_ptr) :: c_string
        character(len=1000,kind=c_char), pointer :: f_string
        
        c_string = nlopt_get_errmsg(this%o)
        if (.not. c_associated(c_string)) then
            errmsg = ""
        else
            call c_f_pointer(c_string,f_string)
            errmsg = f_string(1:index(f_string,c_null_char))
        end if
    end function

    ! potential memory leak!?
    subroutine set_local_optimizer(this,lo,ires)
        class(opt), intent(inout) :: this
        class(opt), intent(in) :: lo
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_set_local_optimizer(this%o,lo%o)
        if (present(ires)) ires = ret
    end subroutine

    integer(c_int) function get_population(this)
        class(opt), intent(in) :: this
        get_population = nlopt_get_population(this%o)
    end function
    subroutine set_population(this,pop,ires)
        class(opt), intent(inout) :: this
        integer(c_int), intent(in) :: pop
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_set_population(this%o,pop)
        if (present(ires)) ires = ret
    end subroutine

    integer(c_int) function get_vector_storage(this)
        class(opt), intent(in) :: this
        get_vector_storage = nlopt_get_vector_storage(this%o)
    end function
    subroutine set_vector_storage(this,dim,ires)
        class(opt), intent(inout) :: this
        integer(c_int), intent(in) :: dim
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_set_vector_storage(this%o,dim)
        if (present(ires)) ires = ret
    end subroutine

    subroutine set_default_initial_step(this,x,ires)
        class(opt), intent(inout) :: this
        real(c_double), intent(in) :: x(nlopt_get_dimension(this%o))
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_set_default_initial_step(this%o,x)
        if (present(ires)) ires = ret
    end subroutine

    subroutine set_initial_step_array(this,dx,ires)
        class(opt), intent(inout) :: this
        real(c_double), intent(in) :: dx(nlopt_get_dimension(this%o))
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_set_initial_step(this%o,dx)
        if (present(ires)) ires = ret
    end subroutine
    subroutine set_initial_step_scalar(this,dx,ires)
        class(opt), intent(inout) :: this
        real(c_double), intent(in) :: dx
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_set_initial_step1(this%o,dx)
        if (present(ires)) ires = ret
    end subroutine
    subroutine get_initial_step(this,x,dx,ires)
        class(opt), intent(in) :: this    
        real(c_double), intent(in) :: x(nlopt_get_dimension(this%o))
        real(c_double), intent(out) :: dx(nlopt_get_dimension(this%o))
        integer(c_int), intent(out), optional :: ires
        integer(c_int) :: ret
        ret = nlopt_get_initial_step(this%o,x,dx)
        if (present(ires)) ires = ret
    end subroutine


    integer(c_int) function nlopt_version_major()
        integer(c_int) :: minor, bugfix
        call nlopt_version(nlopt_version_major,minor,bugfix)
    end function
    integer(c_int) function nlopt_version_minor()
        integer(c_int) :: major, bugfix
        call nlopt_version(major,nlopt_version_minor,bugfix)
    end function
    integer(c_int) function nlopt_version_bugfix()
        integer(c_int) :: major, minor
        call nlopt_version(major,minor,nlopt_version_bugfix)
    end function

    function algorithm_name(a) result(name)
        integer(c_int), intent(in) :: a
        character(len=:,kind=c_char), allocatable :: name
        character(len=256,kind=c_char), pointer :: f_string
        call c_f_pointer(nlopt_algorithm_name(a),f_string)
        name = f_string(1:index(f_string,c_null_char))
    end function

end module nlopt


