# NLopt Modern Fortran Reference

NLopt is written in C and the C NLopt programming interface (API), as described in the [NLopt Reference](NLopt_Reference.md)., is interoperable with Modern Fortran by using appropriate interfaces.

However, we also provide another Modern Fortran (from here on just Fortran) module in file nlopt.f90, that wraps a more natural Fortran interface around the NLopt API, which may be more convenient for Fortran programmers used to using derived types.

The main distinctions of the Modern Fortran API are:

* Use of an `nlopt` module.
* Use of a Fortran derived-type `opt`, with constructors, destructors, and member functions.

The main purpose of this section is to document the syntax and unique features of the Modern Fortran API; for more details on the underlying features, please refer to the C documention in the [NLopt Reference](NLopt_Reference.md).

## Compiling and linking your program to NLopt

To obtain the definitions of the NLopt constants in Fortran, your program should include the following line:
```Fortran
use nlopt_enums
```
An NLopt program in Fortran should import the `nlopt` module:

```Fortran
use nlopt
```

On Unix you would normally link your program exactly as for the C API, with a command like:

*`compiler`*` `*`...source/object` `files...`*` -lnlopt -lm -o myprogram`

where *compiler* is `gfortran` or whatever is appropriate for your machine/system.

A distinction in Fortran is that the path to the .mod files should also be passed.

## The `nlopt_opt` derived-type

The Modern Fortran NLopt API revolves around an object of the derived-type `type(nlopt_opt)`. Internally, this object stores a pointer to the corresponding C object of type `nlopt_opt`. Via methods of this object (type-bound procedures) all of the parameters of the optimization are specified (dimensions, algorithm, stopping criteria, constraints, objective function, etcetera), and then one finally calls `optimize` in order to perform the optimization. The object should normally be created via the constructor:

```Fortran
type(nlopt_opt) :: opt

opt = nlopt_opt(algorithm,n)
```

given an integer `algorithm` (see [NLopt Algorithms]() for possible values, defined in the `nlopt_enums` module) and the dimensionality of the problem (`n`, the number of optimization parameters). Just as in C, the algorithm is specified by constants of the form `NLOPT_MMA`, `NLOPT_COBYLA`, etcetera.

There are also a copy constructor and an assignment operator, both of which make a copy of a given object (equivalent to `nlopt_copy` in the C API):

```Fortran

copt = nlopt_opt(opt) ! copy
copt = opt !
```

When you are finished with your object, any associated storage will deallocated automatically by a [final subroutine]() once the `type(nlopt_opt)` object ceases to exist. This subroutine internally makes a call to the C function `nlopt_destroy`.

The algorithm and dimension parameters of the object are immutable (cannot be changed without constructing a new object), but you can query them for a given object by the subroutines:

```Fortran
call opt%get_algorithm(algorithm)
call opt%get_dimension(n)
```

You can get a Fortran-style string description of the algorithm via the function:
```Fortran
opt%get_algorithm_name()
```

(These accessor methods, along with the other methods below, will throw an exception if you use them on an object initialized with the default no-argument constructor, i.e. if you didn't specify an algorithm or dimensionality yet.)

## Objective function

The objective function is specified by calling one of the methods:

```Fortran
    call opt%set_min_objective(f,f_data [,ires])
    call opt%set_max_objective(f,f_data [,ires])
```

depending on whether one wishes to minimize or maximize the objective function `f`, respectively. The function `f` should be of the form:

```Fortran
function f(n,x,grad,f_data)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n)
    real(dp), intent(out), optional :: grad(n)
    type(c_ptr), value :: f_data
end function
```

The return value should be the value of the function at the point `x`, where `x` is a vector of length `n` of the optimization parameters (the same as the dimension passed to the constructor).

In addition, if the argument `grad` is present, then `grad` is a vector of length `n` which should (upon return) be set to the gradient of the function with respect to the optimization parameters at `x`. That is, `grad(i)` should upon return contain the partial derivative $\partial f / \partial x_i$, for $0 \leq i < n$, if `grad` is present. Not all of the optimization algorithms (below) use the gradient information: for algorithms listed as "derivative-free," the `grad` argument will always be empty and need never be computed. (For algorithms that do use gradient information, however, `grad` may still be empty for some calls.)

## Bound constraints

The [bound constraints](NLopt_Reference#Bound_constraints.md) can be specified by calling the methods:

```Fortran
real(dp) :: lb(n), ub(n)

call opt%set_lower_bounds(lb [,ires])
call opt%set_upper_bounds(ub [,ires])
```
where `lb` and `ub` are arrays of length `n` (the same as the dimension passed to the constructor `opt` constructor), and `ires` is an optional integer [return value](NLopt_Reference#Return_values.md) (positive on success). 

For convenience, these type-bound procedures are overloaded with subroutines that take a single number as arguments, in order to set the lower/upper bounds for all optimization parameters to a constant value:

```Fortran
real(dp) :: lb1, ub1

call opt%set_lower_bounds(lb1 [,ires])
call opt%set_upper_bounds(ub1 [,ires])
```

To retrieve the values of the lower/upper bounds, you can call one of:

```Fortran
call opt%get_lower_bounds(lb [,ires])
call opt%get_upper_bounds(ub [,ires])
```

To specify an unbounded dimension, you can use $\pm$`huge(lb)` in Fortran to specify $\pm \infty$, where `huge` is a Fortran intrinsic function.

## Nonlinear constraints

Just as for [nonlinear constraints in C](NLopt_Reference#Nonlinear_constraints.md), you can specify nonlinear inequality and equality constraints by the methods:

```Fortran

call opt%add_inequality_constraint(fc,tol [,ires])
call opt%add_equality_constraint(h,tol [,ires])
```

where the arguments `fc` and `h` have the same form as the objective function above. The `tol` arguments specify a tolerance in judging the feasibility for the purposes of stopping the optimization, as in C, and `ires` is an optional `integer` [return value](NLopt_Reference#Return_values.md) (positive on sucess).

To remove all of the inequality and/or equality constraints from a given problem, you can call the following methods:

```Fortran
call opt%remove_inequality_constraints([ires])
call opt%remove_equality_constraints([ires])
```

### Vector-valued constraints

Just as for [nonlinear constraints in C](NLopt_Reference#Vector-valued_constraints.md), you can specify nonlinear inequality and equality constraints by the methods:

```Fortran
call opt%add_inequality_mconstraint(m,c,tol [,ires])
call opt%add_equality_mconstraint(m,c,tol [,ires])
```

Here, `m` is the dimensionality of the constraint result and `tol` is a length-`m` array of the tolerances in each constraint dimension. The constraint subroutine `c` must be of the form:

```Fortran
subroutine c(m,result,n,x,grad,f_data)
    integer, intent(in) :: m
    real(dp), intent(out) :: result(m)
    integer, intent(in) :: n
    real(dp), intent(in) :: x(n)

end subroutine
```

(You can add multiple vector-valued constraints and/or scalar constraints in the same problem.)

Stopping criteria
-----------------

As explained in the [C API Reference](NLopt_Reference#Stopping_criteria.md) and the [Introduction](NLopt_Introduction#Termination_conditions.md)), you have multiple options for different stopping criteria that you can specify. (Unspecified stopping criteria are disabled; i.e., they have innocuous defaults.)

For each stopping criteria, there are (at least) two method: a `set` method to specify the stopping criterion, and a `get` method to retrieve the current value for that criterion. The meanings of each criterion are exactly the same as in the C API.

```Fortran
call opt%set_stopval(stopval [,ires])
call opt%get_stopval(stopval)
```
Stop when an objective value of at least `stopval` is found.

```Fortran
call opt%set_ftol_rel(tol [,ires])
call opt%get_ftol_rel(tol)
```
Set relative tolerance on function value.

```Fortran
call opt%set_ftol_abs(tol [,ires])
call opt%get_ftol_abs(tol)
```
Set absolute tolerance on function value.

```Fortran
call opt%set_xtol_rel(tol [,ires])
call opt%get_xtol_rel(tol)
```
Set relative tolerance on optimization parameters.

```Fortran
call opt%set_xtol_abs(tol [,ires])
call opt%set_xtol_abs(tol1 [,ires])
call opt%get_xtol_abs(tol [,ires])
```
Set absolute tolerances on optimization parameters. The `tol` input must be an array of length `n` (the dimension specified in the `opt` constructor). Alternatively, you can pass a single number `tol1` in order to set the same tolerance for all optimization parameters. The type-bound procedure `get_xtol_abs` returns tolerances in the array `tol`.

```Fortran
call opt%set_maxeval(maxeval [,ires])
call opt%get_maxeval(maxeval)
```
Stop when the number of function evaluations exceeds the integer `maxeval` (zero or negative for no limit).

```Fortran
call opt%get_numevals(nevals)
```
Request the number of evaluations.

```Fortran
call opt%set_maxtime(maxtime [,ires])
call opt%get_maxtime(maxtime)
```
Stop when the optimization time (in seconds) exceeds `maxtime` (double precision). (Zero or negative for no limit.)

### Forced termination

In certain cases, the caller may wish to force the optimization to halt, for some reason unknown to NLopt. For example, if the user presses Ctrl-C, or there is an error of some sort in the objective function. In this case, it is possible to tell NLopt to halt the optimization gracefully, returning the best point found so far, by calling the following subroutine from within your objective or constraint functions (exactly analogous to the corresponding C routines):

```Fortran
call opt%force_stop([ires])
```
`ires` is an optional integer return value (positive on success). More generally, you can set and retrieve a force-stop integer code `ival`, where a nonzero value indicates a forced stop.

```Fortran
call opt%set_force_stop(ival [,ires])
call opt%get_force_stop(ival [,ires])
```
The force-stop value is reset to zero at the beginning of `opt%optimize()`. Passing `ival=0` to `opt%set_force_stop()` tells NLopt not to force a halt.

Performing the optimization
---------------------------

Once all the desired optimization parameters have been specified in a given object `myopt`, you can perform the optimization by calling:

```Fortran
real(dp) :: x(n), optf
call opt%optimize(x, optf [,ires])
```
On input, `x` is an array of length `n` (the dimension of the problem from the `opt` constructor) giving an initial guess for the optimization parameters. Upon succesful return, `x` contains the optimized values of the optimization parameters and `optf` contains the corresponding value of the objective function. `ires` is an optional integer return value (positive on sucess).

You can also call the following methods to retrieve the `optf` value from the last `optimize` call, and the integer return value (including negative/failure return values) from the last `optimize` call:

```Fortran
real(dp) :: optf
integer :: ires

optf = opt%last_optimum_value()
ires = opt%last_optimize_result()
```

### Return values

The possible return values are the same as the [return values in the C API](NLopt_Reference#Return_values.md), with the corresponding integer constants defined in the `nlopt_enum` module.


Local/subsidiary optimization algorithm
---------------------------------------

Some of the algorithms, especially MLSL and AUGLAG, use a different optimization algorithm as a subroutine, typically for local optimization. You can change the local search algorithm and its tolerances by calling:

```Fortran
call opt%set_local_optimizer(local_opt [,ires])
```

Here, `local_opt` is another object of `type(opt_type)` whose parameters are used to determine the local search algorithm, its stopping criteria, and other algorithm parameters. (However, the objective function, bounds, and nonlinear-constraint parameters of `local_opt` are ignored.) The dimension `n` of `local_opt` must match that of `opt`. `ires` is an optional integer return value (positive on sucess).

This function makes a copy of the `local_opt` object, so you can freely destroy your original `local_opt` afterwards without affecting `opt`.

Initial step size
-----------------

Just as in the C API, you can [get and set the initial step sizes](NLopt_Reference#Initial_step_size.md) for derivative-free optimization algorithms. The Modern Fortran equivalents of the C functions are the following methods:

```Fortran
real(dp) :: x(n), dx(n), dx1
call opt%set_initial_step(dx [,ires])
call opt%set_initial_step(dx1 [,ires])
call opt%get_initial_step(x,dx [,ires])
```

Here, `dx` is an array of the (nonzero) initial steps for each dimension. For convenience, you can also pass a single number `dx1` to `opt%set_initial_step` if you wish to use the same initial steps for all dimensions. A call to `opt%get_initial_step` sets `dx` to the initial step that will be used for a starting guess of `x` in a call to `opt%optimize(x,optf [,ires])`. `ires` is an optional integer return value (positive on success).

Stochastic population
---------------------

Just as in the C API, you can [get and set the initial population](NLopt_Reference#Stochastic_population.md) for stochastic optimization algorithms, by the methods:

```Fortran
call opt%set_population(ipop [,ires])
call opt%get_population(ipop)
```
where `ipop` is an integer and `ires` is an optional integer return value (positive on success). (An `ipop` of zero implies that the heuristic default will be used.)

Pseudorandom numbers
--------------------

For stochastic optimization algorithms, we use pseudorandom numbers generated by the [Mersenne Twister](https://en.wikipedia.org/wiki/Mersenne_twister) algorithm, based on code from Makoto Matsumoto. By default, the [seed](https://en.wikipedia.org/wiki/Random_seed) for the random numbers is generated from the system time, so that you will get a different sequence of pseudorandom numbers each time you run your program. If you want to use a "deterministic" sequence of pseudorandom numbers, i.e. the same sequence from run to run, you can set the seed by calling:

```Fortran
call nlopt_srand(iseed)
```
where `iseed` is an integer. To reset the seed based on the system time, you can call:

```Fortran
call nlopt_srand_time()
```

(Normally, you don't need to call this as it is called automatically. However, it might be useful if you want to "re-randomize" the pseudorandom numbers after calling `nlopt_srand` to set a deterministic seed.)

Vector storage for limited-memory quasi-Newton algorithms
---------------------------------------------------------

Just as in the C API, you can get and set the [number *M* of stored vectors](NLopt_Reference#Vector_storage_for_limited-memory_quasi-Newton_algorithms.md) for limited-memory quasi-Newton algorithms, via the methods:

```Fortran
call opt%set_vector_storage(M [,ires])
call opt%get_vector_storage(M)
```
(The default is *M*=0, in which case NLopt uses a heuristic nonzero value.)

Version number
--------------

To determine the version number of NLopt at runtime, you can call:

```Fortran
call nlopt_version(major,minor,bugfix)
```

where the three arguments are integers. For example, NLopt version 3.1.4 would return `major=3`, `minor=1`, and `bugfix=4`. You can also retrieve these three values individually by calling:

```Fortran
major = nlopt_version_major()
minor = nlopt_version_minor()
bugfix = nlopt_version_bugfix()
```


[Category:NLopt](index.md)
