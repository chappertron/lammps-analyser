# fix aveforce command

## Syntax

``` LAMMPS
fix ID group-ID aveforce fx fy fz keyword value ...
```

-   ID, group-ID are documented in [fix](fix) command

-   aveforce = style name of this fix command

-   fx,fy,fz = force component values (force units)

        any of fx,fy,fz can be a variable (see below)

-   zero or more keyword/value pairs may be appended to args

-   keyword = *region*

        *region* value = region-ID
          region-ID = ID of region atoms must be in to have added force

## Examples

``` LAMMPS
fix pressdown topwall aveforce 0.0 -1.0 0.0
fix 2 bottomwall aveforce NULL -1.0 0.0 region top
fix 2 bottomwall aveforce NULL -1.0 v_oscillate region top
```

## Description

Apply an additional external force to a group of atoms in such a way
that every atom experiences the same force. This is useful for pushing
on wall or boundary atoms so that the structure of the wall does not
change over time.

The existing force is averaged for the group of atoms, component by
component. The actual force on each atom is then set to the average
value plus the component specified in this command. This means each atom
in the group receives the same force.

Any of the *fx*, *fy*, or *fz* values can be specified as `NULL`, which
means the force in that dimension is not changed. Note that this is not
the same as specifying a 0.0 value, since that sets all forces to the
same average value without adding in any additional force.

Any of the three quantities defining the force components, namely *fx*,
*fy*, and *fz*, can be specified as an equal-style [variable](variable).
If the value is a variable, it should be specified as v_name, where name
is the variable name. In this case, the variable will be evaluated each
timestep, and its value used to determine the average force.

Equal-style variables can specify formulas with various mathematical
functions, and include [thermo_style](thermo_style) command keywords for
the simulation box parameters and timestep and elapsed time. Thus it is
easy to specify a time-dependent average force.

If the *region* keyword is used, the atom must also be in the specified
geometric [region](region) in order to have force added to it.

------------------------------------------------------------------------

Styles with a *gpu*, *intel*, *kk*, *omp*, or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the [Accelerator packages](Speed_packages)
page. The accelerated styles take the same arguments and should produce
the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, INTEL, KOKKOS, OPENMP, and
OPT packages, respectively. They are only enabled if LAMMPS was built
with those packages. See the [Build package](Build_package) page for
more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the [-suffix command-line
switch](Run_options) when you invoke LAMMPS, or you can use the
[suffix](suffix) command in your input script.

See the [Accelerator packages](Speed_packages) page for more
instructions on how to use the accelerated styles effectively.

------------------------------------------------------------------------

## Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to [binary restart
files](restart).

The [fix_modify](fix_modify) *respa* option is supported by this fix.
This allows to set at which level of the [r-RESPA](run_style) integrator
the fix is adding its forces. Default is the outermost level.

This fix computes a global three-vector of forces, which can be accessed
by various [output commands](Howto_output). This is the total force on
the group of atoms before the forces on individual atoms are changed by
the fix. The vector values calculated by this fix are \"extensive\".

No parameter of this fix can be used with the *start/stop* keywords of
the [run](run) command.

The forces due to this fix are imposed during an energy minimization,
invoked by the [minimize](minimize) command. You should not specify
force components with a variable that has time-dependence for use with a
minimizer, since the minimizer increments the timestep as the iteration
count during the minimization.

## Restrictions

> none

## Related commands

[fix setforce](fix_setforce), [fix addforce](fix_addforce)

## Default

none
