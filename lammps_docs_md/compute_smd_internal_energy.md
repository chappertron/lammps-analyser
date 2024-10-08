# compute smd/internal/energy command

## Syntax

    compute ID group-ID smd/internal/energy

-   ID, group-ID are documented in [compute](compute) command
-   smd/smd/internal/energy = style name of this compute command

## Examples

``` LAMMPS
compute 1 all smd/internal/energy
```

## Description

Define a computation which outputs the per-particle enthalpy, i.e., the
sum of potential energy and heat.

See [this PDF guide](PDF/MACHDYN_LAMMPS_userguide.pdf)\_ to use Smooth
Mach Dynamics in LAMMPS.

## Output Info

This compute calculates a per-particle vector, which can be accessed by
any command that uses per-particle values from a compute as input. See
the [Howto output](Howto_output) page for an overview of LAMMPS output
options.

The per-particle vector values will be given in [units](units) of
energy.

## Restrictions

This compute is part of the MACHDYN package. It is only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info. This compute can only be
used for particles which interact via the updated Lagrangian or total
Lagrangian SPH pair styles.

## Related commands

none

## Default
