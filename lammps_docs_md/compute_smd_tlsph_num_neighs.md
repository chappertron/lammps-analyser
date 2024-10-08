# compute smd/tlsph/num/neighs command

## Syntax

    compute ID group-ID smd/tlsph/num/neighs

-   ID, group-ID are documented in [compute](compute) command
-   smd/tlsph/num/neighs = style name of this compute command

## Examples

``` LAMMPS
compute 1 all smd/tlsph/num/neighs
```

## Description

Define a computation that calculates the number of particles inside of
the smoothing kernel radius for particles interacting via the
Total-Lagrangian SPH pair style.

See [this PDF guide](PDF/MACHDYN_LAMMPS_userguide.pdf)\_ to using Smooth
Mach Dynamics in LAMMPS.

## Output info

This compute calculates a per-particle vector, which can be accessed by
any command that uses per-particle values from a compute as input. See
the [Howto output](Howto_output) page for an overview of LAMMPS output
options.

The per-particle values are dimensionless. See [units](units).

## Restrictions

This compute is part of the MACHDYN package. It is only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

This quantity will be computed only for particles which interact with
the Total-Lagrangian pair style.

## Related commands

[smd/ulsph/num/neighs](compute_smd_ulsph_num_neighs)

## Default

none
