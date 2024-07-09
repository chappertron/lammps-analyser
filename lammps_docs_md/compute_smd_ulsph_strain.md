# compute smd/ulsph/strain command

## Syntax

``` LAMMPS
compute ID group-ID smd/ulsph/strain
```

-   ID, group-ID are documented in [compute](compute) command
-   smd/ulsph/strain = style name of this compute command

## Examples

``` LAMMPS
compute 1 all smd/ulsph/strain
```

## Description

Define a computation that outputs the logarithmic strain tensor. for
particles interacting via the updated Lagrangian SPH pair style.

See [this PDF guide](PDF/MACHDYN_LAMMPS_userguide.pdf)\_ to using Smooth
Mach Dynamics in LAMMPS.

## Output info

This compute calculates a per-particle tensor, which can be accessed by
any command that uses per-particle values from a compute as input. See
the [Howto output](Howto_output) page for an overview of LAMMPS output
options.

The per-particle vector has 6 entries, corresponding to the xx, yy, zz,
xy, xz, yz components of the symmetric strain rate tensor.

The per-particle tensor values will be given dimensionless, see
[units](units).

## Restrictions

This compute is part of the MACHDYN package. It is only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info. This compute can only be
used for particles which interact with the updated Lagrangian SPH pair
style.

## Related commands

[compute smd/tlsph/strain](compute_smd_tlsph_strain)

## Default

none
