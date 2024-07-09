# compute ave/sphere/atom command

Accelerator Variants: *ave/sphere/atom/kk*

## Syntax

``` LAMMPS
compute ID group-ID ave/sphere/atom keyword values ...
```

-   ID, group-ID are documented in [compute](compute) command

-   ave/sphere/atom = style name of this compute command

-   one or more keyword/value pairs may be appended

        keyword = *cutoff*
          *cutoff* value = distance cutoff

## Examples

``` LAMMPS
compute 1 all ave/sphere/atom

compute 1 all ave/sphere/atom cutoff 5.0
comm_modify cutoff 5.0
```

## Description

::: versionadded
7Jan2022
:::

Define a computation that calculates the local mass density and
temperature for each atom based on its neighbors inside a spherical
cutoff. If an atom has $M$ neighbors, then its local mass density is
calculated as the sum of its mass and its $M$ neighbor masses, divided
by the volume of the cutoff sphere (or circle in 2d). The local
temperature of the atom is calculated as the temperature of the
collection of $M+1$ atoms, after subtracting the center-of-mass velocity
of the $M+1$ atoms from each of the $M+1$ atom\'s velocities. This is
effectively the thermal velocity of the neighborhood of the central
atom, similar to [compute temp/com](compute_temp_com).

The optional keyword *cutoff* defines the distance cutoff used when
searching for neighbors. The default value is the cutoff specified by
the pair style. If no pair style is defined, then a cutoff must be
defined using this keyword. If the specified cutoff is larger than that
of the pair_style plus neighbor skin (or no pair style is defined), the
*comm_modify cutoff* option must also be set to match that of the
*cutoff* keyword.

The neighbor list needed to compute this quantity is constructed each
time the calculation is performed (i.e. each time a snapshot of atoms is
dumped). Thus it can be inefficient to compute/dump this quantity too
frequently.

:::: note
::: title
Note
:::

If you have a bonded system, then the settings of
[special_bonds](special_bonds) command can remove pairwise interactions
between atoms in the same bond, angle, or dihedral. This is the default
setting for the [special_bonds](special_bonds) command, and means those
pairwise interactions do not appear in the neighbor list. Because this
compute uses the neighbor list, it also means those pairs will not be
included in the order parameter. This difficulty can be circumvented by
writing a dump file, and using the [rerun](rerun) command to compute the
order parameter for snapshots in the dump file. The rerun script can use
a [special_bonds](special_bonds) command that includes all pairs in the
neighbor list.
::::

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

## Output info

This compute calculates a per-atom array with two columns: mass density
in density [units](units) and temperature in temperature [units](units).

These values can be accessed by any command that uses per-atom values
from a compute as input. See the [Howto output](Howto_output) doc page
for an overview of LAMMPS output options.

## Restrictions

This compute is part of the EXTRA-COMPUTE package. It is only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

## Related commands

[comm_modify](comm_modify)

## Default

The option defaults are *cutoff* = pair style cutoff.
