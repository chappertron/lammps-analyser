# compute smd/contact/radius command

## Syntax

    compute ID group-ID smd/contact/radius

-   ID, group-ID are documented in [compute](compute) command
-   smd/contact/radius = style name of this compute command

## Examples

``` LAMMPS
compute 1 all smd/contact/radius
```

## Description

Define a computation which outputs the contact radius, i.e., the radius
used to prevent particles from penetrating each other. The contact
radius is used only to prevent particles belonging to different physical
bodies from penetrating each other. It is used by the contact pair
styles, e.g., smd/hertz and smd/tri_surface.

See [this PDF guide](PDF/MACHDYN_LAMMPS_userguide.pdf)\_ to using Smooth
Mach Dynamics in LAMMPS.

The value of the contact radius will be 0.0 for particles not in the
specified compute group.

## Output info

This compute calculates a per-particle vector, which can be accessed by
any command that uses per-particle values from a compute as input. See
the [Howto output](Howto_output) page for an overview of LAMMPS output
options.

The per-particle vector values will be in distance [units](units).

## Restrictions

This compute is part of the MACHDYN package. It is only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

## Related commands

[dump custom](dump) smd/hertz smd/tri_surface

## Default

none
