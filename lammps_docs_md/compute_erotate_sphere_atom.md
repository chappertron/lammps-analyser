# compute erotate/sphere/atom command

## Syntax

``` LAMMPS
compute ID group-ID erotate/sphere/atom
```

-   ID, group-ID are documented in [compute](compute) command
-   erotate/sphere/atom = style name of this compute command

## Examples

``` LAMMPS
compute 1 all erotate/sphere/atom
```

## Description

Define a computation that calculates the rotational kinetic energy for
each particle in a group.

The rotational energy is computed as $\frac12 I \omega^2$, where $I$ is
the moment of inertia for a sphere and $\omega$ is the particle\'s
angular velocity.

:::: note
::: title
Note
:::

For [2d models](dimension), particles are treated as spheres, not disks,
meaning their moment of inertia will be the same as in 3d.
::::

The value of the rotational kinetic energy will be 0.0 for atoms not in
the specified compute group or for point particles with a radius of 0.0.

## Output info

This compute calculates a per-atom vector, which can be accessed by any
command that uses per-atom values from a compute as input. See the
[Howto output](Howto_output) page for an overview of LAMMPS output
options.

The per-atom vector values will be in energy [units](units).

## Restrictions

> none

## Related commands

[dump custom](dump)

## Default

none
