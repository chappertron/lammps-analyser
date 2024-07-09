# compute erotate/asphere command

## Syntax

``` LAMMPS
compute ID group-ID erotate/asphere
```

-   ID, group-ID are documented in [compute](compute) command
-   erotate/asphere = style name of this compute command

## Examples

``` LAMMPS
compute 1 all erotate/asphere
```

## Description

Define a computation that calculates the rotational kinetic energy of a
group of aspherical particles. The aspherical particles can be
ellipsoids, or line segments, or triangles. See the
[atom_style](atom_style) and [read_data](read_data) commands for
descriptions of these options.

For all 3 types of particles, the rotational kinetic energy is computed
as $\frac12 I \omega^2$, where $I$ is the inertia tensor for the
aspherical particle and $\omega$ is its angular velocity, which is
computed from its angular momentum if needed.

:::: note
::: title
Note
:::

For [2d models](dimension), ellipsoidal particles are treated as
ellipsoids, not ellipses, meaning their moments of inertia will be the
same as in 3d.
::::

## Output info

This compute calculates a global scalar (the KE). This value can be used
by any command that uses a global scalar value from a compute as input.
See the [Howto output](Howto_output) page for an overview of LAMMPS
output options.

The scalar value calculated by this compute is \"extensive\". The scalar
value will be in energy [units](units).

## Restrictions

This compute requires that ellipsoidal particles atoms store a shape and
quaternion orientation and angular momentum as defined by the
[atom_style ellipsoid](atom_style) command.

This compute requires that line segment particles atoms store a length
and orientation and angular velocity as defined by the [atom_style
line](atom_style) command.

This compute requires that triangular particles atoms store a size and
shape and quaternion orientation and angular momentum as defined by the
[atom_style tri](atom_style) command.

All particles in the group must be of finite size. They cannot be point
particles.

## Related commands

none

[compute erotate/sphere](compute_erotate_sphere)

## Default

none
