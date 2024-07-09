# fix nve/asphere command

Accelerator Variants: *nve/asphere/gpu*, *nve/asphere/intel*

## Syntax

    fix ID group-ID nve/asphere

-   ID, group-ID are documented in [fix](fix) command
-   nve/asphere = style name of this fix command

## Examples

``` LAMMPS
fix 1 all nve/asphere
```

## Description

Perform constant NVE integration to update position, velocity,
orientation, and angular velocity for aspherical particles in the group
each timestep. V is volume; E is energy. This creates a system
trajectory consistent with the microcanonical ensemble.

This fix differs from the [fix nve](fix_nve) command, which assumes
point particles and only updates their position and velocity.

## Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to [binary restart
files](restart). None of the [fix_modify](fix_modify) options are
relevant to this fix. No global or per-atom quantities are stored by
this fix for access by various [output commands](Howto_output). No
parameter of this fix can be used with the *start/stop* keywords of the
[run](run) command. This fix is not invoked during [energy
minimization](minimize).

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

## Restrictions

This fix is part of the ASPHERE package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

This fix requires that atoms store torque and angular momentum and a
quaternion as defined by the [atom_style ellipsoid](atom_style) command.

All particles in the group must be finite-size. They cannot be point
particles, but they can be aspherical or spherical as defined by their
shape attribute.

## Related commands

[fix nve](fix_nve), [fix nve/sphere](fix_nve_sphere)

## Default

none
