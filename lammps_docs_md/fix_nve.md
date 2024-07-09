# fix nve command

Accelerator Variants: *nve/gpu*, *nve/intel*, *nve/kk*, *nve/omp*

## Syntax

    fix ID group-ID nve

-   ID, group-ID are documented in [fix](fix) command
-   nve = style name of this fix command

## Examples

``` LAMMPS
fix 1 all nve
```

## Description

Perform plain time integration to update position and velocity for atoms
in the group each timestep. This creates a system trajectory consistent
with the microcanonical ensemble (NVE) provided there are (full)
periodic boundary conditions and no other \"manipulations\" of the
system (e.g. fixes that modify forces or velocities).

This fix invokes the velocity form of the Stoermer-Verlet time
integration algorithm (velocity-Verlet). Other time integration options
can be invoked using the [run_style](run_style) command.

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
files](restart). None of the [fix_modify](fix_modify) options are
relevant to this fix. No global or per-atom quantities are stored by
this fix for access by various [output commands](Howto_output). No
parameter of this fix can be used with the *start/stop* keywords of the
[run](run) command. This fix is not invoked during [energy
minimization](minimize).

## Restrictions

> none

## Related commands

[fix nvt](fix_nh), [fix npt](fix_nh), [run_style](run_style)

## Default

none
