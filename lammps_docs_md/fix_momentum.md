# fix momentum command

Accelerator Variants: *momentum/kk*

# fix momentum/chunk command

## Syntax

    fix ID group-ID momentum N keyword values ...

-   ID, group-ID are documented in [fix](fix) command
-   momentum = style name of this fix command
-   N = adjust the momentum every this many timesteps one or more
    keyword/value pairs may be appended

<!-- -->

    fix ID group-ID momentum/chunk N chunkID keyword values ...

-   ID, group-ID are documented in [fix](fix) command

-   momentum/chunk = style name of this fix command

-   N = adjust the momentum per chunk every this many timesteps

-   chunkID = ID of [compute chunk/atom](compute_chunk_atom) command

    one or more keyword/value settings may be appended to each of the
    fix commands:

-   keyword = *linear* or *angular* or *rescale*

        *linear* values = xflag yflag zflag
          xflag,yflag,zflag = 0/1 to exclude/include each dimension
        *angular* values = none
        *rescale* values = none

## Examples

``` LAMMPS
fix 1 all momentum 1 linear 1 1 0
fix 1 all momentum 1 linear 1 1 1 rescale
fix 1 all momentum 100 linear 1 1 1 angular
fix 1 all momentum/chunk 100 molchunk linear 1 1 1 angular
```

## Description

Fix momentum zeroes the linear and/or angular momentum of the group of
atoms every N timesteps by adjusting the velocities of the atoms. Fix
momentum/chunk works equivalently, but operates on a per-chunk basis.

One (or both) of the *linear* or *angular* keywords **must** be
specified.

If the *linear* keyword is used, the linear momentum is zeroed by
subtracting the center-of-mass velocity of the group or chunk from each
atom. This does not change the relative velocity of any pair of atoms.
One or more dimensions can be excluded from this operation by setting
the corresponding flag to 0.

If the *angular* keyword is used, the angular momentum is zeroed by
subtracting a rotational component from each atom.

This command can be used to ensure the entire collection of atoms (or a
subset of them) does not drift or rotate during the simulation due to
random perturbations (e.g. [fix langevin](fix_langevin) thermostatting).

The *rescale* keyword enables conserving the kinetic energy of the group
or chunk of atoms by rescaling the velocities after the momentum was
removed.

Note that the [velocity](velocity) command can be used to create initial
velocities with zero aggregate linear and/or angular momentum.

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

## Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to [binary restart
files](restart). None of the [fix_modify](fix_modify) options are
relevant to this fix. No global or per-atom quantities are stored by
this fix for access by various [output commands](Howto_output). No
parameter of this fix can be used with the *start/stop* keywords of the
[run](run) command. This fix is not invoked during [energy
minimization](minimize).

## Restrictions

Fix momentum/chunk is part of the EXTRA-FIX package. It is only enabled
if LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

## Related commands

[fix recenter](fix_recenter), [velocity](velocity)

## Default

none
