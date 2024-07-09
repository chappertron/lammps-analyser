# fix sph command

## Syntax

    fix ID group-ID sph

-   ID, group-ID are documented in [fix](fix) command
-   sph = style name of this fix command

## Examples

``` LAMMPS
fix 1 all sph
```

## Description

Perform time integration to update position, velocity, internal energy
and local density for atoms in the group each timestep. This fix is
needed to time-integrate SPH systems where particles carry internal
variables such as internal energy. SPH stands for Smoothed Particle
Hydrodynamics.

See [this PDF guide](PDF/SPH_LAMMPS_userguide.pdf)\_ to using SPH in
LAMMPS.

## Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to [binary restart
files](restart). None of the [fix_modify](fix_modify) options are
relevant to this fix. No global or per-atom quantities are stored by
this fix for access by various [output commands](Howto_output). No
parameter of this fix can be used with the *start/stop* keywords of the
[run](run) command. This fix is not invoked during [energy
minimization](minimize).

## Restrictions

This fix is part of the SPH package. It is only enabled if LAMMPS was
built with that package. See the [Build package](Build_package) page for
more info.

## Related commands

[fix sph/stationary](fix_sph_stationary)

## Default

none
