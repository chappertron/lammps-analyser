# fix nve/asphere/noforce command

## Syntax

    fix ID group-ID nve/asphere/noforce

-   ID, group-ID are documented in [fix](fix) command
-   nve/asphere/noforce = style name of this fix command

## Examples

``` LAMMPS
fix 1 all nve/asphere/noforce
```

## Description

Perform updates of position and orientation, but not velocity or angular
momentum for atoms in the group each timestep. In other words, the force
and torque on the atoms is ignored and their velocity and angular
momentum are not updated. The atom velocities and angular momenta are
used to update their positions and orientation.

This is useful as an implicit time integrator for Fast Lubrication
Dynamics, since the velocity and angular momentum are updated by the
[pair_style lubricuteU](pair_lubricateU) command.

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

This fix is part of the ASPHERE package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

This fix requires that atoms store torque and angular momentum and a
quaternion as defined by the [atom_style ellipsoid](atom_style) command.

All particles in the group must be finite-size. They cannot be point
particles, but they can be aspherical or spherical as defined by their
shape attribute.

## Related commands

[fix nve/noforce](fix_nve_noforce), [fix nve/asphere](fix_nve_asphere)

## Default

none
