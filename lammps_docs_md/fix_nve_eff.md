# fix nve/eff command

## Syntax

    fix ID group-ID nve/eff

-   ID, group-ID are documented in [fix](fix) command
-   nve/eff = style name of this fix command

## Examples

``` LAMMPS
fix 1 all nve/eff
```

## Description

Perform constant NVE integration to update position and velocity for
nuclei and electrons in the group for the [electron force
field](pair_eff) model. V is volume; E is energy. This creates a system
trajectory consistent with the microcanonical ensemble.

The operation of this fix is exactly like that described by the [fix
nve](fix_nve) command, except that the radius and radial velocity of
electrons are also updated.

## Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to [binary restart
files](restart). None of the [fix_modify](fix_modify) options are
relevant to this fix. No global or per-atom quantities are stored by
this fix for access by various [output commands](Howto_output). No
parameter of this fix can be used with the *start/stop* keywords of the
[run](run) command. This fix is not invoked during [energy
minimization](minimize).

## Restrictions

This fix is part of the EFF package. It is only enabled if LAMMPS was
built with that package. See the [Build package](Build_package) page for
more info.

## Related commands

[fix nve](fix_nve), [fix nvt/eff](fix_nh_eff), [fix npt/eff](fix_nh_eff)

## Default

none
