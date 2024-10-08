# fix lineforce command

## Syntax

    fix ID group-ID lineforce x y z

-   ID, group-ID are documented in [fix](fix) command
-   lineforce = style name of this fix command
-   x y z = direction of line as a 3-vector

## Examples

``` LAMMPS
fix hold boundary lineforce 0.0 1.0 1.0
```

## Description

Adjust the forces on each atom in the group so that only the component
of force along the linear direction specified by the vector (x,y,z)
remains. This is done by subtracting out components of force in the
plane perpendicular to the line.

If the initial velocity of the atom is 0.0 (or along the line), then it
should continue to move along the line thereafter.

## Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to [binary restart
files](restart). None of the [fix_modify](fix_modify) options are
relevant to this fix. No global or per-atom quantities are stored by
this fix for access by various [output commands](Howto_output). No
parameter of this fix can be used with the *start/stop* keywords of the
[run](run) command.

The forces due to this fix are imposed during an energy minimization,
invoked by the [minimize](minimize) command.

## Restrictions

> none

## Related commands

[fix planeforce](fix_planeforce)

## Default

none
