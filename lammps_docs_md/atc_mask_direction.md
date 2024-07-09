# fix_modify AtC control mask_direction command

## Syntax

    fix_modify <AtC fixID> control mask_direction <direction> <on|off>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   control mask_direction = name of the AtC sub-command
-   direction = select direction
-   *on* or *off* = Toggles state

## Examples

``` LAMMPS
fix_modify AtC control mask_direction 0 on
```

## Description

Command to mask out certain dimensions from the atomic regulator

## Restrictions

None.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
