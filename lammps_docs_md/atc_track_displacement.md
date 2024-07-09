# fix_modify AtC track_displacement command

## Syntax

    fix_modify <AtC fixID> track_displacement <on|off>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   track_displacement = name of the AtC sub-command
-   *on* or *off* = (undocumented)

## Examples

``` LAMMPS
fix_modify AtC track_displacement on
```

## Description

Determines whether displacement is tracked or not. For solids problems
this is a useful quantity, but for fluids it is not relevant.

## Restrictions

Some constitutive models require the displacement field.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)

## Default

*on*
