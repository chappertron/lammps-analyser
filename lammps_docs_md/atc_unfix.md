# fix_modify AtC unfix command

## Syntax

    fix_modify <AtC fixID> unfix <field> <nodeset>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   unfix = name of the AtC sub-command
-   field = field kind name valid for type of physics: *temperature* or
    *electron_temperature*
-   nodeset = name of set of nodes to apply boundary condition

## Examples

``` LAMMPS
fix_modify AtC unfix temperature groupNAME
```

## Description

Removes constraint on field values for specified nodes.

## Restrictions

The keyword *all* is reserved and thus not available as nodeset name.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix_modify AtC fix](atc_fix)

## Default

None.
