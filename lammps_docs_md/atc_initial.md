# fix_modify AtC initial command

## Syntax

    fix_modify <AtC fixID> initial <field> <nodeset> <constant|function>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   initial = name of the AtC sub-command
-   field = field kind name valid for type of physics: *temperature* or
    *electron_temperature*
-   nodeset = name of set of nodes to apply initial condition
-   *constant* or *function* = value or name of function followed by its
    parameters

## Examples

``` LAMMPS
fix_modify AtC initial temperature groupNAME 10.
```

## Description

Sets the initial values for the specified field at the specified nodes.

## Restrictions

The keyword *all* is reserved and thus not available as nodeset name.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)

## Default

None.
