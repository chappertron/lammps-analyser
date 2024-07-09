# fix_modify AtC output nodeset command

## Syntax

    fix_modify <AtC fixID> output nodeset <nodeset_name> <operation>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   output nodeset = name of the AtC sub-command
-   nodeset_name= name of nodeset to be operated on
-   operation = *sum*
    -   *sum* = creates nodal sum over nodes in specified nodeset

## Examples

``` LAMMPS
fix_modify AtC output nodeset nset1 sum
```

## Description

Performs operation over the nodes belonging to specified nodeset and
outputs resulting variable values to GLOBALS file.

## Restrictions

None.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix atc command](fix_atc)

## Default

None.
