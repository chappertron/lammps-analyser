# fix_modify AtC material command

## Syntax

    fix_modify <AtC fixID> material <elementset_name> <material_id>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   material = name of the AtC sub-command
-   elementset_name = name of the elementset
-   material_id = ID of the material

## Examples

``` LAMMPS
fix_modify AtC material gap_region 1
```

## Description

Sets the material model in *elementset_name* to be of type
*material_id*.

## Restrictions

The element set must already be created and the material must be
specified in the material file given the the atc fix on construction

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)

## Default

All elements default to the first material in the material file.
