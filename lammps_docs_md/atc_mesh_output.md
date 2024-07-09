# fix_modify AtC mesh output command

## Syntax

    fix_modify <AtC fixID> mesh output <file_prefix>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   mesh output = name of the AtC sub-command
-   file_prefix = prefix of various generated output files

## Examples

``` LAMMPS
fix_modify AtC mesh output meshData
```

## Description

Command to output mesh and associated data: nodesets, facesets, and
elementsets. This data is only output once upon initialization since
currently the mesh is static. Creates binary (EnSight, \"gold\" format)
output of mesh data.

## Restrictions

None.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)

## Default

None.
