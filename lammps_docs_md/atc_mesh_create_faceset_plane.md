# fix_modify AtC mesh create_faceset plane command

## Syntax

    fix_modify <AtC fixID> mesh create_faceset <id> plane <x|y|z> <val1> <x|y|z> <lval2> <uval2> [units]

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   mesh create_faceset = name of the AtC sub-command
-   id = id to assign to the collection of FE faces
-   plane = use plane to define faceset
-   \<val1\>,\<lval2\>,\<uval2\> = plane is specified as the xz=val1
    plane bounded by the segments xz = \[lval2,uval2\]
-   units = option to specify real as opposed to lattice units

## Examples

``` LAMMPS
fix_modify AtC mesh create_faceset xyplane plane y 0 x -4 0
```

## Description

Command to assign an id to a set of FE faces.

## Restrictions

Only viable for rectangular grids.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix_modify AtC mesh create_faceset
    box](atc_mesh_create_faceset_box)

## Default

The default options are units = lattice.
