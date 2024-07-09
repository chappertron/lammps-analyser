# fix_modify AtC mesh read command

## Syntax

    fix_modify <AtC fixID> mesh read <f|p> <f|p> <f|p>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   mesh read = name of the AtC sub-command
-   filename = name of the file containing the mesh to be read
-   f or p = periodicity flags for x-, y-, and z-direction (optional)

## Examples

``` LAMMPS
fix_modify AtC mesh read myComponent.mesh p p p
fix_modify AtC mesh read myOtherComponent.exo
```

## Description

Reads a mesh from a text or exodus file, and assigns periodic boundary
conditions if needed.

## Restrictions

None

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix_modify AtC mesh create](atc_mesh_create)
-   [fix_modify AtC mesh write](atc_mesh_write)

## Default

Periodicity flags are set to false (f) by default.
