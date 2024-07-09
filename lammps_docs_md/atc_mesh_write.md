# fix_modify AtC mesh write command

## Syntax

    fix_modify <AtC fixID> mesh write <f|p> <f|p> <f|p>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   mesh write = name of the AtC sub-command
-   filename = name of the file containing the mesh to be write

## Examples

``` LAMMPS
fix_modify AtC mesh write myMesh.mesh
```

## Description

Writes a mesh to a text file.

## Restrictions

None

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix_modify AtC mesh create](atc_mesh_create)
-   [fix_modify AtC mesh read](atc_mesh_read)

## Default

None.
