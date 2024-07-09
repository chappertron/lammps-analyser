# fix_modify AtC mesh delete_elements command

## Syntax

    fix_modify <AtC fixID> mesh delete_elements <id>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   mesh create_elementset = name of the AtC sub-command
-   id = id of the element set

## Examples

``` LAMMPS
fix_modify AtC mesh delete_elements gap
```

## Description

Deletes a group of elements from the mesh.

## Restrictions

None.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix_modify AtC mesh create_elementset](atc_mesh_create_elementset)
-   [fix_modify AtC mesh
    nodeset_to_elementset](atc_mesh_nodeset_to_elementset)

## Default

None.
