# fix_modify AtC mesh add_to_nodeset command

## Syntax

    fix_modify <AtC fixID> mesh add_to_nodeset <id> <xmin> <xmax> <ymin> <ymax> <zmin> <zmax>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   mesh create_nodeset = name of the AtC sub-command
-   id = id to assign to the collection of FE nodes
-   \<xmin\> \<xmax\> \<ymin\> \<ymax\> \<zmin\> \<zmax\> = coordinates
    of the bounding box that contains the desired nodes to be added

## Examples

``` LAMMPS
fix_modify AtC mesh add_to_nodeset lbc -12.1 -11.9 -12 12 -12 12
```

## Description

Command to add nodes to an already existing FE nodeset.

## Restrictions

None

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix_modify AtC mesh create_nodeset](atc_mesh_create_nodeset)

## Default

Coordinates are assumed to be in lattice units.
