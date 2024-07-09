# fix_modify AtC boundary_faceset command

## Syntax

    fix_modify <AtC fixID> boundary_faceset <is|add> <faceset_name>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   boundary_faceset = name of the AtC sub-command
-   *is* or *add* = select whether to select or add a faceset
-   faceset_name = name of the faceset

## Examples

``` LAMMPS
fix_modify AtC boundary_faceset is obndy
```

## Description

This command species the faceset name when using a faceset to compute
the MD/FE boundary fluxes. The faceset must already exist.

## Restrictions

This is only valid when *fe_md_boundary* is set to *faceset*.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix_modify AtC fe_md_boundary](atc_fe_md_boundary)
-   [fix_modify AtC mesh create_faceset
    box](atc_mesh_create_faceset_box)
-   [fix_modify AtC mesh create_faceset
    plane](atc_mesh_create_faceset_plane)

## Default

None.
