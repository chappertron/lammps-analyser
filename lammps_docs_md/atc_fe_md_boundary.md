# fix_modify AtC fe_md_boundary command

## Syntax

    fix_modify <AtC fixID> fe_md_boundary <faceset|interpolate|no_boundary>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   fe_md_boundary = name of the AtC sub-command
-   *faceset* or *interpolate* or *no_boundary*

## Examples

``` LAMMPS
fix_modify AtC fe_md_boundary interpolate
```

## Description

Specifies different methods for computing fluxes between between the MD
and FE integration regions. Faceset defines a faceset separating the MD
and FE regions and uses finite element face quadrature to compute the
flux. Interpolate uses a reconstruction scheme to approximate the flux,
which is more robust but less accurate if the MD/FE boundary does
correspond to a faceset. No boundary results in no fluxes between the
systems being computed.

## Restrictions

If *faceset* is used, all the AtC non-boundary atoms must lie within and
completely fill the domain enclosed by the faceset.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix_modify AtC mesh create_faceset
    box](atc_mesh_create_faceset_box)
-   [fix_modify AtC mesh create_faceset
    plane](atc_mesh_create_faceset_plane)

## Default

*interpolate*
