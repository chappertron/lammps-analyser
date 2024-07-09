# fix_modify AtC mesh create command

## Syntax

    fix_modify <AtC fixID> mesh create <nx> <ny> <nz> <region-ID> <f|p> <f|p> <f|p>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   mesh create = name of the AtC sub-command
-   nx ny nz = number of elements in x-, y-, and z-direction
-   region-ID = ID of region that is to be meshed
-   f or p = periodicity flags for x-, y-, and z-direction

## Examples

``` LAMMPS
fix_modify AtC mesh create 10 1 1 feRegion p p p
```

## Description

Creates a uniform mesh in a rectangular region.

## Restrictions

Creates only uniform rectangular grids in a rectangular region

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix_modify AtC mesh quadrature](atc_mesh_quadrature)

## Default

When created, the mesh defaults to gauss2 (2-point Gaussian) quadrature.
Use the [fix_modify AtC mesh quadrature](atc_mesh_quadrature) command to
change the quadrature style.
