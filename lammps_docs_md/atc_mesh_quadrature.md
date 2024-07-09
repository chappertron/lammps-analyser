# fix_modify AtC mesh quadrature command

## Syntax

    fix_modify <AtC fixID> mesh quatrature <quad>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   mesh quadrature = name of the AtC sub-command
-   quad = *nodal* or *gauss1* or *gauss2* or *gauss3* or *face*

## Examples

``` LAMMPS
fix_modify AtC mesh quadrature face
```

## Description

(Re-)assigns the quadrature style for an existing mesh. When a mesh is
created its quadrature method defaults to gauss2. Use this call to
change it after the fact.

## Restrictions

None.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix_modify AtC mesh create](atc_mesh_create)

## Default

None.
