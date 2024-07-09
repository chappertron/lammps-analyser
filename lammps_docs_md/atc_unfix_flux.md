# fix_modify AtC unfix_flux command

## Syntax

    fix_modify <AtC fixID> unfix_flux <field> <face_set> <value|function>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   unfix_flux = name of the AtC sub-command
-   field = field kind name valid for type of physics: *temperature* or
    *electron_temperature*
-   face_set = name of set of element faces

## Examples

``` LAMMPS
fix_modify AtC unfix_flux temperature faceSet
```

## Description

Command for removing prescribed normal fluxes e.g. heat_flux, stress.

## Restrictions

None.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix_modify AtC fix_flux](atc_fix_flux)

## Default

None.
