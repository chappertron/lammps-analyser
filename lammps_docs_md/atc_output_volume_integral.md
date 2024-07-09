# fix_modify AtC output volume_integral command

## Syntax

    fix_modify <AtC fixID> output volume_integral <elementset_name> <field>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   output volume_integral = name of the AtC sub-command
-   elementset_name= name of elementset to be integrated over
-   fieldname = name of field to integrate

## Examples

``` LAMMPS
fix_modify AtC output volume_integral eset1 mass_density
```

## Description

Performs volume integration of specified field over elementset and
outputs resulting variable values to GLOBALS file.

## Restrictions

None.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix atc command](fix_atc)

## Default

None.
