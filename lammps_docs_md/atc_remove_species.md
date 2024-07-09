# fix_modify AtC remove_species command

## Syntax

    fix_modify <AtC fixID> remove_species <tag>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   remove_species = name of the AtC sub-command
-   tag = tag for tracking a species

## Examples

``` LAMMPS
fix_modify AtC remove_species gold
```

## Description

Removes tag designated for tracking a specified species.

## Restrictions

None.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix_modify AtC add_species](atc_add_species)
-   [fix_modify AtC add_molecule](atc_add_molecule)
-   [fix_modify AtC remove_molecule](atc_remove_molecule)

## Default

None.
