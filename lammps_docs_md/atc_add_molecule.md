# fix_modify AtC add_molecule command

## Syntax

    fix_modify <AtC fixID> add_molecule <small|large> <tag> <group-ID>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   add_molecule = name of the AtC sub-command
-   *small* or *large* = can be *small* if molecule size \< cutoff
    radius, must be *large* otherwise
-   tag = tag for tracking a molecule
-   *group-ID* = LAMMPS defined group-ID

## Examples

``` LAMMPS
group WATERGROUP type 1 2
fix_modify AtC add_molecule small water WATERGROUP
```

## Description

Associates a tag with all molecules corresponding to a specified group.

## Restrictions

None.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix_modify AtC add_species](atc_add_species)
-   [fix_modify AtC remove_species](atc_remove_species)
-   [fix_modify AtC remove_molecule](atc_remove_molecule)

## Default

None.
