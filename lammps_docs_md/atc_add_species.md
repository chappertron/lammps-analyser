# fix_modify AtC add_species command

## Syntax

    fix_modify <AtC fixID> add_species <tag> <group|type> <ID>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   add_species = name of the AtC sub-command
-   tag = tag for tracking a species
-   *group* or *type* = LAMMPS defined group or type of atoms
-   ID = name of group or type number

## Examples

``` LAMMPS
fix_modify AtC add_species gold type 1
group GOLDGROUP type 1
fix_modify AtC add_species gold group GOLDGROUP
```

## Description

Associates a tag with all atoms of a specified type or within a
specified group.

## Restrictions

None.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix_modify AtC add_molecule](atc_add_molecule)
-   [fix_modify AtC remove_species](atc_remove_species)
-   [fix_modify AtC remove_molecule](atc_remove_molecule)

## Default

None.
