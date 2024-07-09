# fix_modify AtC remove_source command

## Syntax

    fix_modify <AtC fixID> remove_source <field> <element_set>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   remove_source = name of the AtC sub-command
-   field = field kind name valid for type of physics: *temperature* or
    *electron_temperature*
-   element_set = name of set of elements

## Examples

``` LAMMPS
fix_modify AtC remove_source temperature groupNAME
```

## Description

Remove a domain source.

## Restrictions

The keyword *all* is reserved and thus not available as element_set
name.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix_modify AtC source](atc_source)

## Default

None.
