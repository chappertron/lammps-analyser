# fix_modify AtC atom_element_map command

## Syntax

    fix_modify <AtC fixID> atom_element_map <eulerian|lagrangian> [<frequency>]

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   atom_element_map = name of the AtC sub-command
-   *eulerian* or *lagrangian* = frame of reference
-   frequency = frequency of updating atom-to-continuum maps based on
    the current configuration - (only for eulerian)

## Examples

``` LAMMPS
fix_modify AtC atom_element_map eulerian 100
```

## Description

Changes frame of reference from *eulerian* to *lagrangian* or vice versa
and sets the frequency for which the map from atoms to elements is
reformed and all the attendant data is recalculated.

## Restrictions

Cannot change map type after initialization.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)

## Default

*lagrangian*
