# fix_modify AtC filter scale command

## Syntax

    fix_modify <AtC fixID> filter scale <scale>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   filter scale = name of the AtC sub-command
-   scale = characteristic times scale of the filter

## Examples

``` LAMMPS
fix_modify AtC filter scale 10.0
```

## Description

Sets the time scale for MD dynamics filter to construct a more
appropriate continuous field.

## Restrictions

Only for use with these specific transfers: thermal, two_temperature

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix_modify AtC filter](atc_time_filter)
-   [fix_modify AtC filter type](atc_filter_type)

## Default

0.0
