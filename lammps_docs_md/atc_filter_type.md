# fix_modify AtC filter type command

## Syntax

    fix_modify <AtC fixID> filter type <exponential|step|no_filter>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   filter type = name of the AtC sub-command
-   *exponential* or *step* or *no_filter* = select type of filter

## Examples

``` LAMMPS
fix_modify AtC filter type exponential
```

## Description

Specifies the type of time filter used.

## Restrictions

Only for use with these specific transfers: thermal, two_temperature

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix_modify AtC filter](atc_time_filter)
-   [fix_modify AtC filter scale](atc_filter_scale)

## Default

None.
