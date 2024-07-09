# fix_modify AtC equilibrium_start command

## Syntax

    fix_modify <AtC fixID> equilibrium_start <on|off>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   equilibrium_start = name of the AtC sub-command
-   *exponential* or *step* or *no_filter* = select type of filter

## Examples

``` LAMMPS
fix_modify AtC equilibrium_start on
```

## Description

Starts filtered calculations assuming they start in equilibrium, i.e.
perfect finite element force balance.

## Restrictions

Only for use with these specific transfers: thermal, two_temperature

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix_modify AtC filter](atc_time_filter)
-   [fix_modify AtC filter scale](atc_filter_scale)

## Default

None.
