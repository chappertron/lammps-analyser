# fix_modify AtC filter command

## Syntax

    fix_modify <AtC fixID> filter <on|off|equilibrate>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   filter = name of the AtC sub-command
-   *on* or *off* or *equilibrate* = Select state of filter

## Examples

``` LAMMPS
fix_modify AtC filter on
```

## Description

Filters the MD dynamics to construct a more appropriate continuous
field. Equilibrating first filters the time derivatives without changing
the dynamics to provide a better initial condition to the filtered
dynamics.

## Restrictions

Only for use with these specific transfers: thermal, two_temperature

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix_modify AtC filter scale](atc_filter_scale)
-   [fix_modify AtC equilibrium_start](atc_equilibrium_start)

## Default

off
