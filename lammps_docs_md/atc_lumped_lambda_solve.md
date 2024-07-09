# fix_modify AtC control lumped_lambda_solve command

## Syntax

    fix_modify <AtC fixID> control lumped_lambda_solve <on|off>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   control lumped_lambda_solve = name of the AtC sub-command
-   *on* or *off* = Toggles state of lumped matrix

## Examples

``` LAMMPS
fix_modify AtC control lumped_lambda_solve on
```

## Description

Command select whether to use or not use lumped matrix for lambda solve.

## Restrictions

None.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)

## Default

off
