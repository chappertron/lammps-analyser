# fix_modify AtC extrinsic exchange command

## Syntax

    fix_modify <AtC fixID> extrinsic exchange <on|off>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   extrinsic exchange = name of the AtC sub-command
-   *on* or *off* = set state of energy exchange

## Examples

``` LAMMPS
fix_modify AtC extrinsic exchange on
```

## Description

Switches energy exchange between the MD system and the electron system
on or off

## Restrictions

For use only with the two_temperature type of the AtC fix (see [fix
atc](fix_atc) command)

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)

## Default

*on*
