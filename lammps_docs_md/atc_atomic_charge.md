# fix_modify AtC atomic_charge command

## Syntax

    fix_modify <AtC fixID> <include|omit> atomic_charge

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   *include* or *omit* = switch to activate/deactivate inclusion of
    intrinsic atomic charge in ATC
-   atomic_charge = name of the AtC sub-command

## Examples

``` LAMMPS
fix_modify AtC include atomic_charge
```

## Description

Determines whether AtC tracks the total charge as a finite element
field.

## Restrictions

Required for: *electrostatics*

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)

## Default

If the atom charge is defined, default is on, otherwise default is off.
