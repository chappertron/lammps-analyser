# fix_modify AtC temperature_definition command

## Syntax

    fix_modify <AtC fixID> temperature_definition <kinetic|total>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   temperature_definition = name of the AtC sub-command
-   *kinetic* or *total* = (undocumented)

## Examples

``` LAMMPS
fix_modify AtC temperature_definition kinetic
```

## Description

Change the definition for the atomic temperature used to create the
finite element temperature. The kinetic option is based only on the
kinetic energy of the atoms while the total option uses the total energy
(kinetic + potential) of an atom.

## Restrictions

This command is only valid when using thermal coupling. Also, while not
a formal restriction, the user should ensure that associating a
potential energy with each atom makes physical sense for the total
option to be meaningful.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)

## Default

*kinetic*
