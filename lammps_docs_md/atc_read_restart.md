# fix_modify AtC read_restart command

## Syntax

    fix_modify <AtC fixID> read_restart <file_name>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   read_restart = name of the AtC sub-command
-   file_name = name of AtC restart file

## Examples

``` LAMMPS
fix_modify AtC read_restart restart.mydata.AtC
```

## Description

Reads the current state of the AtC fields from a named text-based
restart file.

## Restrictions

The restart file only contains fields and their time derivatives. The
reference positions of the atoms and the commands that initialize the
fix are not saved e.g. an identical mesh containing the same atoms will
have to be recreated.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix_modify AtC write_restart](atc_write_restart)

## Default

None.
