# fix_modify AtC write_restart command

## Syntax

    fix_modify <AtC fixID> write_restart <file_name>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   write_restart = name of the AtC sub-command
-   file_name = name of AtC restart file

## Examples

``` LAMMPS
fix_modify AtC write_restart restart.mydata.AtC
```

## Description

Dumps the current state of the fields to a named text-based restart
file. This done when the command is invoked and not repeated, unlike the
otherwise similar LAMMPS command.

## Restrictions

None.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix_modify AtC read_restart](atc_read_restart)

## Default

None.
