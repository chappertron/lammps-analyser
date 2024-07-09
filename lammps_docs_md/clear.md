# clear command

## Syntax

``` LAMMPS
clear
```

## Examples

``` LAMMPS
# (commands for 1st simulation)
clear
# (commands for 2nd simulation)
```

## Description

This command deletes all atoms, restores all settings to their default
values, and frees all memory allocated by LAMMPS. Once a clear command
has been executed, it is almost as if LAMMPS were starting over, with
only the exceptions noted below. This command enables multiple jobs to
be run sequentially from one input script.

These settings are not affected by a clear command: the working
directory ([shell](shell) command), log file status ([log](log)
command), echo status ([echo](echo) command), and input script variables
except for *atomfile* style variables ([variable](variable) command).

## Restrictions

> none

## Related commands

none

## Default

none
