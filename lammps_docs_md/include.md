# include command

## Syntax

    include file

-   file = filename of new input script to switch to

## Examples

``` LAMMPS
include newfile
include in.run2
```

## Description

This command opens a new input script file and begins reading LAMMPS
commands from that file. When the new file is finished, the original
file is returned to. Include files can be nested as deeply as desired.
If input script A includes script B, and B includes A, then LAMMPS could
run for a long time.

If the filename is a variable (see the [variable](variable) command),
different processor partitions can run different input scripts.

## Restrictions

> none

## Related commands

[variable](variable), [jump](jump)

## Default

none
