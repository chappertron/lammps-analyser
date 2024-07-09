# echo command

## Syntax

``` LAMMPS
echo style
```

-   style = *none* or *screen* or *log* or *both*

## Examples

``` LAMMPS
echo both
echo log
```

## Description

This command determines whether LAMMPS echoes each input script command
to the screen and/or log file as it is read and processed. If an input
script has errors, it can be useful to look at echoed output to see the
last command processed.

The [command-line switch](Run_options) -echo can be used in place of
this command.

## Restrictions

> none

## Related commands

none

## Default

``` LAMMPS
echo log
```
