# jump command

## Syntax

    jump file label

-   file = filename of new input script to switch to
-   label = optional label within file to jump to

## Examples

``` LAMMPS
jump newfile
jump in.run2 runloop
jump SELF runloop
```

## Description

This command closes the current input script file, opens the file with
the specified name, and begins reading LAMMPS commands from that file.
Unlike the [include](include) command, the original file is not returned
to, although by using multiple jump commands it is possible to chain
from file to file or back to the original file.

If the word \"SELF\" is used for the filename, then the current input
script is re-opened and read again.

:::: note
::: title
Note
:::

The SELF option is not guaranteed to work when the current input script
is being read through stdin (standard input), e.g.
::::

``` bash
lmp_g++ < in.script
```

since the SELF option invokes the C-library rewind() call, which may not
be supported for stdin on some systems or by some MPI implementations.
This can be worked around by using the [-in command-line
switch](Run_options), e.g.

``` bash
lmp_g++ -in in.script
```

or by using the [-var command-line switch](Run_options) to pass the
script name as a variable to the input script. In the latter case, a
[variable](variable) called \"fname\" could be used in place of SELF,
e.g.

``` bash
lmp_g++ -var fname in.script < in.script
```

The second argument to the jump command is optional. If specified, it is
treated as a label and the new file is scanned (without executing
commands) until the label is found, and commands are executed from that
point forward. This can be used to loop over a portion of the input
script, as in this example. These commands perform 10 runs, each of
10000 steps, and create 10 dump files named file.1, file.2, etc. The
[next](next) command is used to exit the loop after 10 iterations. When
the \"a\" variable has been incremented for the tenth time, it will
cause the next jump command to be skipped.

``` LAMMPS
variable a loop 10
label loop
dump 1 all atom 100 file.$a
run 10000
undump 1
next a
jump in.lj loop
```

If the jump *file* argument is a variable, the jump command can be used
to cause different processor partitions to run different input scripts.
In this example, LAMMPS is run on 40 processors, with 4 partitions of 10
procs each. An in.file containing the example variable and jump command
will cause each partition to run a different simulation.

``` LAMMPS
mpirun -np 40 lmp_ibm -partition 4x10 -in in.file

variable f world script.1 script.2 script.3 script.4
jump $f
```

Here is an example of a loop which checks every 1000 steps if the system
temperature has reached a certain value, and if so, breaks out of the
loop to finish the run. Note that any variable could be checked, so long
as it is current on the timestep when the run completes. As explained on
the [variable](variable) doc page, this can be ensured by including the
variable in thermodynamic output.

``` LAMMPS
variable myTemp equal temp
label loop
variable a loop 1000
run 1000
if "${myTemp} < 300.0" then "jump SELF break"
next a
jump SELF loop
label break
print "ALL DONE"
```

Here is an example of a double loop which uses the if and [jump](jump)
commands to break out of the inner loop when a condition is met, then
continues iterating through the outer loop.

``` LAMMPS
label       loopa
variable    a loop 5
  label     loopb
  variable  b loop 5
    print     "A,B = $a,$b"
    run       10000
    if        "$b > 2" then "jump SELF break"
  next      b
  jump      in.script loopb
  label       break
  variable    b delete
next        a
jump        SELF loopa
```

## Restrictions

If you jump to a file and it does not contain the specified label,
LAMMPS will come to the end of the file and exit.

## Related commands

[variable](variable), [include](include), [label](label), [next](next)

## Default

none
