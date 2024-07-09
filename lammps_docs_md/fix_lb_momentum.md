# fix lb/momentum command

## Syntax

    fix ID group-ID lb/momentum nevery keyword values ...

-   ID, group-ID are documented in the [fix](fix) command

-   lb/momentum = style name of this fix command

-   nevery = adjust the momentum every this many timesteps

-   zero or more keyword/value pairs may be appended

-   keyword = *linear*

        *linear* values = xflag yflag zflag
          xflag,yflag,zflag = 0/1 to exclude/include each dimension.

## Examples

``` LAMMPS
fix 1 sphere lb/momentum
fix 1 all lb/momentum linear 1 1 0
```

## Description

This fix is based on the [fix momentum](fix_momentum) command, and was
created to be used in place of that command, when a lattice-Boltzmann
fluid is present.

Zero the total linear momentum of the system, including both the atoms
specified by group-ID and the lattice-Boltzmann fluid every nevery
timesteps. If there are no atoms specified by group-ID only the fluid
momentum is affected. This is accomplished by adjusting the particle
velocities and the fluid velocities at each lattice site.

:::: note
::: title
Note
:::

This fix only considers the linear momentum of the system.
::::

By default, the subtraction is performed for each dimension. This can be
changed by specifying the keyword *linear*, along with a set of three
flags set to 0/1 in order to exclude/ include the corresponding
dimension.

## Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to [binary restart
files](restart). None of the [fix_modify](fix_modify) options are
relevant to this fix. No global or per-atom quantities are stored by
this fix for access by various [output commands](Howto_output). No
parameter of this fix can be used with the *start/stop* keywords of the
[run](run) command. This fix is not invoked during [energy
minimization](minimize).

## Restrictions

Can only be used if a lattice-Boltzmann fluid has been created via the
[fix lb/fluid](fix_lb_fluid) command, and must come after this command.

This fix is part of the LATBOLTZ package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

## Related commands

[fix momentum](fix_momentum), [fix lb/fluid](fix_lb_fluid)

## Default

Zeros the total system linear momentum in each dimension.
