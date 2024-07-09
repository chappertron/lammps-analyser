# fix drag command

## Syntax

    fix ID group-ID drag x y z fmag delta

-   ID, group-ID are documented in [fix](fix) command
-   drag = style name of this fix command
-   x,y,z = coord to drag atoms towards
-   fmag = magnitude of force to apply to each atom (force units)
-   delta = cutoff distance inside of which force is not applied
    (distance units)

## Examples

``` LAMMPS
fix center small-molecule drag 0.0 10.0 0.0 5.0 2.0
```

## Description

Apply a force to each atom in a group to drag it towards the point
(x,y,z). The magnitude of the force is specified by fmag. If an atom is
closer than a distance delta to the point, then the force is not
applied.

Any of the x,y,z values can be specified as NULL which means do not
include that dimension in the distance calculation or force application.

This command can be used to steer one or more atoms to a new location in
the simulation.

## Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to [binary restart
files](restart).

The [fix_modify](fix_modify) *respa* option is supported by this fix.
This allows to set at which level of the [r-RESPA](run_style) integrator
the fix is adding its forces. Default is the outermost level.

This fix computes a global 3-vector of forces, which can be accessed by
various [output commands](Howto_output). This is the total force on the
group of atoms by the drag force. The vector values calculated by this
fix are \"extensive\".

No parameter of this fix can be used with the *start/stop* keywords of
the [run](run) command. This fix is not invoked during [energy
minimization](minimize).

## Restrictions

This fix is part of the EXTRA-FIX package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

## Related commands

[fix spring](fix_spring), [fix spring/self](fix_spring_self), [fix
spring/rg](fix_spring_rg), [fix smd](fix_smd)

## Default

none
