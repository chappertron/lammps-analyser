# fix nvk command

## Syntax

    fix ID group-ID nvk

-   ID, group-ID are documented in [fix](fix) command
-   nvk = style name of this fix command

## Examples

``` LAMMPS
fix 1 all nvk
```

## Description

Perform constant kinetic energy integration using the Gaussian
thermostat to update position and velocity for atoms in the group each
timestep. V is volume; K is kinetic energy. This creates a system
trajectory consistent with the isokinetic ensemble.

The equations of motion used are those of Minary et al in
[(Minary)](nvk-Minary), a variant of those initially given by Zhang in
[(Zhang)](nvk-Zhang).

The kinetic energy will be held constant at its value given when fix nvk
is initiated. If a different kinetic energy is desired, the
[velocity](velocity) command should be used to change the kinetic energy
prior to this fix.

------------------------------------------------------------------------

## Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to [binary restart
files](restart). None of the [fix_modify](fix_modify) options are
relevant to this fix. No global or per-atom quantities are stored by
this fix for access by various [output commands](Howto_output). No
parameter of this fix can be used with the *start/stop* keywords of the
[run](run) command. This fix is not invoked during [energy
minimization](minimize).

## Restrictions

The Gaussian thermostat only works when it is applied to all atoms in
the simulation box. Therefore, the group must be set to all.

This fix has not yet been implemented to work with the RESPA integrator.

This fix is part of the EXTRA-FIX package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

## Related commands

none

## Default

none

------------------------------------------------------------------------

::: {#nvk-Minary}
**(Minary)** Minary, Martyna, and Tuckerman, J Chem Phys, 18, 2510
(2003).
:::

::: {#nvk-Zhang}
**(Zhang)** Zhang, J Chem Phys, 106, 6102 (1997).
:::
