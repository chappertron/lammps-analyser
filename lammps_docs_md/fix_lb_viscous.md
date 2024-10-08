# fix lb/viscous command

## Syntax

    fix ID group-ID lb/viscous

-   ID, group-ID are documented in [fix](fix) command
-   lb/viscous = style name of this fix command

## Examples

``` LAMMPS
fix 1 flow lb/viscous
```

## Description

This fix is similar to the [fix viscous](fix_viscous) command, and is to
be used in place of that command when a lattice-Boltzmann fluid is
present using the [fix lb/fluid](fix_lb_fluid). This should be used in
conjunction with one of the built-in LAMMPS integrators, such as [fix
NVE](fix_nve) or [fix rigid](fix_rigid).

This fix adds a viscous force to each atom to cause it move with the
same velocity as the fluid (an equal and opposite force is applied to
the fluid via [fix lb/fluid](fix_lb_fluid)). When [fix
lb/fluid](fix_lb_fluid) is called with the noise option, the atoms will
also experience random forces which will thermalize them to the same
temperature as the fluid. In this way, the combination of this fix with
[fix lb/fluid](fix_lb_fluid) and a LAMMPS integrator like [fix
NVE](fix_nve) is analogous to [fix langevin](fix_langevin) except here
the fluid is explicit. The temperature of the particles can be monitored
via the scalar output of [fix lb/fluid](fix_lb_fluid).

------------------------------------------------------------------------

For details of this fix, as well as descriptions and results of several
test runs, see [Denniston et al.](fluid-Denniston2). Please include a
citation to this paper if this fix is used in work contributing to
published research.

------------------------------------------------------------------------

## Restart, fix_modify, output, run start/stop, minimize info

As described in the [fix viscous](fix_viscous) documentation:

\"No information about this fix is written to [binary restart
files](restart). None of the [fix_modify](fix_modify) options are
relevant to this fix. No global or per-atom quantities are stored by
this fix for access by various [output commands](Howto_output). No
parameter of this fix can be used with the *start/stop* keywords of the
[run](run) command.

The forces due to this fix are imposed during an energy minimization,
invoked by the [minimize](minimize) command. This fix should only be
used with damped dynamics minimizers that allow for non-conservative
forces. See the [min_style](min_style) command for details.\"

## Restrictions

This fix is part of the LATBOLTZ package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

Can only be used if a lattice-Boltzmann fluid has been created via the
[fix lb/fluid](fix_lb_fluid) command, and must come after this command.

## Related commands

[fix lb/fluid](fix_lb_fluid)

## Default

none

------------------------------------------------------------------------

::: {#fluid-Denniston2}
**(Denniston et al.)** Denniston, C., Afrasiabian, N., Cole-Andre, M.G.,
Mackay, F. E., Ollila, S.T.T., and Whitehead, T., LAMMPS lb/fluid fix
version 2: Improved Hydrodynamic Forces Implemented into LAMMPS through
a lattice-Boltzmann fluid, Computer Physics Communications 275 (2022)
[108318](https://doi.org/10.1016/j.cpc.2022.108318)\_ .
:::
