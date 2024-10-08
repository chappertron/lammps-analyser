# fix nvt/sllod/eff command

## Syntax

    fix ID group-ID nvt/sllod/eff keyword value ...

-   ID, group-ID are documented in [fix](fix) command

-   nvt/sllod/eff = style name of this fix command

-   zero or more keyword/value pairs may be appended

        keyword = *psllod*
          *psllod* value = *no* or *yes* = use SLLOD or p-SLLOD variant, respectively

-   additional thermostat related keyword/value pairs from the [fix
    nvt/eff](fix_nh_eff) command may be appended, too.

## Examples

``` LAMMPS
fix 1 all nvt/sllod/eff temp 300.0 300.0 0.1
fix 1 all nvt/sllod/eff temp 300.0 300.0 0.1 drag 0.2
```

## Description

Perform constant NVT integration to update positions and velocities each
timestep for nuclei and electrons in the group for the [electron force
field](pair_eff) model, using a Nose/Hoover temperature thermostat. V is
volume; T is temperature. This creates a system trajectory consistent
with the canonical ensemble.

The operation of this fix is exactly like that described by the [fix
nvt/sllod](fix_nvt_sllod) command, except that the radius and radial
velocity of electrons are also updated and thermostatted. Likewise the
temperature calculated by the fix, using the compute it creates (as
discussed in the [fix nvt, npt, and nph](fix_nh) doc page), is performed
with a [compute temp/deform/eff](compute_temp_deform_eff) command that
includes the eFF contribution to the temperature from the electron
radial velocity.

## Restart, fix_modify, output, run start/stop, minimize info

This fix writes the state of the Nose/Hoover thermostat to [binary
restart files](restart). See the [read_restart](read_restart) command
for info on how to re-specify a fix in an input script that reads a
restart file, so that the operation of the fix continues in an
uninterrupted fashion.

The [fix_modify](fix_modify) *temp* option is supported by this fix. You
can use it to assign a [compute](compute) you have defined to this fix
which will be used in its thermostatting procedure.

The cumulative energy change in the system imposed by this fix is
included in the [thermodynamic output](thermo_style) keywords *ecouple*
and *econserve*. See the [thermo_style](thermo_style) doc page for
details.

This fix computes the same global scalar and global vector of quantities
as does the [fix nvt/eff](fix_nh_eff) command.

This fix can ramp its target temperature over multiple runs, using the
*start* and *stop* keywords of the [run](run) command. See the
[run](run) command for details of how to do this.

This fix is not invoked during [energy minimization](minimize).

## Restrictions

This fix is part of the EFF package. It is only enabled if LAMMPS was
built with that package. See the [Build package](Build_package) page for
more info.

This fix works best without Nose-Hoover chain thermostats, i.e. using
tchain = 1. Setting tchain to larger values can result in poor
equilibration.

## Related commands

[fix nve/eff](fix_nve_eff), [fix nvt/eff](fix_nh_eff), [fix
langevin/eff](fix_langevin_eff), [fix nvt/sllod](fix_nvt_sllod),
[fix_modify](fix_modify), [compute
temp/deform/eff](compute_temp_deform_eff)

## Default

Same as [fix nvt/eff](fix_nh_eff), except tchain = 1.

------------------------------------------------------------------------

::: {#Tuckerman2}
**(Tuckerman)** Tuckerman, Mundy, Balasubramanian, Klein, J Chem Phys,
106, 5615 (1997).
:::
