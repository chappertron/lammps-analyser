# fix nvt/sllod command

Accelerator Variants: *nvt/sllod/intel*, *nvt/sllod/omp*, *nvt/sllod/kk*

## Syntax

    fix ID group-ID nvt/sllod keyword value ...

-   ID, group-ID are documented in [fix](fix) command

-   nvt/sllod = style name of this fix command

-   zero or more keyword/value pairs may be appended

        keyword = *psllod*
          *psllod* value = *no* or *yes* = use SLLOD or p-SLLOD variant, respectively

-   additional thermostat related keyword/value pairs from the [fix
    nvt](fix_nh) command can be appended

## Examples

``` LAMMPS
fix 1 all nvt/sllod temp 300.0 300.0 100.0
fix 1 all nvt/sllod temp 300.0 300.0 100.0 drag 0.2
```

## Description

Perform constant NVT integration to update positions and velocities each
timestep for atoms in the group using a Nose/Hoover temperature
thermostat. V is volume; T is temperature. This creates a system
trajectory consistent with the canonical ensemble.

This thermostat is used for a simulation box that is changing size
and/or shape, for example in a non-equilibrium MD (NEMD) simulation. The
size/shape change is induced by use of the [fix deform](fix_deform)
command, so each point in the simulation box can be thought of as having
a \"streaming\" velocity. This position-dependent streaming velocity is
subtracted from each atom\'s actual velocity to yield a thermal velocity
which is used for temperature computation and thermostatting. For
example, if the box is being sheared in x, relative to y, then points at
the bottom of the box (low y) have a small x velocity, while points at
the top of the box (hi y) have a large x velocity. These velocities do
not contribute to the thermal \"temperature\" of the atom.

:::: note
::: title
Note
:::

[Fix deform](fix_deform) has an option for remapping either atom
coordinates or velocities to the changing simulation box. To use fix
nvt/sllod, fix deform should NOT remap atom positions, because fix
nvt/sllod adjusts the atom positions and velocities to create a velocity
profile that matches the changing box size/shape. Fix deform SHOULD
remap atom velocities when atoms cross periodic boundaries since that is
consistent with maintaining the velocity profile created by fix
nvt/sllod. LAMMPS will give an error if this setting is not consistent.
::::

The SLLOD equations of motion, originally proposed by Hoover and Ladd
(see [(Evans and Morriss)](Evans3)), were proven to be equivalent to
Newton\'s equations of motion for shear flow by [(Evans and
Morriss)](Evans3). They were later shown to generate the desired
velocity gradient and the correct production of work by stresses for all
forms of homogeneous flow by [(Daivis and Todd)](Daivis).

::: versionchanged
8Feb2023
:::

For the default (*psllod* = *no*), the LAMMPS implementation adheres to
the standard SLLOD equations of motion, as defined by [(Evans and
Morriss)](Evans3). The option *psllod* = *yes* invokes the slightly
different SLLOD variant first introduced by [(Tuckerman et
al.)](Tuckerman) as g-SLLOD and later by [(Edwards)](Edwards) as
p-SLLOD. In all cases, the equations of motion are coupled to a
Nose/Hoover chain thermostat in a velocity Verlet formulation, closely
following the implementation used for the [fix nvt](fix_nh) command.

:::: note
::: title
Note
:::

A recent (2017) book by [(Todd and Daivis)](Todd-sllod) discusses use of
the SLLOD method and non-equilibrium MD (NEMD) thermostatting generally,
for both simple and complex fluids, e.g. molecular systems. The latter
can be tricky to do correctly.
::::

Additional parameters affecting the thermostat are specified by keywords
and values documented with the [fix nvt](fix_nh) command. See, for
example, discussion of the *temp* and *drag* keywords.

This fix computes a temperature each timestep. To do this, the fix
creates its own compute of style \"temp/deform\", as if this command had
been issued:

``` LAMMPS
compute fix-ID_temp group-ID temp/deform
```

See the [compute temp/deform](compute_temp_deform) command for details.
Note that the ID of the new compute is the fix-ID + underscore +
\"temp\", and the group for the new compute is the same as the fix
group.

Note that this is NOT the compute used by thermodynamic output (see the
[thermo_style](thermo_style) command) with ID = *thermo_temp*. This
means you can change the attributes of this fix\'s temperature (e.g. its
degrees-of-freedom) via the [compute_modify](compute_modify) command or
print this temperature during thermodynamic output via the [thermo_style
custom](thermo_style) command using the appropriate compute-ID. It also
means that changing attributes of *thermo_temp* will have no effect on
this fix.

Like other fixes that perform thermostatting, this fix can be used with
[compute commands](compute) that remove a \"bias\" from the atom
velocities. E.g. to apply the thermostat only to atoms within a spatial
[region](region), or to remove the center-of-mass velocity from a group
of atoms, or to remove the x-component of velocity from the calculation.

This is not done by default, but only if the [fix_modify](fix_modify)
command is used to assign a temperature compute to this fix that
includes such a bias term. See the doc pages for individual [compute
temp commands](compute) to determine which ones include a bias. In this
case, the thermostat works in the following manner: bias is removed from
each atom, thermostatting is performed on the remaining thermal degrees
of freedom, and the bias is added back in.

------------------------------------------------------------------------

Styles with a *gpu*, *intel*, *kk*, *omp*, or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the [Accelerator packages](Speed_packages)
page. The accelerated styles take the same arguments and should produce
the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, INTEL, KOKKOS, OPENMP, and
OPT packages, respectively. They are only enabled if LAMMPS was built
with those packages. See the [Build package](Build_package) page for
more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the [-suffix command-line
switch](Run_options) when you invoke LAMMPS, or you can use the
[suffix](suffix) command in your input script.

See the [Accelerator packages](Speed_packages) page for more
instructions on how to use the accelerated styles effectively.

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
as does the [fix nvt](fix_nh) command.

This fix can ramp its target temperature over multiple runs, using the
*start* and *stop* keywords of the [run](run) command. See the
[run](run) command for details of how to do this.

This fix is not invoked during [energy minimization](minimize).

## Restrictions

This fix works best without Nose-Hoover chain thermostats, i.e. using
*tchain* = 1. Setting *tchain* to larger values can result in poor
equilibration.

## Related commands

[fix nve](fix_nve), [fix nvt](fix_nh), [fix
temp/rescale](fix_temp_rescale), [fix langevin](fix_langevin),
[fix_modify](fix_modify), [compute temp/deform](compute_temp_deform)

## Default

Same as [fix nvt](fix_nh), except *tchain* = 1, psllod = *no*.

------------------------------------------------------------------------

::: {#Evans3}
**(Evans and Morriss)** Evans and Morriss, Phys Rev A, 30, 1528 (1984).
:::

::: {#Daivis}
**(Daivis and Todd)** Daivis and Todd, J Chem Phys, 124, 194103 (2006).
:::

::: {#Todd-sllod}
**(Todd and Daivis)** Todd and Daivis, Nonequilibrium Molecular Dynamics
(book), Cambridge University Press, (2017)
<https://doi.org/10.1017/9781139017848>.
:::

::: {#Tuckerman}
**(Tuckerman et al.)** Tuckerman, Mundy, Balasubramanian, and Klein, J
Chem Phys 106, 5615 (1997).
:::

::: {#Edwards}
**(Edwards)** Edwards, Baig, and Keffer, J Chem Phys 124, 194104 (2006).
:::
