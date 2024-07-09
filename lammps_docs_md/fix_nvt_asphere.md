# fix nvt/asphere command

Accelerator Variants: *nvt/asphere/omp*

## Syntax

    fix ID group-ID nvt/asphere keyword value ...

-   ID, group-ID are documented in [fix](fix) command
-   nvt/asphere = style name of this fix command
-   additional thermostat related keyword/value pairs from the [fix
    nvt](fix_nh) command can be appended

## Examples

``` LAMMPS
fix 1 all nvt/asphere temp 300.0 300.0 100.0
fix 1 all nvt/asphere temp 300.0 300.0 100.0 drag 0.2
```

## Description

Perform constant NVT integration to update position, velocity,
orientation, and angular velocity each timestep for aspherical or
ellipsoidal particles in the group using a Nose/Hoover temperature
thermostat. V is volume; T is temperature. This creates a system
trajectory consistent with the canonical ensemble.

This fix differs from the [fix nvt](fix_nh) command, which assumes point
particles and only updates their position and velocity.

The thermostat is applied to both the translational and rotational
degrees of freedom for the aspherical particles, assuming a compute is
used which calculates a temperature that includes the rotational degrees
of freedom (see below). The translational degrees of freedom can also
have a bias velocity removed from them before thermostatting takes
place; see the description below.

Additional parameters affecting the thermostat are specified by keywords
and values documented with the [fix nvt](fix_nh) command. See, for
example, discussion of the *temp* and *drag* keywords.

This fix computes a temperature each timestep. To do this, the fix
creates its own compute of style \"temp/asphere\", as if this command
had been issued:

``` LAMMPS
compute fix-ID_temp group-ID temp/asphere
```

See the [compute temp/asphere](compute_temp_asphere) command for
details. Note that the ID of the new compute is the fix-ID +
underscore + \"temp\", and the group for the new compute is the same as
the fix group.

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

This fix is part of the ASPHERE package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

This fix requires that atoms store torque and angular momentum and a
quaternion as defined by the [atom_style ellipsoid](atom_style) command.

All particles in the group must be finite-size. They cannot be point
particles, but they can be aspherical or spherical as defined by their
shape attribute.

## Related commands

[fix nvt](fix_nh), [fix nve_asphere](fix_nve_asphere), [fix
npt_asphere](fix_npt_asphere), [fix_modify](fix_modify)

## Default

none
