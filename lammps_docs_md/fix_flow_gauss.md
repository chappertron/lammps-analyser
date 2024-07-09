# fix flow/gauss command

## Syntax

    fix ID group-ID flow/gauss xflag yflag zflag keyword

-   ID, group-ID are documented in [fix](fix) command

-   flow/gauss = style name of this fix command

-   xflag,yflag,zflag = 0 or 1

        0 = do not conserve current in this dimension
        1 = conserve current in this dimension

-   zero or more keyword/value pairs may be appended

-   keyword = *energy*

        *energy* value = *no* or *yes*
          *no* = do not compute work done by this fix
          *yes* = compute work done by this fix

## Examples

``` LAMMPS
fix GD fluid flow/gauss 1 0 0
fix GD fluid flow/gauss 1 1 1 energy yes
```

## Description

This fix implements the Gaussian dynamics (GD) method to simulate a
system at constant mass flux [(Strong)](Strong). GD is a nonequilibrium
molecular dynamics simulation method that can be used to study fluid
flows through pores, pipes, and channels. In its original implementation
GD was used to compute the pressure required to achieve a fixed mass
flux through an opening. The flux can be conserved in any combination of
the directions, x, y, or z, using xflag,yflag,zflag. This fix does not
initialize a net flux through a system, it only conserves the
center-of-mass momentum that is present when the fix is declared in the
input script. Use the [velocity](velocity) command to generate an
initial center-of-mass momentum.

GD applies an external fluctuating gravitational field that acts as a
driving force to keep the system away from equilibrium. To maintain
steady state, a profile-unbiased thermostat must be implemented to
dissipate the heat that is added by the driving force. [Compute
temp/profile](compute_temp_profile) can be used to implement a
profile-unbiased thermostat.

A common use of this fix is to compute a pressure drop across a pipe,
pore, or membrane. The pressure profile can be computed in LAMMPS with
[compute stress/atom](compute_stress_atom) and [fix
ave/chunk](fix_ave_chunk), or with the hardy method in [fix
atc](fix_atc). Note that the simple [compute
stress/atom](compute_stress_atom) method is only accurate away from
inhomogeneities in the fluid, such as fixed wall atoms. Further, the
computed pressure profile must be corrected for the acceleration applied
by GD before computing a pressure drop or comparing it to other methods,
such as the pump method [(Zhu)](Zhu). The pressure correction is
discussed and described in [(Strong)](Strong).

For a complete example including the considerations discussed above, see
the examples/PACKAGES/flow_gauss directory.

:::: note
::: title
Note
:::

Only the flux of the atoms in group-ID will be conserved. If the
velocities of the group-ID atoms are coupled to the velocities of other
atoms in the simulation, the flux will not be conserved. For example, in
a simulation with fluid atoms and harmonically constrained wall atoms,
if a single thermostat is applied to group *all*, the fluid atom
velocities will be coupled to the wall atom velocities, and the flux
will not be conserved. This issue can be avoided by thermostatting the
fluid and wall groups separately.
::::

Adding an acceleration to atoms does work on the system. This added
energy can be optionally subtracted from the potential energy for the
thermodynamic output (see below) to check that the timestep is small
enough to conserve energy. Since the applied acceleration is fluctuating
in time, the work cannot be computed from a potential. As a result,
computing the work is slightly more computationally expensive than
usual, so it is not performed by default. To invoke the work
calculation, use the *energy* keyword. The [fix_modify](fix_modify)
*energy* option also invokes the work calculation, and overrides an
*energy no* setting here. If neither *energy yes* or *fix_modify energy
yes* are set, the global scalar computed by the fix will return zero.

:::: note
::: title
Note
:::

In order to check energy conservation, any other fixes that do work on
the system must have *fix_modify energy yes* set as well. This includes
thermostat fixes and any constraints that hold the positions of wall
atoms fixed, such as [fix spring/self](fix_spring_self).
::::

If this fix is used in a simulation with the [rRESPA](run_style)
integrator, the applied acceleration must be computed and applied at the
same rRESPA level as the interactions between the flowing fluid and the
obstacle. The rRESPA level at which the acceleration is applied can be
changed using the [fix_modify](fix_modify) *respa* option discussed
below. If the flowing fluid and the obstacle interact through multiple
interactions that are computed at different rRESPA levels, then there
must be a separate flow/gauss fix for each level. For example, if the
flowing fluid and obstacle interact through pairwise and long-range
Coulomb interactions, which are computed at rRESPA levels 3 and 4,
respectively, then there must be two separate flow/gauss fixes, one that
specifies *fix_modify respa 3* and one with *fix_modify respa 4*.

------------------------------------------------------------------------

## Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to [binary restart
files](restart).

The [fix_modify](fix_modify) *energy* option is supported by this fix to
add the potential energy added by the fix to the global potential energy
of the system as part of [thermodynamic output](thermo_style). The
default setting for this fix is [fix_modify energy no](fix_modify).

The [fix_modify](fix_modify) *respa* option is supported by this fix.
This allows the user to set at which level of the [rRESPA](run_style)
integrator the fix computes and adds the external acceleration. Default
is the outermost level.

This fix computes a global scalar and a global 3-vector of forces, which
can be accessed by various [output commands](Howto_output). The scalar
is the negative of the work done on the system, see the discussion
above. It is only calculated if the *energy* keyword is enabled or
[fix_modify energy yes](fix_modify) is set.

The vector is the total force that this fix applied to the group of
atoms on the current timestep. The scalar and vector values calculated
by this fix are \"extensive\".

No parameter of this fix can be used with the *start/stop* keywords of
the [run](run) command.

This fix is not invoked during [energy minimization](minimize).

## Restrictions

This fix is part of the EXTRA-FIX package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

## Related commands

[fix addforce](fix_addforce), [compute
temp/profile](compute_temp_profile), [velocity](velocity)

## Default

The option default for the *energy* keyword is energy = no.

------------------------------------------------------------------------

::: {#Strong}
**(Strong)** Strong and Eaves, J. Phys. Chem. B 121, 189 (2017).
:::

::: {#Evans2}
**(Evans)** Evans and Morriss, Phys. Rev. Lett. 56, 2172 (1986).
:::

::: {#Zhu}
**(Zhu)** Zhu, Tajkhorshid, and Schulten, Biophys. J. 83, 154 (2002).
:::
