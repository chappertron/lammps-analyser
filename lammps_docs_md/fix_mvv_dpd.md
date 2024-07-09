# fix mvv/dpd command

# fix mvv/edpd command

# fix mvv/tdpd command

## Syntax

    fix ID group-ID mvv/dpd lambda

    fix ID group-ID mvv/edpd lambda

    fix ID group-ID mvv/tdpd lambda

-   ID, group-ID are documented in [fix](fix) command
-   mvv/dpd, mvv/edpd, mvv/tdpd = style name of this fix command
-   lambda = (optional) relaxation parameter (unitless)

## Examples

``` LAMMPS
fix 1 all mvv/dpd
fix 1 all mvv/dpd 0.5
fix 1 all mvv/edpd
fix 1 all mvv/edpd 0.5
fix 1 all mvv/tdpd
fix 1 all mvv/tdpd 0.5
```

## Description

Perform time integration using the modified velocity-Verlet (MVV)
algorithm to update position and velocity (fix mvv/dpd), or position,
velocity and temperature (fix mvv/edpd), or position, velocity and
concentration (fix mvv/tdpd) for particles in the group each timestep.

The modified velocity-Verlet (MVV) algorithm aims to improve the
stability of the time integrator by using an extrapolated version of the
velocity for the force evaluation:

$$\begin{aligned}
v(t+\frac{\Delta t}{2}) = & v(t) + \frac{\Delta t}{2}\cdot a(t) \\
r(t+\Delta t) = & r(t) + \Delta t\cdot v(t+\frac{\Delta t}{2}) \\
a(t+\Delta t) = & \frac{1}{m}\cdot F\left[ r(t+\Delta t), v(t) +\lambda \cdot \Delta t\cdot a(t)\right] \\
v(t+\Delta t) = & v(t+\frac{\Delta t}{2}) + \frac{\Delta t}{2}\cdot a(t+\Delta t)
\end{aligned}$$

where the parameter $\lambda$ depends on the specific choice of DPD
parameters, and needs to be tuned on a case-by-case basis. Specification
of a *lambda* value is optional. If specified, the setting must be from
0.0 to 1.0. If not specified, a default value of 0.5 is used, which
effectively reproduces the standard velocity-Verlet (VV) scheme. For
more details, see [Groot](Groot2).

Fix *mvv/dpd* updates the position and velocity of each atom. It can be
used with the [pair_style mdpd](pair_mesodpd) command or other pair
styles such as [pair dpd](pair_dpd).

Fix *mvv/edpd* updates the per-atom temperature, in addition to position
and velocity, and must be used with the [pair_style edpd](pair_mesodpd)
command.

Fix *mvv/tdpd* updates the per-atom chemical concentration, in addition
to position and velocity, and must be used with the [pair_style
tdpd](pair_mesodpd) command.

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

These fixes are part of the DPD-MESO package. They are only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

## Related commands

[pair_style mdpd](pair_mesodpd), [pair_style edpd](pair_mesodpd),
[pair_style tdpd](pair_mesodpd)

## Default

The default value for the optional *lambda* parameter is 0.5.

------------------------------------------------------------------------

::: {#Groot2}
**(Groot)** Groot and Warren, J Chem Phys, 107: 4423-4435 (1997). DOI:
10.1063/1.474784
:::
