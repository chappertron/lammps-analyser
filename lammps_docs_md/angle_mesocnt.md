# angle_style mesocnt command

## Syntax

``` LAMMPS
angle_style mesocnt
```

## Examples

``` LAMMPS
angle_style mesocnt
angle_coeff 1 buckling C 10 10 20.0
angle_coeff 4 harmonic C 8 4 10.0
angle_coeff 2 buckling custom 400.0 50.0 5.0
angle_coeff 1 harmonic custom 300.0
```

## Description

::: versionadded
15Sep2022
:::

The *mesocnt* angle style uses the potential

$$\begin{aligned}
E = K_\text{H} \Delta \theta^2, \qquad |\Delta \theta| < \Delta
\theta_\text{B} \\
E = K_\text{H} \Delta \theta_\text{B}^2 +
K_\text{B} (\Delta \theta - \Delta \theta_\text{B}), \qquad |\Delta
\theta| \geq \Delta \theta_\text{B}
\end{aligned}$$

where $\Delta \theta = \theta - \pi$ is the bending angle of the
nanotube, $K_\text{H}$ and $K_\text{B}$ are prefactors for the harmonic
and linear regime respectively and $\Delta
\theta_\text{B}$ is the buckling angle. Note that the usual 1/2 factor
for the harmonic potential is included in $K_\text{H}$.

The style implements parameterization presets of $K_\text{H}$,
$K_\text{B}$ and $\Delta \theta_\text{B}$ for mesoscopic simulations of
carbon nanotubes based on the atomistic simulations of
[(Srivastava)](Srivastava_2) and buckling considerations of
[(Zhigilei)](Zhigilei1_1).

The following coefficients must be defined for each angle type via the
[angle_coeff](angle_coeff) command as in the examples above, or in the
data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   mode = *buckling* or *harmonic*
-   preset = *C* or *custom*
-   additional parameters depending on preset

If mode *harmonic* is chosen, the potential is simply harmonic and does
not switch to the linear term when the buckling angle is reached. In
*buckling* mode, the full piecewise potential is used.

Preset *C* is for carbon nanotubes, and the additional parameters are:

-   chiral index $n$ (unitless)
-   chiral index $m$ (unitless)
-   $r_0$ (distance)

Here, $r_0$ is the equilibrium distance of the bonds included in the
angle, see [bond_style mesocnt](bond_mesocnt).

In harmonic mode with preset *custom*, the additional parameter is:

-   $K_\text{H}$ (energy)

Hence, this setting is simply a wrapper for [bond_style
harmonic](bond_harmonic) with an equilibrium angle of 180 degrees.

In harmonic mode with preset *custom*, the additional parameters are:

-   $K_\text{H}$ (energy)
-   $K_\text{B}$ (energy)
-   $\Delta \theta_\text{B}$ (degrees)

$\Delta \theta_\text{B}$ is specified in degrees, but LAMMPS converts it
to radians internally; hence $K_\text{H}$ is effectively energy per
radian\^2 and $K_\text{B}$ is energy per radian.

------------------------------------------------------------------------

In *buckling* mode, this angle style adds the *buckled* property to all
atoms in the simulation, which is an integer flag indicating whether the
bending angle at a given atom has exceeded $\Delta
\theta_\text{B}$. It can be accessed as an atomic variable, e.g. for
custom dump commands, as *i_buckled*.

:::: note
::: title
Note
:::

If the initial state of the simulation contains buckled nanotubes and
[pair_style mesocnt](pair_mesocnt) is used, the *i_buckled* atomic
variable needs to be initialized before the pair_style is defined by
doing a *run 0* command straight after the angle_style command. See
below for an example.
::::

If CNTs are already buckled at the start of the simulation, this script
will correctly initialize *i_buckled*:

``` LAMMPS
angle_style mesocnt
angle_coeff 1 buckling C 10 10 20.0

run 0

pair_style mesocnt 60.0
pair_coeff * * C_10_10.mesocnt 1
```

## Restrictions

This angle style can only be used if LAMMPS was built with the MOLECULE
and MESONT packages. See the [Build package](Build_package) doc page for
more info.

## Related commands

[angle_coeff](angle_coeff)

## Default

none

------------------------------------------------------------------------

::: {#Srivastava_2}
**(Srivastava)** Zhigilei, Wei, Srivastava, Phys. Rev. B 71, 165417
(2005).
:::

::: {#Zhigilei1_1}
**(Zhigilei)** Volkov and Zhigilei, ACS Nano 4, 6187 (2010).
:::
