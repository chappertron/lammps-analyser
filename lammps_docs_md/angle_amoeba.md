# angle_style amoeba command

## Syntax

``` LAMMPS
angle_style amoeba
```

## Examples

``` LAMMPS
angle_style amoeba
angle_coeff * 75.0 -25.0 1.0 0.3 0.02 0.003
angle_coeff * ba 3.6551 24.895 1.0119 1.5228
angle_coeff * ub -7.6 1.5537
```

## Description

The *amoeba* angle style uses the potential

$$\begin{aligned}
E & = E_a + E_{ba} + E_{ub} \\
E_a & = K_2\left(\theta - \theta_0\right)^2 + K_3\left(\theta - \theta_0\right)^3 + K_4\left(\theta - \theta_0\right)^4 + K_5\left(\theta - \theta_0\right)^5 + K_6\left(\theta - \theta_0\right)^6 \\
E_{ba} & = N_1 (r_{ij} - r_1) (\theta - \theta_0) + N_2(r_{jk} - r_2)(\theta - \theta_0) \\
E_{UB} & = K_{ub} (r_{ik} - r_{ub})^2
\end{aligned}$$

where $E_a$ is the angle term, $E_{ba}$ is a bond-angle term, $E_{UB}$
is a Urey-Bradley bond term, $\theta_0$ is the equilibrium angle, $r_1$
and $r_2$ are the equilibrium bond lengths, and $r_{ub}$ is the
equilibrium Urey-Bradley bond length.

These formulas match how the Tinker MD code performs its angle
calculations for the AMOEBA and HIPPO force fields. See the [Howto
amoeba](Howto_amoeba) page for more information about the implementation
of AMOEBA and HIPPO in LAMMPS.

Note that the $E_a$ and $E_{ba}$ formulas are identical to those used
for the [angle_style class2/p6](angle_class2) command, however there is
no bond-bond cross term formula for $E_{bb}$. Additionally, there is a
$E_{UB}$ term for a Urey-Bradley bond. It is effectively a harmonic bond
between the I and K atoms of angle IJK, even though that bond is not
enumerated in the \"Bonds\" section of the data file.

There are also two ways that Tinker computes the angle $\theta$ in the
$E_a$ formula. The first is the standard way of treating IJK as an
\"in-plane\" angle. The second is an \"out-of-plane\" method which
Tinker may use if the center atom J in the angle is bonded to one
additional atom in addition to I and K. In this case, all 4 atoms are
used to compute the $E_a$ formula, resulting in forces on all 4 atoms.
In the Tinker PRM file, these 2 options are denoted by *angle* versus
*anglep* entries in the \"Angle Bending Parameters\" section of the PRM
force field file. The *pflag* coefficient described below selects
between the 2 options.

------------------------------------------------------------------------

Coefficients for the $E_a$, $E_{bb}$, and $E_{ub}$ formulas must be
defined for each angle type via the [angle_coeff](angle_coeff) command
as in the example above, or in the data file or restart files read by
the [read_data](read_data) or [read_restart](read_restart) commands.

These are the 8 coefficients for the $E_a$ formula:

-   pflag = 0 or 1
-   ubflag = 0 or 1
-   $\theta_0$ (degrees)
-   $K_2$ (energy)
-   $K_3$ (energy)
-   $K_4$ (energy)
-   $K_5$ (energy)
-   $K_6$ (energy)

A pflag value of 0 vs 1 selects between the \"in-plane\" and
\"out-of-plane\" options described above. Ubflag is 1 if there is a
Urey-Bradley term associated with this angle type, else it is 0.
$\theta_0$ is specified in degrees, but LAMMPS converts it to radians
internally; hence the various $K$ values are effectively energy per
radian\^2 or radian\^3 or radian\^4 or radian\^5 or radian\^6.

For the $E_{ba}$ formula, each line in a [angle_coeff](angle_coeff)
command in the input script lists 5 coefficients, the first of which is
\"ba\" to indicate they are BondAngle coefficients. In a data file,
these coefficients should be listed under a \"BondAngle Coeffs\" heading
and you must leave out the \"ba\", i.e. only list 4 coefficients after
the angle type.

-   ba
-   $N_1$ (energy/distance\^2)
-   $N_2$ (energy/distance\^2)
-   $r_1$ (distance)
-   $r_2$ (distance)

The $\theta_0$ value in the $E_{ba}$ formula is not specified, since it
is the same value from the $E_a$ formula.

For the $E_{ub}$ formula, each line in a [angle_coeff](angle_coeff)
command in the input script lists 3 coefficients, the first of which is
\"ub\" to indicate they are UreyBradley coefficients. In a data file,
these coefficients should be listed under a \"UreyBradley Coeffs\"
heading and you must leave out the \"ub\", i.e. only list 2 coefficients
after the angle type.

-   ub
-   $K_{ub}$ (energy/distance\^2)
-   $r_{ub}$ (distance)

------------------------------------------------------------------------

## Restrictions

This angle style can only be used if LAMMPS was built with the AMOEBA
package. See the [Build package](Build_package) doc page for more info.

## Related commands

[angle_coeff](angle_coeff)

## Default

none
