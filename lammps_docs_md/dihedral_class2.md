# dihedral_style class2 command

Accelerator Variants: *class2/omp*, *class2/kk*

## Syntax

``` LAMMPS
dihedral_style class2
```

## Examples

``` LAMMPS
dihedral_style class2
dihedral_coeff 1 100 75 100 70 80 60
dihedral_coeff * mbt 3.5945 0.1704 -0.5490 1.5228
dihedral_coeff * ebt 0.3417 0.3264 -0.9036 0.1368 0.0 -0.8080 1.0119 1.1010
dihedral_coeff 2 at 0.0 -0.1850 -0.7963 -2.0220 0.0 -0.3991 110.2453 105.1270
dihedral_coeff * aat -13.5271 110.2453 105.1270
dihedral_coeff * bb13 0.0 1.0119 1.1010
```

## Description

The *class2* dihedral style uses the potential

$$\begin{aligned}
E        = & E_d + E_{mbt} + E_{ebt} + E_{at} + E_{aat} + E_{bb13} \\
E_d      = & \sum_{n=1}^{3} K_n [ 1 - \cos (n \phi - \phi_n) ] \\
E_{mbt}  = & (r_{jk} - r_2) [ A_1 \cos (\phi) + A_2 \cos (2\phi) + A_3 \cos (3\phi) ] \\
E_{ebt}  = & (r_{ij} - r_1) [ B_1 \cos (\phi) + B_2 \cos (2\phi) + B_3 \cos (3\phi) ] + \\
           & (r_{kl} - r_3) [ C_1 \cos (\phi) + C_2 \cos (2\phi) + C_3 \cos (3\phi) ] \\
E_{at}   = & (\theta_{ijk} - \theta_1) [ D_1 \cos (\phi) + D_2 \cos (2\phi) + D_3 \cos (3\phi) ] + \\
           & (\theta_{jkl} - \theta_2) [ E_1 \cos (\phi) + E_2 \cos (2\phi) + E_3 \cos (3\phi) ] \\
E_{aat}  = & M (\theta_{ijk} - \theta_1) (\theta_{jkl} - \theta_2) \cos (\phi) \\
E_{bb13} = & N (r_{ij} - r_1) (r_{kl} - r_3)
\end{aligned}$$

where $E_d$ is the dihedral term, $E_{mbt}$ is a middle-bond-torsion
term, $E_{ebt}$ is an end-bond-torsion term, $E_{at}$ is an
angle-torsion term, $E_{aat}$ is an angle-angle-torsion term, and
$E_{bb13}$ is a bond-bond-13 term.

$\theta_1$ and $\theta_2$ are equilibrium angles and $r_1$, $r_2$, and
$r_3$ are equilibrium bond lengths.

See [(Sun)](dihedral-Sun) for a description of the COMPASS class2 force
field.

Coefficients for the $E_d$, $E_{mbt}$, $E_{ebt}$, $E_{at}$, $E_{aat}$,
and $E_{bb13}$ formulas must be defined for each dihedral type via the
[dihedral_coeff](dihedral_coeff) command as in the example above, or in
the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands.

These are the 6 coefficients for the $E_d$ formula:

-   $K_1$ (energy)
-   $\phi_1$ (degrees)
-   $K_2$ (energy)
-   $\phi_2$ (degrees)
-   $K_3$ (energy)
-   $phi_3$ (degrees)

For the $E_{mbt}$ formula, each line in a
[dihedral_coeff](dihedral_coeff) command in the input script lists 5
coefficients, the first of which is *mbt* to indicate they are
MiddleBondTorsion coefficients. In a data file, these coefficients
should be listed under a *MiddleBondTorsion Coeffs* heading and you must
leave out the *mbt*, i.e. only list 4 coefficients after the dihedral
type.

-   *mbt*
-   $A_1$ (energy/distance)
-   $A_2$ (energy/distance)
-   $A_3$ (energy/distance)
-   $r_2$ (distance)

For the $E_{ebt}$ formula, each line in a
[dihedral_coeff](dihedral_coeff) command in the input script lists 9
coefficients, the first of which is *ebt* to indicate they are
EndBondTorsion coefficients. In a data file, these coefficients should
be listed under a *EndBondTorsion Coeffs* heading and you must leave out
the *ebt*, i.e. only list 8 coefficients after the dihedral type.

-   *ebt*
-   $B_1$ (energy/distance)
-   $B_2$ (energy/distance)
-   $B_3$ (energy/distance)
-   $C_1$ (energy/distance)
-   $C_2$ (energy/distance)
-   $C_3$ (energy/distance)
-   $r_1$ (distance)
-   $r_3$ (distance)

For the $E_{at}$ formula, each line in a
[dihedral_coeff](dihedral_coeff) command in the input script lists 9
coefficients, the first of which is *at* to indicate they are
AngleTorsion coefficients. In a data file, these coefficients should be
listed under a *AngleTorsion Coeffs* heading and you must leave out the
*at*, i.e. only list 8 coefficients after the dihedral type.

-   *at*
-   $D_1$ (energy)
-   $D_2$ (energy)
-   $D_3$ (energy)
-   $E_1$ (energy)
-   $E_2$ (energy)
-   $E_3$ (energy)
-   $\theta_1$ (degrees)
-   $\theta_2$ (degrees)

$\theta_1$ and $\theta_2$ are specified in degrees, but LAMMPS converts
them to radians internally; hence the various $D$ and $E$ are
effectively energy per radian.

For the $E_{aat}$ formula, each line in a
[dihedral_coeff](dihedral_coeff) command in the input script lists 4
coefficients, the first of which is *aat* to indicate they are
AngleAngleTorsion coefficients. In a data file, these coefficients
should be listed under a *AngleAngleTorsion Coeffs* heading and you must
leave out the *aat*, i.e. only list 3 coefficients after the dihedral
type.

-   *aat*
-   $M$ (energy)
-   $\theta_1$ (degrees)
-   $\theta_2$ (degrees)

$\theta_1$ and $\theta_2$ are specified in degrees, but LAMMPS converts
them to radians internally; hence $M$ is effectively energy per
radian\^2.

For the $E_{bb13}$ formula, each line in a
[dihedral_coeff](dihedral_coeff) command in the input script lists 4
coefficients, the first of which is *bb13* to indicate they are
BondBond13 coefficients. In a data file, these coefficients should be
listed under a *BondBond13 Coeffs* heading and you must leave out the
*bb13*, i.e. only list 3 coefficients after the dihedral type.

-   *bb13*
-   $N$ (energy/distance\^2)
-   $r_1$ (distance)
-   $r_3$ (distance)

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

------------------------------------------------------------------------

## Restrictions

This dihedral style can only be used if LAMMPS was built with the CLASS2
package. See the [Build package](Build_package) doc page for more info.

## Related commands

[dihedral_coeff](dihedral_coeff)

## Default

none

------------------------------------------------------------------------

::: {#dihedral-Sun}
**(Sun)** Sun, J Phys Chem B 102, 7338-7364 (1998).
:::
