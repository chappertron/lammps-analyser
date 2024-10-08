# angle_style class2 command

Accelerator Variants: *class2/kk*, *class2/omp*

# angle_style class2/p6 command

## Syntax

``` LAMMPS
angle_style class2
```

## Examples

``` LAMMPS
angle_style class2
angle_coeff * 75.0 25.0 0.3 0.002
angle_coeff 1 bb 10.5872 1.0119 1.5228
angle_coeff * ba 3.6551 24.895 1.0119 1.5228
```

## Description

The *class2* angle style uses the potential

$$\begin{aligned}
E & = E_a + E_{bb} + E_{ba} \\
E_a & = K_2 (\theta - \theta_0)^2 + K_3 (\theta - \theta_0)^3 + K_4(\theta - \theta_0)^4 \\
E_{bb} & = M (r_{ij} - r_1) (r_{jk} - r_2) \\
E_{ba} & = N_1 (r_{ij} - r_1) (\theta - \theta_0) + N_2(r_{jk} - r_2)(\theta - \theta_0)
\end{aligned}$$

where $E_a$ is the angle term, $E_{bb}$ is a bond-bond term, and
$E_{ba}$ is a bond-angle term. $\theta_0$ is the equilibrium angle and
$r_1$ and $r_2$ are the equilibrium bond lengths.

See [(Sun)](angle-Sun) for a description of the COMPASS class2 force
field.

Coefficients for the $E_a$, $E_{bb}$, and $E_{ba}$ formulas must be
defined for each angle type via the [angle_coeff](angle_coeff) command
as in the example above, or in the data file or restart files read by
the [read_data](read_data) or [read_restart](read_restart) commands.

These are the 4 coefficients for the $E_a$ formula:

-   $\theta_0$ (degrees)
-   $K_2$ (energy)
-   $K_3$ (energy)
-   $K_4$ (energy)

$\theta_0$ is specified in degrees, but LAMMPS converts it to radians
internally; hence the various $K$ are effectively energy per radian\^2
or radian\^3 or radian\^4.

For the $E_{bb}$ formula, each line in a [angle_coeff](angle_coeff)
command in the input script lists 4 coefficients, the first of which is
\"bb\" to indicate they are BondBond coefficients. In a data file, these
coefficients should be listed under a \"BondBond Coeffs\" heading and
you must leave out the \"bb\", i.e. only list 3 coefficients after the
angle type.

-   bb
-   $M$ (energy/distance\^2)
-   $r_1$ (distance)
-   $r_2$ (distance)

For the $E_{ba}$ formula, each line in a [angle_coeff](angle_coeff)
command in the input script lists 5 coefficients, the first of which is
\"ba\" to indicate they are BondAngle coefficients. In a data file,
these coefficients should be listed under a \"BondAngle Coeffs\" heading
and you must leave out the \"ba\", i.e. only list 4 coefficients after
the angle type.

-   ba
-   $N_1$ (energy/distance)
-   $N_2$ (energy/distance)
-   $r_1$ (distance)
-   $r_2$ (distance)

The $\theta_0$ value in the $E_{ba}$ formula is not specified, since it
is the same value from the $E_a$ formula.

:::: note
::: title
Note
:::

It is important that the order of the I,J,K atoms in each angle listed
in the Angles section of the data file read by the
[read_data](read_data) command be consistent with the order of the $r_1$
and $r_2$ BondBond and BondAngle coefficients. This is because the terms
in the formulas for $E_{bb}$ and $E_{ba}$ will use the I,J atoms to
compute $r_{ij}$ and the J,K atoms to compute $r_{jk}$.
::::

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

The *class2/p6* angle style uses the *class2* potential expanded to
sixth order:

$$E_{a} = K_2\left(\theta - \theta_0\right)^2 + K_3\left(\theta - \theta_0\right)^3 + K_4\left(\theta - \theta_0\right)^4 + K_5\left(\theta - \theta_0\right)^5 + K_6\left(\theta - \theta_0\right)^6$$

In this expanded term 6 coefficients for the $E_a$ formula need to be
set:

-   $\theta_0$ (degrees)
-   $K_2$ (energy)
-   $K_3$ (energy)
-   $K_4$ (energy)
-   $K_5$ (energy)
-   $K_6$ (energy)

$\theta_0$ is specified in degrees, but LAMMPS converts it to radians
internally; hence the various $K$ are effectively energy per radian\^2
or radian\^3 or radian\^4 or radian\^5 or radian\^6.

The bond-bond and bond-angle terms remain unchanged.

------------------------------------------------------------------------

## Restrictions

This angle style can only be used if LAMMPS was built with the CLASS2
package. For the *class2/p6* style LAMMPS needs to be built with the
MOFFF package. See the [Build package](Build_package) doc page for more
info.

## Related commands

[angle_coeff](angle_coeff)

## Default

none

------------------------------------------------------------------------

::: {#angle-Sun}
**(Sun)** Sun, J Phys Chem B 102, 7338-7364 (1998).
:::
