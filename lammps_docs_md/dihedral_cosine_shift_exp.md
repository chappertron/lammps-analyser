# dihedral_style cosine/shift/exp command

Accelerator Variants: *cosine/shift/exp/omp*

## Syntax

``` LAMMPS
dihedral_style cosine/shift/exp
```

## Examples

``` LAMMPS
dihedral_style cosine/shift/exp
dihedral_coeff 1 10.0 45.0 2.0
```

## Description

The *cosine/shift/exp* dihedral style uses the potential

$$E = -U_{min}\frac{e^{-a U(\theta,\theta_0)}-1}{e^a-1} \quad\mbox{with}\quad U(\theta,\theta_0)=-0.5 \left(1+\cos(\theta-\theta_0) \right)$$

where $U_{min}$, $\theta$, and $a$ are defined for each dihedral type.

The potential is bounded between $\left[-U_{min}:0\right]$ and the
minimum is located at the angle $\theta_0$. The a parameter can be both
positive or negative and is used to control the spring constant at the
equilibrium.

The spring constant is given by
$k=a e^a \frac{U_{min}}{2 \left(e^a-1\right)}$. For $a>3$ and
$\frac{k}{U_{min}} = \frac{a}{2}$ to better than 5% relative error. For
negative values of the a parameter, the spring constant is essentially
zero, and anharmonic terms takes over. The potential is furthermore well
behaved in the limit $a \rightarrow 0$, where it has been implemented to
linear order in $a$ for $a < 0.001$.

The following coefficients must be defined for each dihedral type via
the [dihedral_coeff](dihedral_coeff) command as in the example above, or
in the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   $U_{min}$ (energy)
-   $\theta$ (angle)
-   $a$ (real number)

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

This dihedral style can only be used if LAMMPS was built with the
MOLECULE package. See the [Build package](Build_package) doc page for
more info.

## Related commands

[dihedral_coeff](dihedral_coeff), [angle_style
cosine/shift/exp](angle_cosine_shift_exp)

## Default

none
