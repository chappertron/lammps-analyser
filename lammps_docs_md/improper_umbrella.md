# improper_style umbrella command

Accelerator Variants: *umbrella/omp*

## Syntax

``` LAMMPS
improper_style umbrella
```

## Examples

``` LAMMPS
improper_style umbrella
improper_coeff 1 100.0 180.0
```

## Description

The *umbrella* improper style uses the following potential, which is
commonly referred to as a classic inversion and used in the
[DREIDING](Howto_bioFF) force field:

$$\begin{aligned}
E = & \frac{1}{2}K\left( \frac{1}{\sin\omega_0}\right) ^2 \left( \cos\omega - \cos\omega_0\right) ^2 \qquad \omega_0 \neq 0^o \\
E = & K\left( 1-cos\omega\right)  \qquad \omega_0 = 0^o
\end{aligned}$$

where $K$ is the force constant and $\omega$ is the angle between the IL
axis and the IJK plane:

![image](JPG/umbrella.jpg){.align-center}

If $\omega_0 = 0$ the potential term has a minimum for the planar
structure. Otherwise it has two minima at $\omega +/- \omega_0$, with a
barrier in between.

See [(Mayo)](umbrella-Mayo) for a description of the DREIDING force
field.

The following coefficients must be defined for each improper type via
the [improper_coeff](improper_coeff) command as in the example above, or
in the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   $K$ (energy)
-   $\omega_0$ (degrees)

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

This improper style can only be used if LAMMPS was built with the
MOLECULE package. See the [Build package](Build_package) doc page for
more info.

## Related commands

[improper_coeff](improper_coeff)

## Default

none

------------------------------------------------------------------------

::: {#umbrella-Mayo}
**(Mayo)** Mayo, Olfason, Goddard III, J Phys Chem, 94, 8897-8909
(1990),
:::
