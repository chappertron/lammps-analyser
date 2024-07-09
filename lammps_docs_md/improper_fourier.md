# improper_style fourier command

Accelerator Variants: *fourier/omp*

## Syntax

``` LAMMPS
improper_style fourier
```

## Examples

``` LAMMPS
improper_style fourier
improper_coeff 1 100.0 0.0 1.0 0.5 1
```

## Description

The *fourier* improper style uses the following potential:

$$E = K [C_0 + C_1 \cos ( \omega) + C_2 \cos( 2 \omega) ]$$

where K is the force constant, C0, C1, C2 are dimensionless
coefficients, and omega is the angle between the IL axis and the IJK
plane:

![image](JPG/umbrella.jpg){.align-center}

If all parameter (see below) is not zero, the all the three possible
angles will taken in account.

The following coefficients must be defined for each improper type via
the [improper_coeff](improper_coeff) command as in the example above, or
in the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   $K$ (energy)
-   $C_0$ (unitless)
-   $C_1$ (unitless)
-   $C_2$ (unitless)
-   all (0 or 1, optional)

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

This angle style can only be used if LAMMPS was built with the MOLECULE
package. See the [Build package](Build_package) doc page for more info.

## Related commands

[improper_coeff](improper_coeff)

## Default

none
