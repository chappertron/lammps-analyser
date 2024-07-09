# angle_style cosine/squared command

Accelerator Variants: *cosine/squared/omp*

## Syntax

``` LAMMPS
angle_style cosine/squared
```

## Examples

``` LAMMPS
angle_style cosine/squared
angle_coeff 2*4 75.0 100.0
```

## Description

The *cosine/squared* angle style uses the potential

$$E = K [\cos(\theta) - \cos(\theta_0)]^2$$

, which is commonly used in the [DREIDING](Howto_bioFF) force field,
where $\theta_0$ is the equilibrium value of the angle, and $K$ is a
prefactor. Note that the usual 1/2 factor is included in $K$.

See [(Mayo)](cosine-Mayo) for a description of the DREIDING force field.

The following coefficients must be defined for each angle type via the
[angle_coeff](angle_coeff) command as in the example above, or in the
data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   $K$ (energy)
-   $\theta_0$ (degrees)

$\theta_0$ is specified in degrees, but LAMMPS converts it to radians
internally.

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

[angle_coeff](angle_coeff)

## Default

none

------------------------------------------------------------------------

::: {#cosine-Mayo}
**(Mayo)** Mayo, Olfason, Goddard III, J Phys Chem, 94, 8897-8909
(1990).
:::
