# bond_style class2 command

Accelerator Variants: *class2/omp*, *class2/kk*

## Syntax

``` LAMMPS
bond_style class2
```

## Examples

``` LAMMPS
bond_style class2
bond_coeff 1 1.0 100.0 80.0 80.0
```

## Description

The *class2* bond style uses the potential

$$E = K_2 (r - r_0)^2 + K_3 (r - r_0)^3 + K_4 (r - r_0)^4$$

where $r_0$ is the equilibrium bond distance.

See [(Sun)](bond-Sun) for a description of the COMPASS class2 force
field.

The following coefficients must be defined for each bond type via the
[bond_coeff](bond_coeff) command as in the example above, or in the data
file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   $r_0$ (distance)
-   $K_2$ (energy/distance\^2)
-   $K_3$ (energy/distance\^3)
-   $K_4$ (energy/distance\^4)

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

This bond style can only be used if LAMMPS was built with the CLASS2
package. See the [Build package](Build_package) page for more info.

## Related commands

[bond_coeff](bond_coeff), [delete_bonds](delete_bonds)

## Default

none

------------------------------------------------------------------------

::: {#bond-Sun}
**(Sun)** Sun, J Phys Chem B 102, 7338-7364 (1998).
:::
