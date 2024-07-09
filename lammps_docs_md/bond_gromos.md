# bond_style gromos command

Accelerator Variants: *gromos/omp*

## Syntax

``` LAMMPS
bond_style gromos
```

## Examples

``` LAMMPS
bond_style gromos
bond_coeff 5 80.0 1.2
```

## Description

The *gromos* bond style uses the potential

$$E = K (r^2 - r_0^2)^2$$

where $r_0$ is the equilibrium bond distance. Note that the usual 1/4
factor is included in $K$.

The following coefficients must be defined for each bond type via the
[bond_coeff](bond_coeff) command as in the example above, or in the data
file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   $K$ (energy/distance\^4)
-   $r_0$ (distance)

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

This bond style can only be used if LAMMPS was built with the MOLECULE
package. See the [Build package](Build_package) page for more info.

## Related commands

[bond_coeff](bond_coeff), [delete_bonds](delete_bonds)

## Default

none
