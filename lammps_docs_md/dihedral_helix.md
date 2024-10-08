# dihedral_style helix command

Accelerator Variants: *helix/omp*

## Syntax

``` LAMMPS
dihedral_style helix
```

## Examples

``` LAMMPS
dihedral_style helix
dihedral_coeff 1 80.0 100.0 40.0
```

## Description

The *helix* dihedral style uses the potential

$$E = A [1 - \cos(\theta)] + B [1 + \cos(3 \theta)] +
    C [1 + \cos(\theta + \frac{\pi}{4})]$$

This coarse-grain dihedral potential is described in [(Guo)](Guo). For
dihedral angles in the helical region, the energy function is
represented by a standard potential consisting of three minima, one
corresponding to the trans (t) state and the other to gauche states (g+
and g-). The paper describes how the $A$, $B$ and, $C$ parameters are
chosen so as to balance secondary (largely driven by local interactions)
and tertiary structure (driven by long-range interactions).

The following coefficients must be defined for each dihedral type via
the [dihedral_coeff](dihedral_coeff) command as in the example above, or
in the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   $A$ (energy)
-   $B$ (energy)
-   $C$ (energy)

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
EXTRA-MOLECULE package. See the [Build package](Build_package) doc page
for more info.

## Related commands

[dihedral_coeff](dihedral_coeff)

## Default

none

------------------------------------------------------------------------

::: {#Guo}
**(Guo)** Guo and Thirumalai, Journal of Molecular Biology, 263, 323-43
(1996).
:::
