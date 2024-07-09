# dihedral_style nharmonic command

Accelerator Variants: *nharmonic/omp*

## Syntax

``` LAMMPS
dihedral_style nharmonic
```

## Examples

``` LAMMPS
dihedral_style nharmonic
dihedral_coeff * 3 10.0 20.0 30.0
```

## Description

The *nharmonic* dihedral style uses the potential:

$$E = \sum_{i=1,n} A_i  \cos^{i-1}(\phi)$$

The following coefficients must be defined for each dihedral type via
the [dihedral_coeff](dihedral_coeff) command as in the example above, or
in the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   $n$ (integer \>=1)
-   $A_1$ (energy)
-   $A_2$ (energy)
-   \...
-   $A_n$ (energy)

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
