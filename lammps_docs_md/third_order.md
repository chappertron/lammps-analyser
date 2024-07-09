# third_order command

Accelerator Variants: third_order/kk

## Syntax

``` LAMMPS
third_order group-ID style delta args keyword value ...
```

-   group-ID = ID of group of atoms to displace

-   style = *regular* or *eskm*

-   delta = finite different displacement length (distance units)

-   one or more keyword/arg pairs may be appended

        keyword = *file* or *binary*
          *file* name = name of output file for the third order tensor
          *binary* arg = *yes* or *no* or *gzip*

## Examples

``` LAMMPS
third_order 1 regular 0.000001
third_order 1 eskm 0.000001
third_order 3 regular 0.00004 file third_order.dat
third_order 5 eskm 0.00000001 file third_order.dat binary yes
```

## Description

Calculate the third order force constant tensor by finite difference of
the selected group,

$$\Phi^{\alpha\beta\gamma}_{ijk} = \frac{\partial^3 U}{\partial x_{i,\alpha} \partial x_{j,\beta} \partial x_{k, \gamma}}$$

where Phi is the third order force constant tensor.

The output of the command is the tensor, three elements at a time. The
three elements correspond to the three gamma elements for a specific
i/alpha/j/beta/k. The initial five numbers are i, alpha, j, beta, and k
respectively.

If the style eskm is selected, the tensor will be using energy units of
10 J/mol. These units conform to eskm style from the dynamical_matrix
command, which will simplify operations using dynamical matrices with
third order tensors.

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

The command collects a 9 times the number of atoms in the group on every
single MPI rank, so the memory requirements can be very significant for
large systems.

This command is part of the PHONON package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

## Related commands

[fix phonon](fix_phonon) [dynamical_matrix](dynamical_matrix)

## Default

The default settings are file = \"third_order.dat\", binary = no
