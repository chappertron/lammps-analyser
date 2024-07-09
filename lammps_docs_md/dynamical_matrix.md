# dynamical_matrix command

Accelerator Variants: dynamical_matrix/kk

## Syntax

``` LAMMPS
dynamical_matrix group-ID style gamma args keyword value ...
```

-   group-ID = ID of group of atoms to displace

-   style = *regular* or *eskm*

-   gamma = finite different displacement length (distance units)

-   one or more keyword/arg pairs may be appended

        keyword = *file* or *binary*
          *file* name = name of output file for the dynamical matrix
          *binary* arg = *yes* or *no* or *gzip*

## Examples

``` LAMMPS
dynamical_matrix 1 regular 0.000001
dynamical_matrix 1 eskm 0.000001
dynamical_matrix 3 regular 0.00004 file dynmat.dat
dynamical_matrix 5 eskm 0.00000001 file dynamical.dat binary yes
```

## Description

Calculate the dynamical matrix by finite difference of the selected
group,

$$D = \frac{\Phi_{ij}^{\alpha\beta}}{\sqrt{M_i M_j}}$$

where D is the dynamical matrix and $\Phi$ is the force constant matrix
defined by

$$\Phi_{ij}^{\alpha\beta} = \frac{\partial^2 U}{\partial x_{i,\alpha} \partial x_{j,\beta}}$$

The output for the dynamical matrix is printed three elements at a time.
The three elements are the three $\beta$ elements for a respective
i/$\alpha$/j combination. Each line is printed in order of j increasing
first, $\alpha$ second, and i last.

If the style eskm is selected, the dynamical matrix will be in units of
inverse squared femtoseconds. These units will then conveniently leave
frequencies in THz.

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

The command collects an array of nine times the number of atoms in a
group on every single MPI rank, so the memory requirements can be very
significant for large systems.

This command is part of the PHONON package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

## Related commands

[fix phonon](fix_phonon), [fix numdiff](fix_numdiff),

[compute hma](compute_hma) uses an analytic formulation of the Hessian
provided by a pair_style\'s Pair::single_hessian() function, if
implemented.

## Default

The default settings are file = \"dynmat.dyn\", binary = no
