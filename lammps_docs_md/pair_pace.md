# pair_style pace command

Accelerator Variants: *pace/kk*, *pace/extrapolation/kk*

# pair_style pace/extrapolation command

## Syntax

``` LAMMPS
pair_style pace ... keyword values ...
```

-   one or more keyword/value pairs may be appended

        keyword = *product* or *recursive* or *chunksize*
          *product* = use product algorithm for basis functions
          *recursive* = use recursive algorithm for basis functions
          *chunksize* value = number of atoms in each pass

``` LAMMPS
pair_style pace/extrapolation
```

## Examples

``` LAMMPS
pair_style pace
pair_style pace product chunksize 2048
pair_coeff * * Cu-PBE-core-rep.ace Cu

pair_style pace/extrapolation
pair_coeff * * Cu.yaml Cu.asi Cu
```

## Description

Pair style *pace* computes interactions using the Atomic Cluster
Expansion (ACE), which is a general expansion of the atomic energy in
multi-body basis functions. [(Drautz19)](Drautz20191). The *pace* pair
style provides an efficient implementation that is described in this
paper [(Lysogorskiy21)](Lysogorskiy20211).

In ACE, the total energy is decomposed into a sum over atomic energies.
The energy of atom *i* is expressed as a linear or non-linear function
of one or more density functions. By projecting the density onto a local
atomic base, the lowest order contributions to the energy can be
expressed as a set of scalar polynomials in basis function contributions
summed over neighbor atoms.

Only a single pair_coeff command is used with the *pace* style which
specifies an ACE coefficient file followed by N additional arguments
specifying the mapping of ACE elements to LAMMPS atom types, where N is
the number of LAMMPS atom types:

-   ACE coefficient file
-   N element names = mapping of ACE elements to atom types

Only a single pair_coeff command is used with the *pace* style which
specifies an ACE file that fully defines the potential. Note that unlike
for other potentials, cutoffs are not set in the pair_style or
pair_coeff command; they are specified in the ACE file.

The pair_style *pace* command may be followed by the optional keyword
*product* or *recursive*, which determines which of two algorithms is
used for the calculation of basis functions and derivatives. The default
is *recursive*.

The keyword *chunksize* is only applicable when using the pair style
*pace* with the KOKKOS package on GPUs and is ignored otherwise. This
keyword controls the number of atoms in each pass used to compute the
atomic cluster expansion and is used to avoid running out of memory. For
example if there are 8192 atoms in the simulation and the *chunksize* is
set to 4096, the ACE calculation will be broken up into two passes
(running on a single GPU).

## Extrapolation grade

Calculation of extrapolation grade in PACE is implemented in [pair_style
pace/extrapolation]{.title-ref}. It is based on the MaxVol algorithm
similar to Moment Tensor Potential (MTP) by Shapeev et al. and is
described in [(Lysogorskiy23)](Lysogorskiy2023). In order to compute
extrapolation grade one needs to provide:

1.  ACE potential in B-basis form ([.yaml]{.title-ref} format) and
2.  Active Set Inverted (ASI) file for corresponding potential
    ([.asi]{.title-ref} format)

Calculation of extrapolation grades requires matrix-vector
multiplication for each atom and is slower than the usual [pair_style
pace recursive]{.title-ref}, therefore it is *not* computed by default.
Extrapolation grade calculation is involved by [fix pair]{.title-ref},
which requests to compute [gamma]{.title-ref}, as shown in example
below:

``` LAMMPS
pair_style  pace/extrapolation
pair_coeff  * * Cu.yaml Cu.asi Cu

fix pace_gamma all pair 10 pace/extrapolation gamma 1

compute max_pace_gamma all reduce max f_pace_gamma
variable dump_skip equal "c_max_pace_gamma < 5"

dump pace_dump all custom 20 extrapolative_structures.dump id type x y z f_pace_gamma
dump_modify pace_dump skip v_dump_skip

variable max_pace_gamma equal c_max_pace_gamma
fix extreme_extrapolation all halt 10 v_max_pace_gamma > 25
```

Here extrapolation grade gamma is computed every 10 steps and is stored
in [f_pace_gamma]{.title-ref} per-atom variable. The largest value of
extrapolation grade among all atoms in a structure is reduced to
[c_max_pace_gamma]{.title-ref} variable. Only if this value exceeds
extrapolation threshold 5, then the structure will be dumped into
[extrapolative_structures.dump]{.title-ref} file, but not more often
than every 20 steps.

On all other steps [pair_style pace recursive]{.title-ref} will be used.

When using the pair style *pace/extrapolation* with the KOKKOS package
on GPUs product B-basis evaluator is always used and only *linear* ASI
is supported.

------------------------------------------------------------------------

See the [pair_coeff](pair_coeff) page for alternate ways to specify the
path for the ACE coefficient file.

## Mixing, shift, table, tail correction, restart, rRESPA info

For atom type pairs I,J and I != J, where types I and J correspond to
two different element types, mixing is performed by LAMMPS with
user-specifiable parameters as described above. You never need to
specify a pair_coeff command with I != J arguments for this style.

This pair style does not support the [pair_modify](pair_modify) shift,
table, and tail options.

This pair style does not write its information to [binary restart
files](restart), since it is stored in potential files. Thus, you need
to re-specify the pair_style and pair_coeff commands in an input script
that reads a restart file.

This pair style can only be used via the *pair* keyword of the
[run_style respa](run_style) command. It does not support the *inner*,
*middle*, *outer* keywords.

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

This pair style is part of the ML-PACE package. It is only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

## Related commands

[pair_style snap](pair_snap), [fix pair](fix_pair)

## Default

recursive, chunksize = 4096,

::: {#Drautz20191}
**(Drautz19)** Drautz, Phys Rev B, 99, 014104 (2019).
:::

::: {#Lysogorskiy20211}
**(Lysogorskiy21)** Lysogorskiy, van der Oord, Bochkarev, Menon,
Rinaldi, Hammerschmidt, Mrovec, Thompson, Csanyi, Ortner, Drautz, npj
Comp Mat, 7, 97 (2021).
:::

::: {#Lysogorskiy2023}
**(Lysogorskiy23)** Lysogorskiy, Bochkarev, Mrovec, Drautz, Phys Rev
Mater, 7, 043801 (2023) / arXiv:2212.08716 (2022).
:::
