# pair_style snap command

Accelerator Variants: *snap/kk*

## Syntax

``` LAMMPS
pair_style snap
```

## Examples

``` LAMMPS
pair_style snap
pair_coeff * * InP.snapcoeff InP.snapparam In In P P
```

## Description

Pair style *snap* defines the spectral neighbor analysis potential
(SNAP), a machine-learning interatomic potential
[(Thompson)](Thompson20142). Like the GAP framework of Bartok et al.
[(Bartok2010)](Bartok20102), SNAP uses bispectrum components to
characterize the local neighborhood of each atom in a very general way.
The mathematical definition of the bispectrum calculation and its
derivatives w.r.t. atom positions is identical to that used by [compute
snap](compute_sna_atom), which is used to fit SNAP potentials to *ab
initio* energy, force, and stress data. In SNAP, the total energy is
decomposed into a sum over atom energies. The energy of atom *i* is
expressed as a weighted sum over bispectrum components.

$$E^i_{SNAP}(B_1^i,...,B_K^i) = \beta^{\mu_i}_0 + \sum_{k=1}^K \beta_k^{\mu_i} B_k^i$$

where $B_k^i$ is the *k*-th bispectrum component of atom *i*, and
$\beta_k^{\mu_i}$ is the corresponding linear coefficient that depends
on $\mu_i$, the SNAP element of atom *i*. The number of bispectrum
components used and their definitions depend on the value of *twojmax*
and other parameters defined in the SNAP parameter file described below.
The bispectrum calculation is described in more detail in [compute
sna/atom](compute_sna_atom).

Note that unlike for other potentials, cutoffs for SNAP potentials are
not set in the pair_style or pair_coeff command; they are specified in
the SNAP potential files themselves.

Only a single pair_coeff command is used with the *snap* style which
specifies a SNAP coefficient file followed by a SNAP parameter file and
then N additional arguments specifying the mapping of SNAP elements to
LAMMPS atom types, where N is the number of LAMMPS atom types:

-   SNAP coefficient file
-   SNAP parameter file
-   N element names = mapping of SNAP elements to atom types

As an example, if a LAMMPS indium phosphide simulation has 4 atoms
types, with the first two being indium and the third and fourth being
phophorous, the pair_coeff command would look like this:

``` LAMMPS
pair_coeff * * snap InP.snapcoeff InP.snapparam In In P P
```

The first 2 arguments must be \* \* so as to span all LAMMPS atom types.
The two filenames are for the coefficient and parameter files,
respectively. The two trailing \'In\' arguments map LAMMPS atom types 1
and 2 to the SNAP \'In\' element. The two trailing \'P\' arguments map
LAMMPS atom types 3 and 4 to the SNAP \'P\' element.

If a SNAP mapping value is specified as NULL, the mapping is not
performed. This can be used when a *snap* potential is used as part of
the *hybrid* pair style. The NULL values are placeholders for atom types
that will be used with other potentials.

The name of the SNAP coefficient file usually ends in the \".snapcoeff\"
extension. It may contain coefficients for many SNAP elements. The only
requirement is that each of the unique element names appearing in the
LAMMPS pair_coeff command appear exactly once in the SNAP coefficient
file. It is okay if the SNAP coefficient file contains additional
elements not in the pair_coeff command, except when using *chemflag*
(see below). The name of the SNAP parameter file usually ends in the
\".snapparam\" extension. It contains a small number of parameters that
define the overall form of the SNAP potential. See the
[pair_coeff](pair_coeff) page for alternate ways to specify the path for
these files.

SNAP potentials are quite commonly combined with one or more other
LAMMPS pair styles using the *hybrid/overlay* pair style. As an example,
the SNAP tantalum potential provided in the LAMMPS potentials directory
combines the *snap* and *zbl* pair styles. It is invoked by the
following commands:

``` LAMMPS
variable zblcutinner equal 4
variable zblcutouter equal 4.8
variable zblz equal 73
pair_style hybrid/overlay &
zbl ${zblcutinner} ${zblcutouter} snap
pair_coeff * * zbl 0.0
pair_coeff 1 1 zbl ${zblz}
pair_coeff * * snap Ta06A.snapcoeff Ta06A.snapparam Ta
```

It is convenient to keep these commands in a separate file that can be
inserted in any LAMMPS input script using the [include](include)
command.

The top of the SNAP coefficient file can contain any number of blank and
comment lines (start with #), but follows a strict format after that.
The first non-blank non-comment line must contain two integers:

-   nelem = Number of elements
-   ncoeff = Number of coefficients

This is followed by one block for each of the *nelem* elements. The
first line of each block contains three entries:

-   Element name (text string)
-   R = Element radius (distance units)
-   w = Element weight (dimensionless)

This line is followed by *ncoeff* coefficients, one per line.

The SNAP parameter file can contain blank and comment lines (start with
#) anywhere. Each non-blank non-comment line must contain one
keyword/value pair. The required keywords are *rcutfac* and *twojmax*.
Optional keywords are *rfac0*, *rmin0*, *switchflag*, *bzeroflag*,
*quadraticflag*, *chemflag*, *bnormflag*, *wselfallflag*,
*switchinnerflag*, *sinner*, *dinner*, *chunksize*, and
*parallelthresh*.

The default values for these keywords are

-   *rfac0* = 0.99363
-   *rmin0* = 0.0
-   *switchflag* = 1
-   *bzeroflag* = 1
-   *quadraticflag* = 0
-   *chemflag* = 0
-   *bnormflag* = 0
-   *wselfallflag* = 0
-   *switchinnerflag* = 0
-   *chunksize* = 32768
-   *parallelthresh* = 8192

For detailed definitions of all of these keywords, see the [compute
sna/atom](compute_sna_atom) doc page.

If *quadraticflag* is set to 1, then the SNAP energy expression includes
additional quadratic terms that have been shown to increase the overall
accuracy of the potential without much increase in computational cost
[(Wood)](Wood20182).

$$E^i_{SNAP}(\mathbf{B}^i) = \beta^{\mu_i}_0 + \boldsymbol{\beta}^{\mu_i} \cdot \mathbf{B}_i + \frac{1}{2}\mathbf{B}^t_i \cdot \boldsymbol{\alpha}^{\mu_i} \cdot \mathbf{B}_i$$

where $\mathbf{B}_i$ is the *K*-vector of bispectrum components,
$\boldsymbol{\beta}^{\mu_i}$ is the *K*-vector of linear coefficients
for element $\mu_i$, and $\boldsymbol{\alpha}^{\mu_i}$ is the symmetric
*K* by *K* matrix of quadratic coefficients. The SNAP coefficient file
should contain *K*(\*K\*+1)/2 additional coefficients in each element
block, the upper-triangular elements of $\boldsymbol{\alpha}^{\mu_i}$.

If *chemflag* is set to 1, then the energy expression is written in
terms of explicit multi-element bispectrum components indexed on ordered
triplets of elements, which has been shown to increase the ability of
the SNAP potential to capture energy differences in chemically complex
systems, at the expense of a significant increase in computational cost
[(Cusentino)](Cusentino20202).

$$E^i_{SNAP}(\mathbf{B}^i) = \beta^{\mu_i}_0 + \sum_{\kappa,\lambda,\mu} \boldsymbol{\beta}^{\kappa\lambda\mu}_{\mu_i} \cdot \mathbf{B}^{\kappa\lambda\mu}_i$$

where $\mathbf{B}^{\kappa\lambda\mu}_i$ is the *K*-vector of bispectrum
components for neighbors of elements $\kappa$, $\lambda$, and $\mu$ and
$\boldsymbol{\beta}^{\kappa\lambda\mu}_{\mu_i}$ is the corresponding
*K*-vector of linear coefficients for element $\mu_i$. The SNAP
coefficient file should contain a total of $K N_{elem}^3$ coefficients
in each element block, where $N_{elem}$ is the number of elements in the
SNAP coefficient file, which must equal the number of unique elements
appearing in the LAMMPS pair_coeff command, to avoid ambiguity in the
number of coefficients.

The keyword *switchinnerflag* activates an additional switching function
that smoothly turns off contributions to the SNAP potential from
neighbor atoms at short separations. If *switchinnerflag* is set to 1
then the additional keywords *sinner* and *dinner* must also be
provided. Each of these is followed by *nelements* values, where
*nelements* is the number of unique elements appearing in appearing in
the LAMMPS pair_coeff command. The element order should correspond to
the order in which elements first appear in the pair_coeff command
reading from left to right.

The keywords *chunksize* and *parallelthresh* are only applicable when
using the pair style *snap* with the KOKKOS package on GPUs and are
ignored otherwise. The *chunksize* keyword controls the number of atoms
in each pass used to compute the bispectrum components and is used to
avoid running out of memory. For example if there are 8192 atoms in the
simulation and the *chunksize* is set to 4096, the bispectrum
calculation will be broken up into two passes (running on a single GPU).
The *parallelthresh* keyword controls a crossover threshold for
performing extra parallelism. For small systems, exposing additional
parallelism can be beneficial when there is not enough work to fully
saturate the GPU threads otherwise. However, the extra parallelism also
leads to more divergence and can hurt performance when the system is
already large enough to saturate the GPU threads. Extra parallelism will
be performed if the *chunksize* (or total number of atoms per GPU) is
smaller than *parallelthresh*.

:::: note
::: title
Note
:::

The previously used *diagonalstyle* keyword was removed in 2019, since
all known SNAP potentials use the default value of 3.
::::

------------------------------------------------------------------------

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

This style is part of the ML-SNAP package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

## Related commands

[compute sna/atom](compute_sna_atom), [compute
snad/atom](compute_sna_atom), [compute snav/atom](compute_sna_atom),
[compute snap](compute_sna_atom)

## Default

none

------------------------------------------------------------------------

::: {#Thompson20142}
**(Thompson)** Thompson, Swiler, Trott, Foiles, Tucker, J Comp Phys,
285, 316 (2015).
:::

::: {#Bartok20102}
**(Bartok2010)** Bartok, Payne, Kondor, Csanyi, Phys Rev Lett, 104,
136403 (2010).
:::

::: {#Wood20182}
**(Wood)** Wood and Thompson, J Chem Phys, 148, 241721, (2018)
:::

::: {#Cusentino20202}
**(Cusentino)** Cusentino, Wood, and Thompson, J Phys Chem A, xxx,
xxxxx, (2020)
:::
