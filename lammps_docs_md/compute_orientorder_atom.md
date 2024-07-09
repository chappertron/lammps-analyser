# compute orientorder/atom command

Accelerator Variants: *orientorder/atom/kk*

## Syntax

``` LAMMPS
compute ID group-ID orientorder/atom keyword values ...
```

-   ID, group-ID are documented in [compute](compute) command

-   orientorder/atom = style name of this compute command

-   one or more keyword/value pairs may be appended

        keyword = *cutoff* or *nnn* or *degrees* or *wl* or *wl/hat* or *components* or *chunksize*
          *cutoff* value = distance cutoff
          *nnn* value = number of nearest neighbors
          *degrees* values = nlvalues, l1, l2,...
          *wl* value = *yes* or *no*
          *wl/hat* value = *yes* or *no*
          *components* value = ldegree
          *chunksize* value = number of atoms in each pass

## Examples

``` LAMMPS
compute 1 all orientorder/atom
compute 1 all orientorder/atom degrees 5 4 6 8 10 12 nnn NULL cutoff 1.5
compute 1 all orientorder/atom wl/hat yes
compute 1 all orientorder/atom components 6
```

## Description

Define a computation that calculates a set of bond-orientational order
parameters $Q_\ell$ for each atom in a group. These order parameters
were introduced by [Steinhardt et al.](Steinhardt) as a way to
characterize the local orientational order in atomic structures. For
each atom, $Q_\ell$ is a real number defined as follows:

$$\begin{aligned}
\bar{Y}_{\ell m} = & \frac{1}{nnn}\sum_{j = 1}^{nnn} Y_{\ell m}\bigl( \theta( {\bf r}_{ij} ), \phi( {\bf r}_{ij} ) \bigr) \\
Q_\ell  = & \sqrt{\frac{4 \pi}{2 \ell  + 1} \sum_{m = -\ell }^{m = \ell } \bar{Y}_{\ell m} \bar{Y}^*_{\ell m}}
\end{aligned}$$

The first equation defines the local order parameters as averages of the
spherical harmonics $Y_{\ell m}$ for each neighbor. These are complex
number components of the 3D analog of the 2D order parameter $q_n$,
which is implemented as LAMMPS compute
[hexorder/atom](compute_hexorder_atom). The summation is over the *nnn*
nearest neighbors of the central atom. The angles $\theta$ and $\phi$
are the standard spherical polar angles defining the direction of the
bond vector $r_{ij}$. The phase and sign of $Y_{\ell m}$ follow the
standard conventions, so that
$\mathrm{sign}(Y_{\ell\ell}(0,0)) = (-1)^\ell$. The second equation
defines $Q_\ell$, which is a rotationally invariant non-negative
amplitude obtained by summing over all the components of degree $\ell$.

The optional keyword *cutoff* defines the distance cutoff used when
searching for neighbors. The default value, also the maximum allowable
value, is the cutoff specified by the pair style.

The optional keyword *nnn* defines the number of nearest neighbors used
to calculate $Q_\ell$. The default value is 12. If the value is NULL,
then all neighbors up to the specified distance cutoff are used.

The optional keyword *degrees* defines the list of order parameters to
be computed. The first argument *nlvalues* is the number of order
parameters. This is followed by that number of non-negative integers
giving the degree of each order parameter. Because $Q_2$ and all
odd-degree order parameters are zero for atoms in cubic crystals (see
[Steinhardt](Steinhardt)), the default order parameters are $Q_4$,
$Q_6$, $Q_8$, $Q_{10}$, and $Q_{12}$. For the FCC crystal with *nnn*
=12,

$$Q_4 = \sqrt{\frac{7}{192}} \approx 0.19094$$

The numerical values of all order parameters up to $Q_{12}$ for a range
of commonly encountered high-symmetry structures are given in Table I of
[Mickel et al.](Mickel), and these can be reproduced with this compute.

The optional keyword *wl* will output the third-order invariants
$W_\ell$ (see Eq. 1.4 in [Steinhardt](Steinhardt)) for the same degrees
as for the $Q_\ell$ parameters. For the FCC crystal with *nnn* = 12,

$$W_4 = -\sqrt{\frac{14}{143}} \left(\frac{49}{4096}\right) \pi^{-3/2} \approx -0.0006722136$$

The optional keyword *wl/hat* will output the normalized third-order
invariants $\hat{W}_\ell$ (see Eq. 2.2 in [Steinhardt](Steinhardt)) for
the same degrees as for the $Q_\ell$ parameters. For the FCC crystal
with *nnn* =12,

$$\hat{W}_4 = -\frac{7}{3} \sqrt{\frac{2}{429}} \approx -0.159317$$

The numerical values of $\hat{W}_\ell$ for a range of commonly
encountered high-symmetry structures are given in Table I of
[Steinhardt](Steinhardt), and these can be reproduced with this keyword.

The optional keyword *components* will output the components of the
*normalized* complex vector
$\hat{Y}_{\ell m} = \bar{Y}_{\ell m}/|\bar{Y}_{\ell m}|$ of degree
*ldegree*, which must be included in the list of order parameters to be
computed. This option can be used in conjunction with [compute
coord_atom](compute_coord_atom) to calculate the ten Wolde\'s criterion
to identify crystal-like particles, as discussed in [ten
Wolde](tenWolde2).

The optional keyword *chunksize* is only applicable when using the the
KOKKOS package and is ignored otherwise. This keyword controls the
number of atoms in each pass used to compute the bond-orientational
order parameters and is used to avoid running out of memory. For example
if there are 32768 atoms in the simulation and the *chunksize* is set to
16384, the parameter calculation will be broken up into two passes.

The value of $Q_\ell$ is set to zero for atoms not in the specified
compute group, as well as for atoms that have less than *nnn* neighbors
within the distance cutoff, unless *nnn* is NULL.

The neighbor list needed to compute this quantity is constructed each
time the calculation is performed (i.e., each time a snapshot of atoms
is dumped). Thus it can be inefficient to compute/dump this quantity too
frequently.

:::: note
::: title
Note
:::

If you have a bonded system, then the settings of
[special_bonds](special_bonds) command can remove pairwise interactions
between atoms in the same bond, angle, or dihedral. This is the default
setting for the [special_bonds](special_bonds) command, and means those
pairwise interactions do not appear in the neighbor list. Because this
fix uses the neighbor list, it also means those pairs will not be
included in the order parameter. This difficulty can be circumvented by
writing a dump file, and using the [rerun](rerun) command to compute the
order parameter for snapshots in the dump file. The rerun script can use
a [special_bonds](special_bonds) command that includes all pairs in the
neighbor list.
::::

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

## Output info

This compute calculates a per-atom array with *nlvalues* columns, giving
the $Q_\ell$ values for each atom, which are real numbers in the range
$0 \le Q_\ell \le 1$.

If the keyword *wl* is set to yes, then the $W_\ell$ values for each
atom will be added to the output array, which are real numbers.

If the keyword *wl/hat* is set to yes, then the $\hat{W}_\ell$ values
for each atom will be added to the output array, which are real numbers.

If the keyword *components* is set, then the real and imaginary parts of
each component of *normalized* $\hat{Y}_{\ell m}$ will be added to the
output array in the following order: $\Re(\hat{Y}_{-m}),$
$\Im(\hat{Y}_{-m}),$ $\Re(\hat{Y}_{-m+1}),$
$\Im(\hat{Y}_{-m+1}), \dotsc,$ $\Re(\hat{Y}_m),$ $\Im(\hat{Y}_m).$

In summary, the per-atom array will contain *nlvalues* columns, followed
by an additional *nlvalues* columns if *wl* is set to yes, followed by
an additional *nlvalues* columns if *wl/hat* is set to yes, followed by
an additional 2\*(2\* *ldegree*+1) columns if the *components* keyword
is set.

These values can be accessed by any command that uses per-atom values
from a compute as input. See the [Howto output](Howto_output) doc page
for an overview of LAMMPS output options.

## Restrictions

none

## Related commands

[compute coord/atom](compute_coord_atom), [compute
centro/atom](compute_centro_atom), [compute
hexorder/atom](compute_hexorder_atom)

## Default

The option defaults are *cutoff* = pair style cutoff, *nnn* = 12,
*degrees* = 5 4 6 8 10 12 (i.e., $Q_4$, $Q_6$, $Q_8$, $Q_{10}$, and
$Q_{12}$), *wl* = no, *wl/hat* = no, *components* off, and *chunksize* =
16384

------------------------------------------------------------------------

::: {#Steinhardt}
**(Steinhardt)** P. Steinhardt, D. Nelson, and M. Ronchetti, Phys. Rev.
B 28, 784 (1983).
:::

::: {#Mickel}
**(Mickel)** W. Mickel, S. C. Kapfer, G. E. Schroeder-Turkand, K. Mecke,
J. Chem. Phys. 138, 044501 (2013).
:::

::: {#tenWolde2}
**(tenWolde)** P. R. ten Wolde, M. J. Ruiz-Montero, D. Frenkel, J. Chem.
Phys. 104, 9932 (1996).
:::
