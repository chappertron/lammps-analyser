# pair_style local/density command

## Syntax

``` LAMMPS
pair_style style arg
```

-   style = *local/density*
-   arg = name of file containing tabulated values of local density and
    the potential

## Examples

``` LAMMPS
pair_style local/density benzene_water.localdensity.table

pair_style hybrid/overlay table spline 500 local/density
pair_coeff * * local/density  benzene_water.localdensity.table
```

## Description

The local density (LD) potential is a mean-field manybody potential,
and, in some way, a generalization of embedded atom models (EAM). The
name \"local density potential\" arises from the fact that it assigns an
energy to an atom depending on the number of neighboring atoms of a
given type around it within a predefined spherical volume (i.e., within
the cutoff). The bottom-up coarse-graining (CG) literature suggests that
such potentials can be widely useful in capturing effective multibody
forces in a computationally efficient manner and thus improve the
quality of CG models of implicit solvation [(Sanyal1)](Sanyal1) and
phase-segregation in liquid mixtures [(Sanyal2)](Sanyal2), and provide
guidelines to determine the extent of manybody correlations present in a
CG model [(Rosenberger)](Rosenberger). The LD potential in LAMMPS is
primarily intended to be used as a corrective potential over traditional
pair potentials in bottom-up CG models via [hybrid/overlay pair
style](pair_hybrid) with other explicit pair interaction terms (e.g.,
tabulated, Lennard-Jones, Morse etc.). Because the LD potential is not a
pair potential per se, it is implemented simply as a single auxiliary
file with all specifications that will be read upon initialization.

:::: note
::: title
Note
:::

Thus when used as the only interaction in the system, there is no
corresponding pair_coeff command and when used with other pair styles
using the hybrid/overlay option, the corresponding pair_coeff command
must be supplied \* \* as placeholders for the atom types.
::::

------------------------------------------------------------------------

**System with a single CG atom type:**

A system of a single atom type (e.g., LJ argon) with a single local
density (LD) potential would have an energy given by:

$$U_{LD} = \sum_i F(\rho_i)$$

where $\rho_i$ is the LD at atom *i* and $F(\rho)$ is similar in spirit
to the embedding function used in EAM potentials. The LD at atom *i* is
given by the sum

$$\rho_i = \sum_{j \neq i} \varphi(r_{ij})$$

where $\varphi$ is an indicator function that is one at r=0 and zero
beyond a cutoff distance R2. The choice of the functional form of
$\varphi$ is somewhat arbitrary, but the following piecewise cubic
function has proven sufficiently general: [(Sanyal1)](Sanyal1),
[(Sanyal2)](Sanyal2) [(Rosenberger)](Rosenberger)

$$\begin{aligned}
\varphi(r) =
\begin{cases}
1 & r \le R_1 \\
c_0 + c_2r^2 + c_4r^4 + c_6r^6  & r \in (R_1, R_2) \\
0 & r \ge R_2
\end{cases}
\end{aligned}$$

The constants *c* are chosen so that the indicator function smoothly
interpolates between 1 and 0 between the distances R1 and R2, which are
called the inner and outer cutoffs, respectively. Thus phi satisfies
phi(R1) = 1, phi(R2) = dphi/dr @ (r=R1) = dphi/dr @ (r=R2) = 0. The
embedding function F(rho) may or may not have a closed-form expression.
To maintain generality, it is practically represented with a
spline-interpolated table over a predetermined range of rho. Outside of
that range it simply adopts zero values at the endpoints.

It can be shown that the total force between two atoms due to the LD
potential takes the form of a pair force, which motivates its
designation as a LAMMPS pair style. Please see [(Sanyal1)](Sanyal1) for
details of the derivation.

------------------------------------------------------------------------

**Systems with arbitrary numbers of atom types:**

The potential is easily generalized to systems involving multiple atom
types:

$$U_{LD} = \sum_i a_\alpha F(\rho_i)$$

with the LD expressed as

$$\rho_i = \sum_{j \neq i} b_\beta \varphi(r_{ij})$$

where $\alpha$ gives the type of atom *i*, $\beta$ the type of atom *j*,
and the coefficients *a* and *b* filter for atom types as specified by
the user. *a* is called the central atom filter as it determines to
which atoms the potential applies; $a_{\alpha} =
1$ if the LD potential applies to atom type $\alpha$ else zero. On the
other hand, *b* is called the neighbor atom filter because it specifies
which atom types to use in the calculation of the LD; $b_{\beta} =
1$ if atom type $\beta$ contributes to the LD and zero otherwise.

:::: note
::: title
Note
:::

Note that the potentials need not be symmetric with respect to atom
types, which is the reason for two distinct sets of coefficients *a* and
*b*. An atom type may contribute to the LD but not the potential, or to
the potential but not the LD. Such decisions are made by the user and
should (ideally) be motivated on physical grounds for the problem at
hand.
::::

------------------------------------------------------------------------

**General form for implementation in LAMMPS:**

Of course, a system with many atom types may have many different
possible LD potentials, each with their own atom type filters, cutoffs,
and embedding functions. The most general form of this potential as
implemented in the pair_style local/density is:

$$U_{LD} = \sum_k U_{LD}^{(k)} = \sum_i \left[ \sum_k a_\alpha^{(k)} F^{(k)} \left(\rho_i^{(k)}\right) \right]$$

where, *k* is an index that spans the (arbitrary) number of applied LD
potentials N_LD. Each LD is calculated as before with:

$$\rho_i^{(k)} = \sum_j b_\beta^{(k)} \varphi^{(k)} (r_{ij})$$

The superscript on the indicator function phi simply indicates that it
is associated with specific values of the cutoff distances R1(k) and
R2(k). In summary, there may be N_LD distinct LD potentials. With each
potential type (k), one must specify:

-   the inner and outer cutoffs as R1 and R2
-   the central type filter a(k), where k = 1,2,\...N_LD
-   the neighbor type filter b(k), where k = 1,2,\...N_LD
-   the LD potential function F(k)(rho), typically as a table that is
    later spline-interpolated

------------------------------------------------------------------------

**Tabulated input file format:**

    Line 1:             comment or blank (ignored)
    Line 2:             comment or blank (ignored)
    Line 3:             N_LD N_rho (# of LD potentials and # of tabulated values, single space separated)
    Line 4:             blank (ignored)
    Line 5:             R1(k) R2(k) (lower and upper cutoffs, single space separated)
    Line 6:             central-types (central atom types, single space separated)
    Line 7:             neighbor-types (neighbor atom types single space separated)
    Line 8:             rho_min rho_max drho (min, max and diff. in tabulated rho values, single space separated)
    Line 9:             F(k)(rho_min + 0.drho)
    Line 10:            F(k)(rho_min + 1.drho)
    Line 11:            F(k)(rho_min + 2.drho)
    ...
    Line 9+N_rho:       F(k)(rho_min + N_rho . drho)
    Line 10+N_rho:      blank (ignored)

    Block 2

    Block 3

    Block N_LD

Lines 5 to 9+N_rho constitute the first block. Thus the input file is
separated (by blank lines) into N_LD blocks each representing a separate
LD potential and each specifying its own upper and lower cutoffs,
central and neighbor atoms, and potential. In general, blank lines
anywhere are ignored.

------------------------------------------------------------------------

## Mixing, shift, table, tail correction, restart, rRESPA info

This pair style does not support automatic mixing. For atom type pairs
$\alpha$, $\beta$ and $\alpha$ != $\beta$, even if LD potentials of type
($\alpha$, $\alpha$) and ($\beta$, $\beta$) are provided, you will need
to explicitly provide LD potential types ($\alpha$, $\beta$) and
($\beta$, $\alpha$) if need be (Here, the notation ($\alpha$, $\beta$)
means that $\alpha$ is the central atom to which the LD potential is
applied and $\beta$ is the neighbor atom which contributes to the LD
potential on $\alpha$).

This pair style does not support the [pair_modify](pair_modify) shift,
table, and tail options.

The local/density pair style does not write its information to [binary
restart files](restart), since it is stored in tabulated potential
files. Thus, you need to re-specify the pair_style and pair_coeff
commands in an input script that reads a restart file.

------------------------------------------------------------------------

## Restrictions

The local/density pair style is a part of the MANYBODY package. It is
only enabled if LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

## Related commands

[pair_coeff](pair_coeff)

## Default

none

------------------------------------------------------------------------

:::: {#Sanyal1}
::: {#Sanyal2}
**(Sanyal1)** Sanyal and Shell, Journal of Chemical Physics, 2016, 145
(3), 034109.
:::
::::

**(Sanyal2)** Sanyal and Shell, Journal of Physical Chemistry B, 122
(21), 5678-5693.

::: {#Rosenberger}
**(Rosenberger)** Rosenberger, Sanyal, Shell and van der Vegt, Journal
of Chemical Physics, 2019, 151 (4), 044111.
:::
