# pair_style lj/long/coul/long command

Accelerator Variants: *lj/long/coul/long/intel*,
*lj/long/coul/long/omp*, *lj/long/coul/long/opt*

# pair_style lj/long/tip4p/long command

Accelerator Variants: *lj/long/tip4p/long/omp*

## Syntax

``` LAMMPS
pair_style style args
```

-   style = *lj/long/coul/long* or *lj/long/tip4p/long*
-   args = list of arguments for a particular style

<!-- -->

    *lj/long/coul/long* args = flag_lj flag_coul cutoff (cutoff2)
      flag_lj = *long* or *cut* or *off*
        *long* = use Kspace long-range summation for dispersion 1/r\^6 term
        *cut* = use a cutoff on dispersion 1/r\^6 term
        *off* = omit disperion 1/r\^6 term entirely
      flag_coul = *long* or *off*
        *long* = use Kspace long-range summation for Coulombic 1/r term
        *off* = omit Coulombic term
      cutoff = global cutoff for LJ (and Coulombic if only 1 arg) (distance units)
      cutoff2 = global cutoff for Coulombic (optional) (distance units)
    *lj/long/tip4p/long* args = flag_lj flag_coul otype htype btype atype qdist cutoff (cutoff2)
      flag_lj = *long* or *cut*
        *long* = use Kspace long-range summation for dispersion 1/r\^6 term
        *cut* = use a cutoff
      flag_coul = *long* or *off*
        *long* = use Kspace long-range summation for Coulombic 1/r term
        *off* = omit Coulombic term
      otype,htype = atom types for TIP4P O and H
      btype,atype = bond and angle types for TIP4P waters
      qdist = distance from O atom to massless charge (distance units)
      cutoff = global cutoff for LJ (and Coulombic if only 1 arg) (distance units)
      cutoff2 = global cutoff for Coulombic (optional) (distance units)

## Examples

``` LAMMPS
pair_style lj/long/coul/long cut off 2.5
pair_style lj/long/coul/long cut long 2.5 4.0
pair_style lj/long/coul/long long long 2.5 4.0
pair_coeff * * 1 1
pair_coeff 1 1 1 3 4

pair_style lj/long/tip4p/long long long 1 2 7 8 0.15 12.0
pair_style lj/long/tip4p/long long long 1 2 7 8 0.15 12.0 10.0
pair_coeff * * 100.0 3.0
pair_coeff 1 1 100.0 3.5 9.0
```

## Description

Style *lj/long/coul/long* computes the standard 12/6 Lennard-Jones
potential:

$$\begin{aligned}
E = 4 \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} -
                    \left(\frac{\sigma}{r}\right)^6 \right]
                    \qquad r < r_c \\
\end{aligned}$$

with $\epsilon$ and $\sigma$ being the usual Lennard-Jones potential
parameters, plus the Coulomb potential, given by:

$$E = \frac{C q_i q_j}{\epsilon  r} \qquad r < r_c$$

where C is an energy-conversion constant, $q_i$ and $q_j$ are the
charges on the 2 atoms, $\epsilon$ is the dielectric constant which can
be set by the [dielectric](dielectric) command, and $r_c$ is the cutoff.
If one cutoff is specified in the pair_style command, it is used for
both the LJ and Coulombic terms. If two cutoffs are specified, they are
used as cutoffs for the LJ and Coulombic terms respectively.

The purpose of this pair style is to capture long-range interactions
resulting from both attractive 1/r\^6 Lennard-Jones and Coulombic 1/r
interactions. This is done by use of the *flag_lj* and *flag_coul*
settings. The [In \'t Veld](Veld2) paper has more details on when it is
appropriate to include long-range 1/r\^6 interactions, using this
potential.

Style *lj/long/tip4p/long* implements the TIP4P water model of
[(Jorgensen)](Jorgensen4), which introduces a massless site located a
short distance away from the oxygen atom along the bisector of the HOH
angle. The atomic types of the oxygen and hydrogen atoms, the bond and
angle types for OH and HOH interactions, and the distance to the
massless charge site are specified as pair_style arguments.

:::: note
::: title
Note
:::

For each TIP4P water molecule in your system, the atom IDs for the O and
2 H atoms must be consecutive, with the O atom first. This is to enable
LAMMPS to \"find\" the 2 H atoms associated with each O atom. For
example, if the atom ID of an O atom in a TIP4P water molecule is 500,
then its 2 H atoms must have IDs 501 and 502.
::::

See the [Howto tip4p](Howto_tip4p) page for more information on how to
use the TIP4P pair style. Note that the neighbor list cutoff for Coulomb
interactions is effectively extended by a distance 2\*qdist when using
the TIP4P pair style, to account for the offset distance of the
fictitious charges on O atoms in water molecules. Thus it is typically
best in an efficiency sense to use a LJ cutoff \>= Coulombic cutoff +
2\*qdist, to shrink the size of the neighbor list. This leads to
slightly larger cost for the long-range calculation, so you can test the
trade-off for your model.

If *flag_lj* is set to *long*, no cutoff is used on the LJ 1/r\^6
dispersion term. The long-range portion can be calculated by using the
[kspace_style ewald/disp or pppm/disp](kspace_style) commands. The
specified LJ cutoff then determines which portion of the LJ interactions
are computed directly by the pair potential versus which part is
computed in reciprocal space via the Kspace style. If *flag_lj* is set
to *cut*, the LJ interactions are simply cutoff, as with [pair_style
lj/cut](pair_lj).

If *flag_coul* is set to *long*, no cutoff is used on the Coulombic
interactions. The long-range portion can calculated by using any of
several [kspace_style](kspace_style) command options such as *pppm* or
*ewald*. Note that if *flag_lj* is also set to long, then the
*ewald/disp* or *pppm/disp* Kspace style needs to be used to perform the
long-range calculations for both the LJ and Coulombic interactions. If
*flag_coul* is set to *off*, Coulombic interactions are not computed.

The following coefficients must be defined for each pair of atoms types
via the [pair_coeff](pair_coeff) command as in the examples above, or in
the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands, or by mixing as described below:

-   $\epsilon$ (energy units)
-   $\sigma$ (distance units)
-   cutoff1 (distance units)
-   cutoff2 (distance units)

Note that sigma is defined in the LJ formula as the zero-crossing
distance for the potential, not as the energy minimum at 2\^(1/6) sigma.

The latter 2 coefficients are optional. If not specified, the global LJ
and Coulombic cutoffs specified in the pair_style command are used. If
only one cutoff is specified, it is used as the cutoff for both LJ and
Coulombic interactions for this type pair. If both coefficients are
specified, they are used as the LJ and Coulombic cutoffs for this type
pair.

Note that if you are using *flag_lj* set to *long*, you cannot specify a
LJ cutoff for an atom type pair, since only one global LJ cutoff is
allowed. Similarly, if you are using *flag_coul* set to *long*, you
cannot specify a Coulombic cutoff for an atom type pair, since only one
global Coulombic cutoff is allowed.

For *lj/long/tip4p/long* only the LJ cutoff can be specified since a
Coulombic cutoff cannot be specified for an individual I,J type pair.
All type pairs use the same global Coulombic cutoff specified in the
pair_style command.

------------------------------------------------------------------------

A version of these styles with a soft core, *lj/cut/soft*, suitable for
use in free energy calculations, is part of the FEP package and is
documented with the [pair_style \*/soft](pair_fep_soft) styles. The
version with soft core is only available if LAMMPS was built with that
package. See the [Build package](Build_package) page for more info.

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

## Mixing, shift, table, tail correction, restart, rRESPA info

For atom type pairs I,J and I != J, the epsilon and sigma coefficients
and cutoff distance for all of the lj/long pair styles can be mixed. The
default mix value is *geometric*. See the \"pair_modify\" command for
details.

These pair styles support the [pair_modify](pair_modify) shift option
for the energy of the Lennard-Jones portion of the pair interaction,
assuming *flag_lj* is *cut*.

These pair styles support the [pair_modify](pair_modify) table and
table/disp options since they can tabulate the short-range portion of
the long-range Coulombic and dispersion interactions.

Thes pair styles do not support the [pair_modify](pair_modify) tail
option for adding a long-range tail correction to the Lennard-Jones
portion of the energy and pressure.

These pair styles write their information to [binary restart
files](restart), so pair_style and pair_coeff commands do not need to be
specified in an input script that reads a restart file.

The pair lj/long/coul/long styles support the use of the *inner*,
*middle*, and *outer* keywords of the [run_style respa](run_style)
command, meaning the pairwise forces can be partitioned by distance at
different levels of the rRESPA hierarchy. See the [run_style](run_style)
command for details.

------------------------------------------------------------------------

## Restrictions

These styles are part of the KSPACE package. They are only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

## Related commands

[pair_coeff](pair_coeff)

## Default

none

------------------------------------------------------------------------

::: {#Veld2}
**(In \'t Veld)** In \'t Veld, Ismail, Grest, J Chem Phys, 127, 144711
(2007).
:::

::: {#Jorgensen4}
**(Jorgensen)** Jorgensen, Chandrasekhar, Madura, Impey, Klein, J Chem
Phys, 79, 926 (1983).
:::
