# pair_style lj/charmm/coul/charmm command

Accelerator Variants: *lj/charmm/coul/charmm/gpu*,
*lj/charmm/coul/charmm/intel*, *lj/charmm/coul/charmm/kk*,
*lj/charmm/coul/charmm/omp*

# pair_style lj/charmm/coul/charmm/implicit command

Accelerator Variants: *lj/charmm/coul/charmm/implicit/kk*,
*lj/charmm/coul/charmm/implicit/omp*

# pair_style lj/charmm/coul/long command

Accelerator Variants: *lj/charmm/coul/long/gpu*,
*lj/charmm/coul/long/intel*, *lj/charmm/coul/long/kk*,
*lj/charmm/coul/long/opt*, *lj/charmm/coul/long/omp*

# pair_style lj/charmm/coul/msm command

Accelerator Variants: *lj/charmm/coul/msm/omp*

# pair_style lj/charmmfsw/coul/charmmfsh command

# pair_style lj/charmmfsw/coul/long command

## Syntax

``` LAMMPS
pair_style style args
```

-   style = *lj/charmm/coul/charmm* or *lj/charmm/coul/charmm/implicit*
    or *lj/charmm/coul/long* or *lj/charmm/coul/msm* or
    *lj/charmmfsw/coul/charmmfsh* or *lj/charmmfsw/coul/long*
-   args = list of arguments for a particular style

<!-- -->

    *lj/charmm/coul/charmm* args = inner outer (inner2) (outer2)
      inner, outer = global switching cutoffs for Lennard Jones (and Coulombic if only 2 args)
      inner2, outer2 = global switching cutoffs for Coulombic (optional)
    *lj/charmm/coul/charmm/implicit* args = inner outer (inner2) (outer2)
      inner, outer = global switching cutoffs for LJ (and Coulombic if only 2 args)
      inner2, outer2 = global switching cutoffs for Coulombic (optional)
    *lj/charmm/coul/long* args = inner outer (cutoff)
      inner, outer = global switching cutoffs for LJ (and Coulombic if only 2 args)
      cutoff = global cutoff for Coulombic (optional, outer is Coulombic cutoff if only 2 args)
    *lj/charmm/coul/msm* args = inner outer (cutoff)
      inner, outer = global switching cutoffs for LJ (and Coulombic if only 2 args)
      cutoff = global cutoff for Coulombic (optional, outer is Coulombic cutoff if only 2 args)
    *lj/charmmfsw/coul/charmmfsh* args = inner outer (cutoff)
      inner, outer = global cutoffs for LJ (and Coulombic if only 2 args)
      cutoff = global cutoff for Coulombic (optional, outer is Coulombic cutoff if only 2 args)
    *lj/charmmfsw/coul/long* args = inner outer (cutoff)
      inner, outer = global cutoffs for LJ (and Coulombic if only 2 args)
      cutoff = global cutoff for Coulombic (optional, outer is Coulombic cutoff if only 2 args)

## Examples

``` LAMMPS
pair_style lj/charmm/coul/charmm 8.0 10.0
pair_style lj/charmm/coul/charmm 8.0 10.0 7.0 9.0
pair_style lj/charmmfsw/coul/charmmfsh 10.0 12.0
pair_style lj/charmmfsw/coul/charmmfsh 10.0 12.0 9.0
pair_coeff * * 100.0 2.0
pair_coeff 1 1 100.0 2.0 150.0 3.5

pair_style lj/charmm/coul/charmm/implicit 8.0 10.0
pair_style lj/charmm/coul/charmm/implicit 8.0 10.0 7.0 9.0
pair_coeff * * 100.0 2.0
pair_coeff 1 1 100.0 2.0 150.0 3.5

pair_style lj/charmm/coul/long 8.0 10.0
pair_style lj/charmm/coul/long 8.0 10.0 9.0
pair_style lj/charmmfsw/coul/long 8.0 10.0
pair_style lj/charmmfsw/coul/long 8.0 10.0 9.0
pair_coeff * * 100.0 2.0
pair_coeff 1 1 100.0 2.0 150.0 3.5

pair_style lj/charmm/coul/msm 8.0 10.0
pair_style lj/charmm/coul/msm 8.0 10.0 9.0
pair_coeff * * 100.0 2.0
pair_coeff 1 1 100.0 2.0 150.0 3.5
```

## Description

These pair styles compute Lennard Jones (LJ) and Coulombic interactions
with additional switching or shifting functions that ramp the energy
and/or force smoothly to zero between an inner and outer cutoff. They
are implementations of the widely used CHARMM force field used in the
[CHARMM](https://www.charmm.org)\_ MD code (and others). See
[(MacKerell)](pair-MacKerell) for a description of the CHARMM force
field.

The styles with *charmm* (not *charmmfsw* or *charmmfsh*) in their name
are the older, original LAMMPS implementations. They compute the LJ and
Coulombic interactions with an energy switching function (esw, shown in
the formula below as S(r)), which ramps the energy smoothly to zero
between the inner and outer cutoff. This can cause irregularities in
pairwise forces (due to the discontinuous second derivative of energy at
the boundaries of the switching region), which in some cases can result
in detectable artifacts in an MD simulation.

The newer styles with *charmmfsw* or *charmmfsh* in their name replace
the energy switching with force switching (fsw) and force shifting (fsh)
functions, for LJ and Coulombic interactions respectively. These follow
the formulas and description given in [(Steinbach)](Steinbach) and
[(Brooks)](Brooks1) to minimize these artifacts.

:::: note
::: title
Note
:::

The newer *charmmfsw* or *charmmfsh* styles were released in March 2017.
We recommend they be used instead of the older *charmm* styles. This
includes the newer [dihedral_style charmmfsw](dihedral_charmm) command.
Eventually code from the new styles will propagate into the related pair
styles (e.g. implicit, accelerator, free energy variants).
::::

:::: note
::: title
Note
:::

The newest CHARMM pair styles reset the Coulombic energy conversion
factor used internally in the code, from the LAMMPS value to the CHARMM
value, as if it were effectively a parameter of the force field. This is
because the CHARMM code uses a slightly different value for the this
conversion factor in [real units](units) (kcal/mol), namely CHARMM =
332.0716, LAMMPS = 332.06371. This is to enable more precise agreement
by LAMMPS with the CHARMM force field energies and forces, when using
one of these two CHARMM pair styles.
::::

$$\begin{aligned}
E = & LJ(r) \qquad \qquad \qquad r < r_{\rm in} \\
  = & S(r) * LJ(r) \qquad \qquad r_{\rm in} < r < r_{\rm out} \\
  = & 0 \qquad \qquad \qquad \qquad r > r_{\rm out} \\
E = & C(r) \qquad \qquad \qquad r < r_{\rm in} \\
  = & S(r) * C(r) \qquad \qquad r_{\rm in} < r < r_{\rm out} \\
  = & 0 \qquad \qquad \qquad \qquad r > r_{\rm out} \\
LJ(r) = & 4 \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} -
        \left(\frac{\sigma}{r}\right)^6 \right] \\
C(r) = & \frac{C q_i q_j}{ \epsilon r} \\
S(r) = & \frac{ \left[r_{\rm out}^2 - r^2\right]^2
  \left[r_{\rm out}^2 + 2r^2 - 3{r_{\rm in}^2}\right]}
{ \left[r_{\rm out}^2 - {r_{\rm in}}^2\right]^3 }
\end{aligned}$$

where S(r) is the energy switching function mentioned above for the
*charmm* styles. See the [(Steinbach)](Steinbach) paper for the
functional forms of the force switching and force shifting functions
used in the *charmmfsw* and *charmmfsh* styles.

When using the *lj/charmm/coul/charmm styles*, both the LJ and Coulombic
terms require an inner and outer cutoff. They can be the same for both
formulas or different depending on whether 2 or 4 arguments are used in
the pair_style command. For the *lj/charmmfsw/coul/charmmfsh* style, the
LJ term requires both an inner and outer cutoff, while the Coulombic
term requires only one cutoff. If the Coulombic cutoff is not specified
(2 instead of 3 arguments), the LJ outer cutoff is used for the
Coulombic cutoff. In all cases where an inner and outer cutoff are
specified, the inner cutoff distance must be less than the outer cutoff.
It is typical to make the difference between the inner and outer cutoffs
about 2.0 Angstroms.

Style *lj/charmm/coul/charmm/implicit* computes the same formulas as
style *lj/charmm/coul/charmm* except that an additional 1/r term is
included in the Coulombic formula. The Coulombic energy thus varies as
1/r\^2. This is effectively a distance-dependent dielectric term which
is a simple model for an implicit solvent with additional screening. It
is designed for use in a simulation of an unsolvated biomolecule (no
explicit water molecules).

Styles *lj/charmm/coul/long* and *lj/charmm/coul/msm* compute the same
formulas as style *lj/charmm/coul/charmm* and style
*lj/charmmfsw/coul/long* computes the same formulas as style
*lj/charmmfsw/coul/charmmfsh*, except that an additional damping factor
is applied to the Coulombic term, so it can be used in conjunction with
the [kspace_style](kspace_style) command and its *ewald* or *pppm* or
*msm* option. Only one Coulombic cutoff is specified for these styles;
if only 2 arguments are used in the pair_style command, then the outer
LJ cutoff is used as the single Coulombic cutoff. The Coulombic cutoff
specified for these styles means that pairwise interactions within this
distance are computed directly; interactions outside that distance are
computed in reciprocal space.

The following coefficients must be defined for each pair of atoms types
via the [pair_coeff](pair_coeff) command as in the examples above, or in
the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands, or by mixing as described below:

-   $\epsilon$ (energy units)
-   $\sigma$ (distance units)
-   $\epsilon_{14}$ (energy units)
-   $\sigma_{14}$ (distance units)

Note that $\sigma$ is defined in the LJ formula as the zero-crossing
distance for the potential, not as the energy minimum at
$2^{1/6} \sigma$.

The latter 2 coefficients are optional. If they are specified, they are
used in the LJ formula between 2 atoms of these types which are also
first and fourth atoms in any dihedral. No cutoffs are specified because
the CHARMM force field does not allow varying cutoffs for individual
atom pairs; all pairs use the global cutoff(s) specified in the
pair_style command.

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

For atom type pairs I,J and I != J, the epsilon, sigma, epsilon_14, and
sigma_14 coefficients for all of the lj/charmm pair styles can be mixed.
The default mix value is *arithmetic* to coincide with the usual
settings for the CHARMM force field. See the \"pair_modify\" command for
details.

None of the *lj/charmm* or *lj/charmmfsw* pair styles support the
[pair_modify](pair_modify) shift option, since the Lennard-Jones portion
of the pair interaction is smoothed to 0.0 at the cutoff.

The *lj/charmm/coul/long* and *lj/charmmfsw/coul/long* styles support
the [pair_modify](pair_modify) table option since they can tabulate the
short-range portion of the long-range Coulombic interaction.

None of the *lj/charmm* or *lj/charmmfsw* pair styles support the
[pair_modify](pair_modify) tail option for adding long-range tail
corrections to energy and pressure, since the Lennard-Jones portion of
the pair interaction is smoothed to 0.0 at the cutoff.

All of the *lj/charmm* and *lj/charmmfsw* pair styles write their
information to [binary restart files](restart), so pair_style and
pair_coeff commands do not need to be specified in an input script that
reads a restart file.

The *lj/charmm/coul/long* and *lj/charmmfsw/coul/long* pair styles
support the use of the *inner*, *middle*, and *outer* keywords of the
[run_style respa](run_style) command, meaning the pairwise forces can be
partitioned by distance at different levels of the rRESPA hierarchy. The
other styles only support the *pair* keyword of run_style respa. See the
[run_style](run_style) command for details.

------------------------------------------------------------------------

## Restrictions

All the styles with *coul/charmm* or *coul/charmmfsh* styles are part of
the MOLECULE package. All the styles with *coul/long* style are part of
the KSPACE package. They are only enabled if LAMMPS was built with those
packages. See the [Build package](Build_package) doc page for more info.

## Related commands

[pair_coeff](pair_coeff)

## Default

none

------------------------------------------------------------------------

::: {#Brooks1}
**(Brooks)** Brooks, et al, J Comput Chem, 30, 1545 (2009).
:::

::: {#pair-MacKerell}
**(MacKerell)** MacKerell, Bashford, Bellott, Dunbrack, Evanseck, Field,
Fischer, Gao, Guo, Ha, et al, J Phys Chem, 102, 3586 (1998).
:::

::: {#Steinbach}
**(Steinbach)** Steinbach, Brooks, J Comput Chem, 15, 667 (1994).
:::
