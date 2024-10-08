# dihedral_style charmm command

Accelerator Variants: *charmm/intel*, *charmm/kk*, *charmm/omp*

# dihedral_style charmmfsw command

## Syntax

``` LAMMPS
dihedral_style style
```

-   style = *charmm* or *charmmfsw*

## Examples

``` LAMMPS
dihedral_style charmm
dihedral_style charmmfsw
dihedral_coeff  1 0.2 1 180 1.0
dihedral_coeff  2 1.8 1   0 1.0
dihedral_coeff  1 3.1 2 180 0.5
```

## Description

The *charmm* and *charmmfsw* dihedral styles use the potential

$$E = K [ 1 + \cos (n \phi - d) ]$$

See [(MacKerell)](dihedral-MacKerell) for a description of the CHARMM
force field. This dihedral style can also be used for the AMBER force
field (see comment on weighting factors below). See
[(Cornell)](dihedral-Cornell) for a description of the AMBER force
field.

:::: note
::: title
Note
:::

The newer *charmmfsw* style was released in March 2017. We recommend it
be used instead of the older *charmm* style when running a simulation
with the CHARMM force field, either with long-range Coulombics or a
Coulombic cutoff, via the [pair_style
lj/charmmfsw/coul/long](pair_charmm) and [pair_style
lj/charmmfsw/coul/charmmfsh](pair_charmm) commands respectively.
Otherwise the older *charmm* style is fine to use. See the discussion
below and more details on the [pair_style charmm](pair_charmm) doc page.
::::

The following coefficients must be defined for each dihedral type via
the [dihedral_coeff](dihedral_coeff) command as in the example above, or
in the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   $K$ (energy)
-   $n$ (integer \>= 0)
-   $d$ (integer value of degrees)
-   weighting factor (1.0, 0.5, or 0.0)

The weighting factor is required to correct for double counting pairwise
non-bonded Lennard-Jones interactions in cyclic systems or when using
the CHARMM dihedral style with non-CHARMM force fields. With the CHARMM
dihedral style, interactions between the first and fourth atoms in a
dihedral are skipped during the normal non-bonded force computation and
instead evaluated as part of the dihedral using special epsilon and
sigma values specified with the [pair_coeff](pair_charmm) command of
pair styles that contain \"lj/charmm\" (e.g. [pair_style
lj/charmm/coul/long](pair_charmm)) In 6-membered rings, the same 1-4
interaction would be computed twice (once for the clockwise 1-4 pair in
dihedral 1-2-3-4 and once in the counterclockwise dihedral 1-6-5-4) and
thus the weighting factor has to be 0.5 in this case. In 4-membered or
5-membered rings, the 1-4 dihedral also is counted as a 1-2 or 1-3
interaction when going around the ring in the opposite direction and
thus the weighting factor is 0.0, as the 1-2 and 1-3 exclusions take
precedence.

Note that this dihedral weighting factor is unrelated to the scaling
factor specified by the [special bonds](special_bonds) command which
applies to all 1-4 interactions in the system. For CHARMM force fields,
the special_bonds 1-4 interaction scaling factor should be set to 0.0.
Since the corresponding 1-4 non-bonded interactions are computed with
the dihedral. This means that if any of the weighting factors defined as
dihedral coefficients (fourth coeff above) are non-zero, then you must
use a pair style with \"lj/charmm\" and set the special_bonds 1-4
scaling factor to 0.0 (which is the default). Otherwise 1-4 non-bonded
interactions in dihedrals will be computed twice.

For simulations using the CHARMM force field with a Coulombic cutoff,
the difference between the *charmm* and *charmmfsw* styles is in the
computation of the 1-4 non-bond interactions, though only if the
distance between the two atoms is within the switching region of the
pairwise potential defined by the corresponding CHARMM pair style, i.e.
within the outer cutoff specified for the pair style. The *charmmfsw*
style should only be used when using the corresponding [pair_style
lj/charmmfsw/coul/charmmfsw](pair_charmm) or [pair_style
lj/charmmfsw/coul/long](pair_charmm) commands. Use the *charmm* style
with the older [pair_style](pair_charmm) commands that have just
\"charmm\" in their style name. See the discussion on the [CHARMM
pair_style](pair_charmm) page for details.

Note that for AMBER force fields, which use pair styles with \"lj/cut\",
the special_bonds 1-4 scaling factor should be set to the AMBER defaults
(1/2 and 5/6) and all the dihedral weighting factors (fourth coeff
above) must be set to 0.0. In this case, you can use any pair style you
wish, since the dihedral does not need any Lennard-Jones parameter
information and will not compute any 1-4 non-bonded interactions.
Likewise the *charmm* or *charmmfsw* styles are identical in this case
since no 1-4 non-bonded interactions are computed.

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

When using run_style [respa](run_style), these dihedral styles must be
assigned to the same r-RESPA level as *pair* or *outer*.

When used in combination with CHARMM pair styles, the 1-4
[special_bonds](special_bonds) scaling factors must be set to 0.0.
Otherwise non-bonded contributions for these 1-4 pairs will be computed
multiple times.

These dihedral styles can only be used if LAMMPS was built with the
MOLECULE package. See the [Build package](Build_package) doc page for
more info.

## Related commands

[dihedral_coeff](dihedral_coeff)

## Default

none

------------------------------------------------------------------------

::: {#dihedral-Cornell}
**(Cornell)** Cornell, Cieplak, Bayly, Gould, Merz, Ferguson,
Spellmeyer, Fox, Caldwell, Kollman, JACS 117, 5179-5197 (1995).
:::

::: {#dihedral-MacKerell}
**(MacKerell)** MacKerell, Bashford, Bellott, Dunbrack, Evanseck, Field,
Fischer, Gao, Guo, Ha, et al, J Phys Chem B, 102, 3586 (1998).
:::
