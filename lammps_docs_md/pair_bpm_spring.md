# pair_style bpm/spring command

## Syntax

``` LAMMPS
pair_style bpm/spring
```

## Examples

``` LAMMPS
pair_style bpm/spring
pair_coeff * * 1.0 1.0 1.0
pair_coeff 1 1 1.0 1.0 1.0
```

## Description

::: versionadded
4May2022
:::

Style *bpm/spring* computes pairwise forces with the formula

$$F = k (r - r_c)$$

where $k$ is a stiffness and $r_c$ is the cutoff length. An additional
damping force is also applied to interacting particles. The force is
proportional to the difference in the normal velocity of particles

$$F_D = - \gamma w (\hat{r} \bullet \vec{v})$$

where $\gamma$ is the damping strength, $\hat{r}$ is the radial normal
vector, $\vec{v}$ is the velocity difference between the two particles,
and $w$ is a smoothing factor. This smoothing factor is constructed such
that damping forces go to zero as particles come out of contact to avoid
discontinuities. It is given by

$$w = 1.0 - \left( \frac{r}{r_c} \right)^8 .$$

This pair style is designed for use in a spring-based bonded particle
model. It mirrors the construction of the [bpm/spring](bond_bpm_spring)
bond style.

This pair interaction is always applied to pairs of non-bonded particles
that are within the interaction distance. For pairs of bonded particles
that are within the interaction distance, there is the option to either
include this pair interaction and overlay the pair force over the bond
force or to exclude this pair interaction such that the two particles
only interact via the bond force. See discussion of the *overlay/pair*
option for BPM bond styles and the [special_bonds](special_bonds)
command in the [how to](Howto_bpm) page on BPMs for more details.

The following coefficients must be defined for each pair of atom types
via the [pair_coeff](pair_coeff) command as in the examples above, or in
the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands, or by mixing as described below:

-   $k$ (force/distance units)
-   $r_c$ (distance units)
-   $\gamma$ (force/velocity units)

------------------------------------------------------------------------

## Mixing, shift, table, tail correction, restart, rRESPA info

For atom type pairs I,J and I != J, the A coefficient and cutoff
distance for this pair style can be mixed. A is always mixed via a
*geometric* rule. The cutoff is mixed according to the pair_modify mix
value. The default mix value is *geometric*. See the \"pair_modify\"
command for details.

This pair style does not support the [pair_modify](pair_modify) shift
option, since the pair interaction goes to 0.0 at the cutoff.

The [pair_modify](pair_modify) table and tail options are not relevant
for this pair style.

This pair style writes its information to [binary restart
files](restart), so pair_style and pair_coeff commands do not need to be
specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
[run_style respa](run_style) command. It does not support the *inner*,
*middle*, *outer* keywords.

------------------------------------------------------------------------

## Restrictions

This pair style is part of the BPM package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

## Related commands

[pair_coeff](pair_coeff), [bond bpm/spring](bond_bpm_spring)

## Default

none
