# pair_style lj/expand/sphere command

Accelerator Variant: *lj/expand/sphere/omp*

## Syntax

``` LAMMPS
pair_style style args
```

-   style = *lj/expand/sphere*
-   args = list of arguments for a particular style

<!-- -->

    *lj/expand/sphere* args = cutoff
      cutoff = global cutoff for Lennard Jones interactions (distance units)

## Examples

``` LAMMPS
pair_style lj/expand/sphere 2.5
pair_coeff * * 1.0 1.0
pair_coeff 1 1 1.1 0.4 2.8
```

## Description

::: versionadded
15Jun2023
:::

The *lj/expand/sphere* style compute a 12/6 Lennard-Jones potential with
a distance shifted by $\Delta = \frac{1}{2} (d_i + d_j)$, the average
diameter of both atoms. This can be used to model particles of different
sizes but same interactions, which is different from using different
sigma values as in [pair style lj/cut/sphere](pair_lj_cut_sphere).

$$E = 4 \epsilon \left[ \left(\frac{\sigma}{r - \Delta}\right)^{12} -
  \left(\frac{\sigma}{r - \Delta}\right)^6 \right]
  \qquad r < r_c + \Delta$$

$r_c$ is the cutoff which does not include the distance $\Delta$. I.e.
the actual force cutoff is the sum $r_c + \Delta$.

This is the same potential function used by the [lj/expand](pair_lj)
pair style, but the $\Delta$ parameter is not set as a per-type
parameter via the [pair_coeff command](pair_coeff). Instead it is
calculated individually for each pair using the per-atom diameter
attribute of [atom_style sphere](atom_style) for the two atoms as the
average diameter, $\Delta = \frac{1}{2} (d_i + d_j)$

Note that $\sigma$ is defined in the LJ formula above as the
zero-crossing distance for the potential, *not* as the energy minimum
which is at $2^{\frac{1}{6}} \sigma$.

::: {.admonition .note}
Notes on cutoffs, neighbor lists, and efficiency

If your system is mildly polydisperse, meaning the ratio of the diameter
of the largest particle to the smallest is less than 2.0, then the
neighbor lists built by the code should be reasonably efficient. Which
means they will not contain too many particle pairs that do not
interact. However, if your system is highly polydisperse (ratio \> 2.0),
the neighbor list build and force computations may be inefficient. There
are two ways to try and speed up the simulations.

The first is to assign atoms to different atom types so that atoms of
each type are similar in size. E.g. if particle diameters range from 1
to 5 use 4 atom types, ensuring atoms of type 1 have diameters from
1.0-2.0, type 2 from 2.0-3.0, etc. This will reduce the number of
non-interacting pairs in the neighbor lists and thus reduce the time
spent on computing pairwise interactions.

The second is to use the [neighbor multi](neighbor) command which
enabled a different algorithm for building neighbor lists. This will
also require that you assign multiple atom types according to diameters,
but will in addition use a more efficient size-dependent strategy to
construct the neighbor lists and thus reduce the time spent on building
neighbor lists.

Here are example input script commands using the first option for a
highly polydisperse system:

``` c++
units           lj
atom_style      sphere
lattice         fcc 0.8442
region          box block 0 10 0 10 0 10
create_box      2 box
create_atoms    1 box

# create atoms with random diameters from bimodal distribution
variable switch atom random(0.0,1.0,345634)
variable diam atom (v_switch<0.75)*normal(0.2,0.04,325)+(v_switch>=0.7)*normal(0.6,0.2,453)
set group all diameter v_diam

# assign type 2 to atoms with diameter > 0.35
variable large atom (2.0*radius)>0.35
group large variable large
set group large type 2

pair_style      lj/expand/sphere 2.0
pair_coeff      * * 1.0 0.5

neighbor 0.3 bin
```

Using multiple atom types speeds up the calculation for this example by
more than 30 percent, but using the multi-style neighbor list does not
provide a speedup.
:::

## Coefficients

The following coefficients must be defined for each pair of atoms types
via the [pair_coeff](pair_coeff) command as in the examples above, or in
the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands, or by mixing as described below:

-   $\epsilon$ (energy units)
-   $\sigma$ (distance units)
-   LJ cutoff (distance units) (optional)

The last coefficient is optional. If not specified, the global LJ cutoff
specified in the [pair_style command](pair_style) is used.

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

For atom type pairs I,J and I != J, the epsilon, sigma, and cutoff
coefficients for the *lj/expand/sphere* pair style can be mixed. The
default mixing style is *geometric*. See the [pair_modify
command](pair_modify) for details.

The *lj/expand/sphere* pair style supports the [pair_modify
shift](pair_modify) option for the energy of the Lennard-Jones portion
of the pair interaction.

The *lj/expand/sphere* pair style does *not* support the
[pair_modify](pair_modify) tail option for adding a long-range tail
corrections to the energy and pressure.

The *lj/expand/sphere* pair style writes its information to [binary
restart files](restart), so pair_style and pair_coeff commands do not
need to be specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
[run_style respa](run_style) command. It does *not* support the *inner*,
*middle*, *outer* keywords.

------------------------------------------------------------------------

## Restrictions

The *lj/expand/sphere* pair style is only enabled if LAMMPS was built
with the EXTRA-PAIR package. See the [Build package](Build_package) page
for more info.

------------------------------------------------------------------------

## Related commands

-   [pair_coeff](pair_coeff)
-   [pair_style lj/cut](pair_lj)
-   [pair_style lj/cut/sphere](pair_lj_cut_sphere)

## Default

none
