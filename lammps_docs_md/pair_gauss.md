# pair_style gauss command

Accelerator Variants: *gauss/gpu*, *gauss/omp*

# pair_style gauss/cut command

Accelerator Variants: *gauss/cut/omp*

## Syntax

``` LAMMPS
pair_style gauss cutoff
pair_style gauss/cut cutoff
```

-   cutoff = global cutoff for Gauss interactions (distance units)

## Examples

``` LAMMPS
pair_style gauss 12.0
pair_coeff * * 1.0 0.9
pair_coeff 1 4 1.0 0.9 10.0

pair_style gauss/cut 3.5
pair_coeff 1 4 0.2805 1.45 0.112
```

## Description

Style *gauss* computes a tethering potential of the form

$$E = - A \exp(-B r^2) \qquad r < r_c$$

between an atom and its corresponding tether site which will typically
be a frozen atom in the simulation. $r_c$ is the cutoff.

The following coefficients must be defined for each pair of atom types
via the [pair_coeff](pair_coeff) command as in the examples above, or in
the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   A (energy units)
-   B (1/distance\^2 units)
-   cutoff (distance units)

The last coefficient is optional. If not specified, the global cutoff is
used.

Style *gauss/cut* computes a generalized Gaussian interaction potential
between pairs of particles:

$$E = \frac{H}{\sigma_h\sqrt{2\pi}} \exp\left[-\frac{(r-r_{mh})^2}{2\sigma_h^2}\right]$$

where H determines together with the standard deviation $\sigma_h$ the
peak height of the Gaussian function, and $r_{mh}$ the peak position.
Examples of the use of the Gaussian potentials include implicit solvent
simulations of salt ions [(Lenart)](Lenart2) and of surfactants
[(Jusufi)](Jusufi2). In these instances the Gaussian potential mimics
the hydration barrier between a pair of particles. The hydration barrier
is located at $r_{mh}$ and has a width of $\sigma_h$. The prefactor
determines the height of the potential barrier.

The following coefficients must be defined for each pair of atom types
via the [pair_coeff](pair_coeff) command as in the example above, or in
the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   H (energy \* distance units)
-   $r_{mh}$ (distance units)
-   $\sigma_h$ (distance units)
-   cutoff (distance units)

The last coefficient is optional. If not specified, the global cutoff is
used.

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

For atom type pairs I,J and I != J, the A, B, H, sigma_h, r_mh
parameters, and the cutoff distance for these pair styles can be mixed:

-   A (energy units)
-   $\sqrt{\frac{1}{B}}$ (distance units, see below)
-   H (energy units)
-   $r_{mh}$ (distance units)
-   $\sigma_h$ (distance units)
-   cutoff (distance units)

The default mix value is *geometric*. Only *arithmetic* and *geometric*
mix values are supported. See the \"pair_modify\" command for details.

The A and H parameters are mixed using the same rules normally used to
mix the \"epsilon\" parameter in a Lennard Jones interaction. The
sigma_h, r_mh, and the cutoff distance are mixed using the same rules
used to mix the \"sigma\" parameter in a Lennard Jones interaction. The
B parameter is converted to a distance (sigma), before mixing (using
sigma=B\^-0.5), and converted back to a coefficient afterwards (using
B=sigma\^2). Negative A values are converted to positive A values (using
abs(A)) before mixing, and converted back after mixing (by multiplying
by min(sign(Ai),sign(Aj))). This way, if either particle is repulsive
(if Ai\<0 or Aj\<0), then the default interaction between both particles
will be repulsive.

For the *gauss* style there is no effect due to the Gaussian well beyond
the cutoff; hence reasonable cutoffs need to be specified.

The *gauss/cut* style supports the [pair_modify](pair_modify) shift
option for the energy of the Gauss-potential portion of the pair
interaction.

The [pair_modify](pair_modify) table and tail options are not relevant
for these pair styles.

These pair styles write their information to [binary restart
files](restart), so pair_style and pair_coeff commands do not need to be
specified in an input script that reads a restart file.

These pair styles can only be used via the *pair* keyword of the
[run_style respa](run_style) command. They do not support the *inner*,
*middle*, *outer* keywords.

The *gauss* pair style tallies an \"occupancy\" count of how many
Gaussian-well sites have an atom within the distance at which the force
is a maximum = sqrt(0.5/b). This quantity can be accessed via the
[compute pair](compute_pair) command as a vector of values of length 1.

To print this quantity to the log file (with a descriptive column
heading) the following commands could be included in an input script:

``` LAMMPS
compute gauss all pair gauss
variable occ equal c_gauss[1]
thermo_style custom step temp epair v_occ
```

------------------------------------------------------------------------

## Restrictions

The *gauss* and *gauss/cut* styles are part of the EXTRA-PAIR package.
They are only enabled if LAMMPS is build with that package. See the
[Build package](Build_package) page for more info.

::: versionchanged
28Mar2023
:::

Prior to this version, the *gauss* pair style did not apply
[special_bonds](special_bonds) factors.

## Related commands

[pair_coeff](pair_coeff), [pair_style coul/diel](pair_coul_diel)

## Default

none

::: {#Lenart2}
**(Lenart)** Lenart , Jusufi, and Panagiotopoulos, J Chem Phys, 126,
044509 (2007).
:::

::: {#Jusufi2}
**(Jusufi)** Jusufi, Hynninen, and Panagiotopoulos, J Phys Chem B, 112,
13783 (2008).
:::
