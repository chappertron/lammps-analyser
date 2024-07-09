# pair_style spin/dmi command

## Syntax

``` LAMMPS
pair_style spin/dmi cutoff
```

-   cutoff = global cutoff pair (distance in metal units)

## Examples

``` LAMMPS
pair_style spin/dmi 4.0
pair_coeff * * dmi 2.6 0.001 1.0 0.0 0.0
pair_coeff 1 2 dmi 4.0 0.00109 0.0 0.0 1.0
```

## Description

Style *spin/dmi* computes the Dzyaloshinskii-Moriya (DM) interaction
between pairs of magnetic spins. According to the expression reported in
[(Rohart)](Rohart), one has the following DM energy:

$$\mathbf{H}_{dm} = \sum_{{ i,j}=1,i\neq j}^{N}
\left( \vec{e}_{ij} \times \vec{D} \right)
\cdot\left(\vec{s}_{i}\times \vec{s}_{j}\right),$$

where $\vec{s}_i$ and $\vec{s}_j$ are two neighboring magnetic spins of
two particles, $\vec{e}_ij = \frac{r_i - r_j}{\left| r_i - r_j \right|}$
is the unit vector between sites *i* and *j*, and $\vec{D}$ is the DM
vector defining the intensity (in eV) and the direction of the
interaction.

In [(Rohart)](Rohart), $\vec{D}$ is defined as the direction normal to
the film oriented from the high spin-orbit layer to the magnetic
ultra-thin film.

The application of a spin-lattice Poisson bracket to this energy (as
described in [(Tranchida)](Tranchida5)) allows to derive a magnetic
torque omega, and a mechanical force F (for spin-lattice calculations
only) for each magnetic particle i:

$$\vec{\omega}_i = -\frac{1}{\hbar} \sum_{j}^{Neighb} \vec{s}_{j}\times \left(\vec{e}_{ij}\times \vec{D} \right)
~~{\rm and}~~
\vec{F}_i = -\sum_{j}^{Neighb} \frac{1}{r_{ij}} \vec{D} \times \left( \vec{s}_{i}\times \vec{s}_{j} \right)$$

More details about the derivation of these torques/forces are reported
in [(Tranchida)](Tranchida5).

For the *spin/dmi* pair style, the following coefficients must be
defined for each pair of atoms types via the [pair_coeff](pair_coeff)
command as in the examples above, or in the data file or restart files
read by the [read_data](read_data) or [read_restart](read_restart)
commands, and set in the following order:

-   rc (distance units)
-   \|D\| (energy units)
-   Dx, Dy, Dz (direction of D)

Note that rc is the radius cutoff of the considered DM interaction,
\|D\| is the norm of the DM vector (in eV), and Dx, Dy and Dz define its
direction.

None of those coefficients is optional. If not specified, the *spin/dmi*
pair style cannot be used.

------------------------------------------------------------------------

## Restrictions

All the *pair/spin* styles are part of the SPIN package. These styles
are only enabled if LAMMPS was built with this package, and if the
atom_style \"spin\" was declared. See the [Build package](Build_package)
page for more info.

## Related commands

[atom_style spin](atom_style), [pair_coeff](pair_coeff),
[pair_eam](pair_eam),

## Default

none

------------------------------------------------------------------------

:::: {#Rohart}
::: {#Tranchida5}
**(Rohart)** Rohart and Thiaville, Physical Review B, 88(18), 184422.
(2013).
:::
::::

**(Tranchida)** Tranchida, Plimpton, Thibaudeau and Thompson, Journal of
Computational Physics, 372, 406-425, (2018).
