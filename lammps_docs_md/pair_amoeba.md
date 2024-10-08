# pair_style amoeba command

Accelerator Variants: *amoeba/gpu*

# pair_style hippo command

Accelerator Variants: *hippo/gpu*

## Syntax

``` LAMMPS
pair_style style
```

-   style = *amoeba* or *hippo*

## Examples

``` LAMMPS
pair_style amoeba
pair_coeff * * protein.prm.amoeba protein.key.amoeba
```

``` LAMMPS
pair_style hippo
pair_coeff * * water.prm.hippo water.key.hippo
```

## Additional info

-   [Howto amoeba](Howto_amoeba)
-   examples/amoeba
-   tools/amoeba
-   potentials/\*.amoeba
-   potentials/\*.hippo

## Description

The *amoeba* style computes the AMOEBA polarizable field formulated by
Jay Ponder\'s group at the U Washington at St Louis [(Ren)](amoeba-Ren),
[(Shi)](amoeba-Shi). The *hippo* style computes the HIPPO polarizable
force field, an extension to AMOEBA, formulated by Josh Rackers and
collaborators in the Ponder group [(Rackers)](amoeba-Rackers).

These force fields can be used when polarization effects are desired in
simulations of water, organic molecules, and biomolecules including
proteins, provided that parameterizations (Tinker PRM force field files)
are available for the systems you are interested in. Files in the LAMMPS
potentials directory with a \"amoeba\" or \"hippo\" suffix can be used.
The Tinker distribution and website have additional force field files as
well.

As discussed on the [Howto amoeba](Howto_amoeba) doc page, the
intermolecular (non-bonded) portion of the AMOEBA force field contains
these terms:

$$U_{amoeba} = U_{multipole} + U_{polar} + U_{hal}$$

while the HIPPO force field contains these terms:

$$U_{hippo} = U_{multipole} + U_{polar} + U_{qxfer} + U_{repulsion} + U_{dispersion}$$

Conceptually, these terms compute the following interactions:

-   $U_{hal}$ = buffered 14-7 van der Waals with offsets applied to
    hydrogen atoms
-   $U_{repulsion}$ = Pauli repulsion due to rearrangement of electron
    density
-   $U_{dispersion}$ = dispersion between correlated, instantaneous
    induced dipole moments
-   $U_{multipole}$ = electrostatics between permanent point charges,
    dipoles, and quadrupoles
-   $U_{polar}$ = electronic polarization between induced point dipoles
-   $U_{qxfer}$ = charge transfer effects

Note that the AMOEBA versus HIPPO force fields typically compute the
same term differently using their own formulas. The references on this
doc page give full details for both force fields.

The formulas for the AMOEBA energy terms are:

$$\begin{aligned}
U_{hal} = & \epsilon_{ij} \left( \frac{1.07}{\rho_{ij} + 0.07} \right)^7 \left( \frac{1.12}{\rho_{ij}^7 + 0.12} - 2 \right) \\
U_{multipole} = & \vec{M}_i\boldsymbol{T_{ij}}\vec{M}_j, \quad \mbox{with} \quad
   \vec{M} = \left(q, \vec{\mu}_{perm}, \boldsymbol{\Theta} \right) \\
U_{polar} = & \frac{1}{2}\vec{\mu}_i^{ind} \vec{E}_i^{perm}
\end{aligned}$$

The formulas for the HIPPO energy terms are:

$$\begin{aligned}
U_{multipole} = & Z_i \frac{1}{r_{ij}} Z_j + Z_i T_{ij}^{damp} \vec{M}_j + Z_j T_{ji}^{damp} \vec{M}_i + \vec{M}_i T_{ij}^{damp} \vec{M}_j, \quad \mbox{with} \quad
   \vec{M} = \left(q, \vec{\mu}_{perm}, \boldsymbol{\Theta} \right) \\
U_{polar} = & \frac{1}{2}\vec{\mu}_i^{ind} \vec{E}_i^{perm} \\
U_{qxfer} = & \epsilon_i e^{-\eta_j r_{ij}} + \epsilon_j e^{-\eta_i r_{ij}} \\
U_{repulsion} = & \frac{K_i K_j}{r_{ij}} S^2
   S^2 = \left( \int{\phi_i \phi_j} dv \right)^2 = \vec{M}_i\boldsymbol{T_{ij}^{repulsion}}\vec{M}_j \\
U_{dispersion} = & -\frac{C_6^iC_6^j}{r_{ij}^6} \left( f_{damp}^{dispersion} \right)_{ij}^2
\end{aligned}$$

:::: note
::: title
Note
:::

The AMOEBA and HIPPO force fields compute long-range charge, dipole, and
quadrupole interactions as well as long-range dispersion effects.
However, unlike other models with long-range interactions in LAMMPS,
this does not require use of a KSpace style via the
[kspace_style](kspace_style) command. That is because for AMOEBA and
HIPPO the long-range computations are intertwined with the pairwise
computations. So these pair style include both short-and long-range
computations. This means the energy and virial computed by the pair
style as well as the \"Pair\" timing reported by LAMMPS will include the
long-range calculations.
::::

The implementation of the AMOEBA and HIPPO force fields in LAMMPS was
done using F90 code provided by the Ponder group from their [Tinker MD
code](https://dasher.wustl.edu/tinker/)\_.

The current implementation (July 2022) of AMOEBA in LAMMPS matches the
version discussed in [(Ponder)](amoeba-Ponder), [(Ren)](amoeba-Ren), and
[(Shi)](amoeba-Shi). Likewise the current implementation of HIPPO in
LAMMPS matches the version discussed in [(Rackers)](amoeba-Rackers).

::: versionadded
8Feb2023
:::

Accelerator support via the GPU package is available.

------------------------------------------------------------------------

Only a single pair_coeff command is used with either the *amoeba* and
*hippo* styles which specifies two Tinker files, a PRM and KEY file.

``` LAMMPS
pair_coeff * * ../potentials/protein.prm.amoeba ../potentials/protein.key.amoeba
pair_coeff * * ../potentials/water.prm.hippo ../potentials/water.key.hippo
```

Examples of the PRM files are in the potentials directory with an
\*.amoeba or \*.hippo suffix. The examples/amoeba directory has examples
of both PRM and KEY files.

A Tinker PRM file is composed of sections, each of which has multiple
lines. A Tinker KEY file is composed of lines, each of which has a
keyword followed by zero or more parameters.

The list of PRM sections and KEY keywords which LAMMPS recognizes are
listed on the [Howto amoeba](Howto_amoeba) doc page. If not recognized,
the section or keyword is skipped.

Note that if the KEY file is specified as NULL, then no file is
required; default values for various AMOEBA/HIPPO settings are used. The
[Howto amoeba](Howto_amoeba) doc page also gives the default settings.

------------------------------------------------------------------------

::: versionadded
3Nov2022
:::

The *amoeba* and *hippo* pair styles support extraction of two per-atom
quantities by the [fix pair](fix_pair) command. This allows the
quantities to be output to files by the [dump](dump) or otherwise
processed by other LAMMPS commands.

The names of the two quantities are \"uind\" and \"uinp\" for the
induced dipole moments for each atom. Neither quantity needs to be
triggered by the [fix pair](fix_pair) command in order for these pair
styles to calculate it.

------------------------------------------------------------------------

## Mixing, shift, table, tail correction, restart, rRESPA info

These pair styles do not support the [pair_modify](pair_modify) mix,
shift, table, and tail options.

These pair styles do not write their information to [binary restart
files](restart), since it is stored in potential files. Thus, you need
to re-specify the pair_style and pair_coeff commands in an input script
that reads a restart file.

These pair styles can only be used via the *pair* keyword of the
[run_style respa](run_style) command. They do not support the *inner*,
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

:::: note
::: title
Note
:::

Using the GPU accelerated pair styles \'amoeba/gpu\' or \'hippo/gpu\'
when compiling the GPU package for OpenCL has a few known issues when
running on integrated GPUs and the calculation may crash.

The GPU accelerated pair styles are also not (yet) compatible with
single precision FFTs.
::::

------------------------------------------------------------------------

## Restrictions

These pair styles are part of the AMOEBA package. They are only enabled
if LAMMPS was built with that package. See the [Build
package](Build_package) doc page for more info.

The AMOEBA and HIPPO potential (PRM) and KEY files provided with LAMMPS
in the potentials and examples/amoeba directories are Tinker files
parameterized for Tinker units. Their numeric parameters are converted
by LAMMPS to its real units [units](units). Thus you can only use these
pair styles with real units.

These potentials do not yet calculate per-atom energy or virial
contributions.

As explained on the [AMOEBA and HIPPO howto](Howto_amoeba) page, use of
these pair styles to run a simulation with the AMOEBA or HIPPO force
fields requires several things.

The first is a data file generated by the tools/tinker/tinker2lmp.py
conversion script which uses Tinker file force field file input to
create a data file compatible with LAMMPS.

The second is use of these commands:

-   [atom_style amoeba](atom_style)
-   [fix property/atom](fix_property_atom)
-   [special_bonds one/five](special_bonds)

And third, depending on the model being simulated, these commands for
intramolecular interactions may also be required:

-   [bond_style class2](bond_class2)
-   [angle_style amoeba](angle_amoeba)
-   [dihedral_style fourier](dihedral_fourier)
-   [improper_style amoeba](improper_amoeba)
-   [fix amoeba/pitorsion](fix_amoeba_pitorsion)
-   [fix amoeba/bitorsion](fix_amoeba_bitorsion)

------------------------------------------------------------------------

## Related commands

[atom_style amoeba](atom_style), [bond_style class2](bond_class2),
[angle_style amoeba](angle_amoeba), [dihedral_style
fourier](dihedral_fourier), [improper_style amoeba](improper_amoeba),
[fix amoeba/pitorsion](fix_amoeba_pitorsion), [fix
amoeba/bitorsion](fix_amoeba_bitorsion), [special_bonds
one/five](special_bonds), [fix property/atom](fix_property_atom)

## Default

none

------------------------------------------------------------------------

::: {#amoeba-Ponder}
**(Ponder)** Ponder, Wu, Ren, Pande, Chodera, Schnieders, Haque, Mobley,
Lambrecht, DiStasio Jr, M. Head-Gordon, Clark, Johnson, T. Head-Gordon,
J Phys Chem B, 114, 2549-2564 (2010).
:::

::: {#amoeba-Rackers}
**(Rackers)** Rackers, Silva, Wang, Ponder, J Chem Theory Comput, 17,
7056-7084 (2021).
:::

::: {#amoeba-Ren}
**(Ren)** Ren and Ponder, J Phys Chem B, 107, 5933 (2003).
:::

::: {#amoeba-Shi}
**(Shi)** Shi, Xia, Zhang, Best, Wu, Ponder, Ren, J Chem Theory Comp, 9,
4046, 2013.
:::
