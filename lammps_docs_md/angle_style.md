# angle_style command

## Syntax

``` LAMMPS
angle_style style
```

-   style = *none* or *zero* or *hybrid* or *amoeba* or *charmm* or
    *class2* or *class2/p6* or *cosine* or *cosine/buck6d* or
    *cosine/delta* or *cosine/periodic* or *cosine/shift* or
    *cosine/shift/exp* or *cosine/squared* or *cross* or *dipole* or
    *fourier* or *fourier/simple* or *gaussian* or *harmonic* or
    *lepton* or *mm3* or *quartic* or *spica* or *table*

## Examples

``` LAMMPS
angle_style harmonic
angle_style charmm
angle_style hybrid harmonic cosine
```

## Description

Set the formula(s) LAMMPS uses to compute angle interactions between
triplets of atoms, which remain in force for the duration of the
simulation. The list of angle triplets is read in by a
[read_data](read_data) or [read_restart](read_restart) command from a
data or restart file.

Hybrid models where angles are computed using different angle potentials
can be setup using the *hybrid* angle style.

The coefficients associated with a angle style can be specified in a
data or restart file or via the [angle_coeff](angle_coeff) command.

All angle potentials store their coefficient data in binary restart
files which means angle_style and [angle_coeff](angle_coeff) commands do
not need to be re-specified in an input script that restarts a
simulation. See the [read_restart](read_restart) command for details on
how to do this. The one exception is that angle_style *hybrid* only
stores the list of sub-styles in the restart file; angle coefficients
need to be re-specified.

:::: note
::: title
Note
:::

When both an angle and pair style is defined, the
[special_bonds](special_bonds) command often needs to be used to turn
off (or weight) the pairwise interaction that would otherwise exist
between 3 bonded atoms.
::::

In the formulas listed for each angle style, *theta* is the angle
between the 3 atoms in the angle.

------------------------------------------------------------------------

Here is an alphabetic list of angle styles defined in LAMMPS. Click on
the style to display the formula it computes and coefficients specified
by the associated [angle_coeff](angle_coeff) command.

Click on the style to display the formula it computes, any additional
arguments specified in the angle_style command, and coefficients
specified by the associated [angle_coeff](angle_coeff) command.

There are also additional accelerated pair styles included in the LAMMPS
distribution for faster performance on CPUs, GPUs, and KNLs. The
individual style names on the [Commands angle](angle) page are followed
by one or more of (g,i,k,o,t) to indicate which accelerated styles
exist.

-   [none](angle_none) - turn off angle interactions
-   [zero](angle_zero) - topology but no interactions
-   [hybrid](angle_hybrid) - define multiple styles of angle
    interactions
-   [amoeba](angle_amoeba) - AMOEBA angle
-   [charmm](angle_charmm) - CHARMM angle
-   [class2](angle_class2) - COMPASS (class 2) angle
-   [class2/p6](angle_class2) - COMPASS (class 2) angle expanded to 6th
    order
-   [cosine](angle_cosine) - angle with cosine term
-   [cosine/buck6d](angle_cosine_buck6d) - same as cosine with
    Buckingham term between 1-3 atoms
-   [cosine/delta](angle_cosine_delta) - angle with difference of
    cosines
-   [cosine/periodic](angle_cosine_periodic) - DREIDING angle
-   [cosine/shift](angle_cosine_shift) - angle cosine with a shift
-   [cosine/shift/exp](angle_cosine_shift_exp) - cosine with shift and
    exponential term in spring constant
-   [cosine/squared](angle_cosine_squared) - angle with cosine squared
    term
-   [cross](angle_cross) - cross term coupling angle and bond lengths
-   [dipole](angle_dipole) - angle that controls orientation of a point
    dipole
-   [fourier](angle_fourier) - angle with multiple cosine terms
-   [fourier/simple](angle_fourier_simple) - angle with a single cosine
    term
-   [gaussian](angle_gaussian) - multi-centered Gaussian-based angle
    potential
-   [harmonic](angle_harmonic) - harmonic angle
-   [lepton](angle_lepton) - angle potential from evaluating a string
-   [mesocnt](angle_mesocnt) - piecewise harmonic and linear angle for
    bending-buckling of nanotubes
-   [mm3](angle_mm3) - anharmonic angle
-   [quartic](angle_quartic) - angle with cubic and quartic terms
-   [spica](angle_spica) - harmonic angle with repulsive SPICA pair
    style between 1-3 atoms
-   [table](angle_table) - tabulated by angle

------------------------------------------------------------------------

## Restrictions

Angle styles can only be set for atom_styles that allow angles to be
defined.

Most angle styles are part of the MOLECULE package. They are only
enabled if LAMMPS was built with that package. See the [Build
package](Build_package) page for more info. The doc pages for individual
bond potentials tell if it is part of a package.

## Related commands

[angle_coeff](angle_coeff)

## Default

``` LAMMPS
angle_style none
```
