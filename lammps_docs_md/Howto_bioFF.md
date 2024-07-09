# CHARMM, AMBER, COMPASS, and DREIDING force fields

A force field has 2 parts: the formulas that define it and the
coefficients used for a particular system. Here we only discuss formulas
implemented in LAMMPS that correspond to formulas commonly used in the
CHARMM, AMBER, COMPASS, and DREIDING force fields. Setting coefficients
is done either from special sections in an input data file via the
[read_data](read_data) command or in the input script with commands like
[pair_coeff](pair_coeff) or [bond_coeff](bond_coeff) and so on. See the
[Tools](Tools) doc page for additional tools that can use CHARMM, AMBER,
or Materials Studio generated files to assign force field coefficients
and convert their output into LAMMPS input.

See [(MacKerell)](howto-MacKerell) for a description of the CHARMM force
field. See [(Cornell)](howto-Cornell) for a description of the AMBER
force field. See [(Sun)](howto-Sun) for a description of the COMPASS
force field.

The interaction styles listed below compute force field formulas that
are consistent with common options in CHARMM or AMBER. See each
command\'s documentation for the formula it computes.

-   [bond_style](bond_harmonic) harmonic
-   [angle_style](angle_charmm) charmm
-   [dihedral_style](dihedral_charmm) charmmfsh
-   [dihedral_style](dihedral_charmm) charmm
-   [pair_style](pair_charmm) lj/charmmfsw/coul/charmmfsh
-   [pair_style](pair_charmm) lj/charmmfsw/coul/long
-   [pair_style](pair_charmm) lj/charmm/coul/charmm
-   [pair_style](pair_charmm) lj/charmm/coul/charmm/implicit
-   [pair_style](pair_charmm) lj/charmm/coul/long
-   [special_bonds](special_bonds) charmm
-   [special_bonds](special_bonds) amber

:::: note
::: title
Note
:::

For CHARMM, newer *charmmfsw* or *charmmfsh* styles were released in
March 2017. We recommend they be used instead of the older *charmm*
styles. See discussion of the differences on the [pair
charmm](pair_charmm) and [dihedral charmm](dihedral_charmm) doc pages.
::::

COMPASS is a general force field for atomistic simulation of common
organic molecules, inorganic small molecules, and polymers which was
developed using ab initio and empirical parameterization techniques. See
the [Tools](Tools) page for the msi2lmp tool for creating LAMMPS
template input and data files from BIOVIA\'s Materials Studio files.
Please note that the msi2lmp tool is very old and largely unmaintained,
so it does not support all features of Materials Studio provided force
field files, especially additions during the last decade. You should
watch the output carefully and compare results, where possible. See
[(Sun)](howto-Sun) for a description of the COMPASS force field.

These interaction styles listed below compute force field formulas that
are consistent with the COMPASS force field. See each command\'s
documentation for the formula it computes.

-   [bond_style](bond_class2) class2
-   [angle_style](angle_class2) class2
-   [dihedral_style](dihedral_class2) class2
-   [improper_style](improper_class2) class2
-   [pair_style](pair_class2) lj/class2
-   [pair_style](pair_class2) lj/class2/coul/cut
-   [pair_style](pair_class2) lj/class2/coul/long
-   [special_bonds](special_bonds) lj/coul 0 0 1

DREIDING is a generic force field developed by the [Goddard
group](http://www.wag.caltech.edu)\_ at Caltech and is useful for
predicting structures and dynamics of organic, biological and main-group
inorganic molecules. The philosophy in DREIDING is to use general force
constants and geometry parameters based on simple hybridization
considerations, rather than individual force constants and geometric
parameters that depend on the particular combinations of atoms involved
in the bond, angle, or torsion terms. DREIDING has an [explicit hydrogen
bond term](pair_hbond_dreiding) to describe interactions involving a
hydrogen atom on very electronegative atoms (N, O, F).

See [(Mayo)](howto-Mayo) for a description of the DREIDING force field

The interaction styles listed below compute force field formulas that
are consistent with the DREIDING force field. See each command\'s
documentation for the formula it computes.

-   [bond_style](bond_harmonic) harmonic
-   [bond_style](bond_morse) morse
-   [angle_style](angle_cosine_squared) cosine/squared
-   [angle_style](angle_harmonic) harmonic
-   [angle_style](angle_cosine) cosine
-   [angle_style](angle_cosine_periodic) cosine/periodic
-   [dihedral_style](dihedral_charmm) charmm
-   [improper_style](improper_umbrella) umbrella
-   [pair_style](pair_buck) buck
-   [pair_style](pair_buck) buck/coul/cut
-   [pair_style](pair_buck) buck/coul/long
-   [pair_style](pair_lj) lj/cut
-   [pair_style](pair_lj_cut_coul) lj/cut/coul/cut
-   [pair_style](pair_lj_cut_coul) lj/cut/coul/long
-   [pair_style](pair_hbond_dreiding) hbond/dreiding/lj
-   [pair_style](pair_hbond_dreiding) hbond/dreiding/morse
-   [special_bonds](special_bonds) dreiding

------------------------------------------------------------------------

::: {#howto-MacKerell}
**(MacKerell)** MacKerell, Bashford, Bellott, Dunbrack, Evanseck, Field,
Fischer, Gao, Guo, Ha, et al, J Phys Chem, 102, 3586 (1998).
:::

::: {#howto-Cornell}
**(Cornell)** Cornell, Cieplak, Bayly, Gould, Merz, Ferguson,
Spellmeyer, Fox, Caldwell, Kollman, JACS 117, 5179-5197 (1995).
:::

::: {#howto-Sun}
**(Sun)** Sun, J. Phys. Chem. B, 102, 7338-7364 (1998).
:::

::: {#howto-Mayo}
**(Mayo)** Mayo, Olfason, Goddard III, J Phys Chem, 94, 8897-8909
(1990).
:::
