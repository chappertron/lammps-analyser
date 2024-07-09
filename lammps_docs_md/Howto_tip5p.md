# TIP5P water model

The five-point TIP5P rigid water model extends the [three-point TIP3P
model](Howto_tip3p) by adding two additional sites L, usually massless,
where the charge associated with the oxygen atom is placed. These sites
L are located at a fixed distance away from the oxygen atom, forming a
tetrahedral angle that is rotated by 90 degrees from the HOH plane.
Those sites thus somewhat approximate lone pairs of the oxygen and
consequently improve the water structure to become even more
\"tetrahedral\" in comparison to the [four-point TIP4P
model](Howto_tip4p).

A suitable pair style with cutoff Coulomb would be:

-   [pair_style lj/cut/coul/cut](pair_lj_cut_coul)

or these commands for a long-range model:

-   [pair_style lj/cut/coul/long](pair_lj_cut_coul)
-   [pair_style lj/cut/coul/long/soft](pair_fep_soft)
-   [kspace_style pppm](kspace_style)
-   [kspace_style pppm/disp](kspace_style)

A TIP5P model *must* be run using a [rigid fix](fix_rigid) since there
is no other option to keep this kind of structure rigid in LAMMPS. In
order to avoid that LAMMPS produces an error due to the massless L
sites, those need to be assigned a tiny non-zero mass.

The table below lists the force field parameters (in real
[units](units)) to for a the TIP5P model with a cutoff
[(Mahoney)](Mahoney) and the TIP5P-E model [(Rick)](Rick) for use with a
long-range Coulombic solver (e.g. Ewald or PPPM in LAMMPS).

> 
>   Parameter                                         TIP5P              TIP5P-E
>   ------------------------------------------------- ------------------ ------------------
>   O mass (amu)                                      15.9994            15.9994
>   H mass (amu)                                      1.008              1.008
>   O charge ($e$)                                    0.0                0.0
>   L charge ($e$)                                    -0.241             -0.241
>   H charge ($e$)                                    0.241              0.241
>   LJ $\epsilon$ of OO (kcal/mole)                   0.1600             0.1780
>   LJ $\sigma$ of OO ($\AA$)                         3.1200             3.0970
>   LJ $\epsilon$ of HH, LL, OH, OL, HL (kcal/mole)   0.0                0.0
>   LJ $\sigma$ of HH, LL, OH, OL, HL ($\AA$)         1.0                1.0
>   $r_0$ of OH bond ($\AA$)                          0.9572             0.9572
>   $\theta_0$ of HOH angle                           104.52$^{\circ}$   104.52$^{\circ}$
>   OL distance ($\AA$)                               0.70               0.70
>   $\theta_0$ of LOL angle                           109.47$^{\circ}$   109.47$^{\circ}$

Below is the code for a LAMMPS input file for setting up a simulation of
TIP5P water with a molecule file. Because of using [fix
rigid/nvt/small](fix_rigid) no bonds need to be defined and thus no
extra storage needs to be reserved for them, but we need to switch to
atom style full or use [fix property/atom mol](fix_property_atom) so
that fix rigid/nvt/small can identify rigid bodies by their molecule ID:

``` LAMMPS
units real
atom_style charge
atom_modify map array
region box block -5 5 -5 5 -5 5
create_box 3 box

mass 1 15.9994
mass 2 1.008
mass 3 1.0e-100

pair_style lj/cut/coul/cut 8.0
pair_coeff 1 1 0.160  3.12
pair_coeff 2 2 0.0    1.0
pair_coeff 3 3 0.0    1.0

fix mol all property/atom mol
molecule water tip5p.mol
create_atoms 0 random 33 34564 NULL mol water 25367 overlap 1.33

timestep 0.5
fix integrate all rigid/nvt/small molecule temp 300.0 300.0 100.0
reset_timestep 0
velocity all create 300.0 5463576

thermo_style custom step temp press etotal density pe ke
thermo 1000
run 20000
write_data tip5p.data nocoeff




# Water molecule. Explicit TIP5P geometry for use with fix rigid

5 atoms

Coords

1    0.00000  -0.06556   0.00000
2    0.75695   0.52032   0.00000
3   -0.75695   0.52032   0.00000
4    0.00000  -0.46971   0.57154
5    0.00000  -0.46971  -0.57154

Types

1        1   # O
2        2   # H
3        2   # H
4        3   # L
5        3   # L

Charges

1        0.000
2        0.241
3        0.241
4       -0.241
5       -0.241
```

Wikipedia also has a nice article on [water
models](https://en.wikipedia.org/wiki/Water_model)\_.

------------------------------------------------------------------------

::: {#Mahoney}
**(Mahoney)** Mahoney, Jorgensen, J Chem Phys 112, 8910 (2000)
:::

::: {#Rick}
**(Rick)** Rick, J Chem Phys 120, 6085 (2004)
:::
