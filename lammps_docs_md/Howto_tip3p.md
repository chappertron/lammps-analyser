# TIP3P water model

The TIP3P water model as implemented in CHARMM
[(MacKerell)](howto-tip3p) specifies a 3-site rigid water molecule with
charges and Lennard-Jones parameters assigned to each of the 3 atoms.

A suitable pair style with cutoff Coulomb would be:

-   [pair_style lj/cut/coul/cut](pair_lj_cut_coul)

or these commands for a long-range Coulomb model:

-   [pair_style lj/cut/coul/long](pair_lj_cut_coul)
-   [pair_style lj/cut/coul/long/soft](pair_fep_soft)
-   [kspace_style pppm](kspace_style)
-   [kspace_style pppm/disp](kspace_style)

In LAMMPS the [fix shake or fix rattle](fix_shake) command can be used
to hold the two O-H bonds and the H-O-H angle rigid. A bond style of
[harmonic](bond_harmonic) and an angle style of
[harmonic](angle_harmonic) or [charmm](angle_charmm) should also be
used. In case of rigid bonds also bond style [zero](bond_zero) and angle
style [zero](angle_zero) can be used.

The table below lists the force field parameters (in real
[units](units)) to for the water molecule atoms to run a rigid or
flexible TIP3P-CHARMM model with a cutoff, the original 1983 TIP3P model
[(Jorgensen)](Jorgensen1), or a TIP3P model with parameters optimized
for a long-range Coulomb solver (e.g. Ewald or PPPM in LAMMPS)
[(Price)](Price1). The K values can be used if a flexible TIP3P model
(without fix shake) is desired, for rigid bonds/angles they are ignored.

> 
>   Parameter                          TIP3P-CHARMM       TIP3P (original)   TIP3P (Ewald)
>   ---------------------------------- ------------------ ------------------ ------------------
>   O mass (amu)                       15.9994            15.9994            15.9994
>   H mass (amu)                       1.008              1.008              1.008
>   O charge ($e$)                     -0.834             -0.834             -0.834
>   H charge ($e$)                     0.417              0.417              0.417
>   LJ $\epsilon$ of OO (kcal/mole)    0.1521             0.1521             0.1020
>   LJ $\sigma$ of OO ($\AA$)          3.1507             3.1507             3.188
>   LJ $\epsilon$ of HH (kcal/mole)    0.0460             0.0                0.0
>   LJ $\sigma$ of HH ($\AA$)          0.4                1.0                1.0
>   LJ $\epsilon$ of OH (kcal/mole)    0.0836             0.0                0.0
>   LJ $\sigma$ of OH ($\AA$)          1.7753             1.0                1.0
>   K of OH bond (kcal/mole/$\AA^2$)   450                450                450
>   $r_0$ of OH bond ($\AA$)           0.9572             0.9572             0.9572
>   K of HOH angle (kcal/mole)         55.0               55.0               55.0
>   $\theta_0$ of HOH angle            104.52$^{\circ}$   104.52$^{\circ}$   104.52$^{\circ}$

Below is the code for a LAMMPS input file and a molecule file
(`tip3p.mol`) of TIP3P water for use with the [molecule
command](molecule) demonstrating how to set up a small bulk water system
for TIP3P with rigid bonds.

``` LAMMPS
units real
atom_style full
region box block -5 5 -5 5 -5 5
create_box 2 box bond/types 1 angle/types 1 &
            extra/bond/per/atom 2 extra/angle/per/atom 1 extra/special/per/atom 2

mass 1 15.9994
mass 2 1.008

pair_style lj/cut/coul/cut 8.0
pair_coeff 1 1 0.1521 3.1507
pair_coeff 2 2 0.0    1.0

bond_style zero
bond_coeff 1 0.9574

angle_style zero
angle_coeff 1 104.52

molecule water tip3p.mol
create_atoms 0 random 33 34564 NULL mol water 25367 overlap 1.33

fix rigid all shake 0.001 10 10000 b 1 a 1
minimize 0.0 0.0 1000 10000

reset_timestep 0
timestep 1.0
velocity all create 300.0 5463576
fix integrate all nvt temp 300 300 100.0

thermo_style custom step temp press etotal pe

thermo 1000
run 20000
write_data tip3p.data nocoeff




# Water molecule. TIP3P geometry

3 atoms
2 bonds
1 angles

Coords

1    0.00000  -0.06556   0.00000
2    0.75695   0.52032   0.00000
3   -0.75695   0.52032   0.00000

Types

1        1   # O
2        2   # H
3        2   # H

Charges

1       -0.834
2        0.417
3        0.417

Bonds

1   1      1      2
2   1      1      3

Angles

1   1      2      1      3

Shake Flags

1 1
2 1
3 1

Shake Atoms

1 1 2 3
2 1 2 3
3 1 2 3

Shake Bond Types

1 1 1 1
2 1 1 1
3 1 1 1

Special Bond Counts

1 2 0 0
2 1 1 0
3 1 1 0

Special Bonds

1 2 3
2 1 3
3 1 2
```

Wikipedia also has a nice article on [water
models](https://en.wikipedia.org/wiki/Water_model)\_.

------------------------------------------------------------------------

::: {#howto-tip3p}
**(MacKerell)** MacKerell, Bashford, Bellott, Dunbrack, Evanseck, Field,
Fischer, Gao, Guo, Ha, et al, J Phys Chem, 102, 3586 (1998).
:::

::: {#Jorgensen1}
**(Jorgensen)** Jorgensen, Chandrasekhar, Madura, Impey, Klein, J Chem
Phys, 79, 926 (1983).
:::

::: {#Price1}
**(Price)** Price and Brooks, J Chem Phys, 121, 10096 (2004).
:::
