# bond_style hybrid command

## Syntax

``` LAMMPS
bond_style hybrid style1 style2 ...
```

-   style1,style2 = list of one or more bond styles

## Examples

> bond_style hybrid harmonic fene bond_coeff 1 harmonic 80.0 1.2
> bond_coeff 2\* fene 30.0 1.5 1.0 1.0

## Description

The *hybrid* style enables the use of multiple bond styles in one
simulation. A bond style is assigned to each bond type. For example,
bonds in a polymer flow (of bond type 1) could be computed with a *fene*
potential and bonds in the wall boundary (of bond type 2) could be
computed with a *harmonic* potential. The assignment of bond type to
style is made via the [bond_coeff](bond_coeff) command or in the data
file.

In the bond_coeff commands, the name of a bond style must be added after
the bond type, with the remaining coefficients being those appropriate
to that style. In the example above, the 2 bond_coeff commands set bonds
of bond type 1 to be computed with a *harmonic* potential with
coefficients 80.0, 1.2 for $K$, $r_0$. All other bond types (2-N) are
computed with a *fene* potential with coefficients 30.0, 1.5, 1.0, 1.0
for $K$, $R_0$, $\epsilon$, $\sigma$.

If bond coefficients are specified in the data file read via the
[read_data](read_data) command, then the same rule applies. E.g.
\"harmonic\" or \"fene\" must be added after the bond type, for each
line in the \"Bond Coeffs\" section, e.g.

    Bond Coeffs

    1 harmonic 80.0 1.2
    2 fene 30.0 1.5 1.0 1.0
    ...

A bond style of *none* with no additional coefficients can be used in
place of a bond style, either in a input script bond_coeff command or in
the data file, if you desire to turn off interactions for specific bond
types.

------------------------------------------------------------------------

## Restrictions

This bond style can only be used if LAMMPS was built with the MOLECULE
package. See the [Build package](Build_package) page for more info.

Unlike other bond styles, the hybrid bond style does not store bond
coefficient info for individual sub-styles in a [binary restart
files](restart). Thus when restarting a simulation from a restart file,
you need to re-specify bond_coeff commands.

## Related commands

[bond_coeff](bond_coeff), [delete_bonds](delete_bonds)

## Default

none
