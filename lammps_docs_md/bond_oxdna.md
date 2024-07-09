# bond_style oxdna/fene command

# bond_style oxdna2/fene command

# bond_style oxrna2/fene command

## Syntax

``` LAMMPS
bond_style oxdna/fene

bond_style oxdna2/fene

bond_style oxrna2/fene
```

## Examples

``` LAMMPS
bond_style oxdna/fene
bond_coeff * 2.0 0.25 0.7525

bond_style oxdna2/fene
bond_coeff * 2.0 0.25 0.7564

bond_style oxrna2/fene
bond_coeff * 2.0 0.25 0.76107
```

## Description

The *oxdna/fene*, *oxdna2/fene*, and *oxrna2/fene* bond styles use the
potential

$$E = - \frac{\epsilon}{2} \ln \left[ 1 - \left(\frac{r-r_0}{\Delta}\right)^2\right]$$

to define a modified finite extensible nonlinear elastic (FENE)
potential [(Ouldridge)](Ouldridge0) to model the connectivity of the
phosphate backbone in the oxDNA/oxRNA force field for coarse-grained
modelling of DNA/RNA.

The following coefficients must be defined for the bond type via the
[bond_coeff](bond_coeff) command as given in the above example, or in
the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   $\epsilon$ (energy)
-   $\Delta$ (distance)
-   $r_0$ (distance)

:::: note
::: title
Note
:::

The oxDNA bond style has to be used together with the corresponding
oxDNA pair styles for excluded volume interaction *oxdna/excv* ,
stacking *oxdna/stk* , cross-stacking *oxdna/xstk* and coaxial stacking
interaction *oxdna/coaxstk* as well as hydrogen-bonding interaction
*oxdna/hbond* (see also documentation of [pair_style
oxdna/excv](pair_oxdna)). For the oxDNA2 [(Snodin)](Snodin0) bond style
the analogous pair styles *oxdna2/excv* , *oxdna2/stk* , *oxdna2/xstk* ,
*oxdna2/coaxstk* , *oxdna2/hbond* and an additional Debye-Hueckel pair
style *oxdna2/dh* have to be defined. The same applies to the oxRNA2
[(Sulc1)](Sulc01) styles. The coefficients in the above example have to
be kept fixed and cannot be changed without reparameterizing the entire
model.
::::

:::: note
::: title
Note
:::

This bond style has to be used with the *atom_style hybrid bond
ellipsoid oxdna* (see documentation of [atom_style](atom_style)). The
*atom_style oxdna* stores the 3\'-to-5\' polarity of the nucleotide
strand, which is set through the bond topology in the data file. The
first (second) atom in a bond definition is understood to point towards
the 3\'-end (5\'-end) of the strand.
::::

:::: warning
::: title
Warning
:::

If data files are produced with [write_data](write_data), then the
[newton](newton) command should be set to *newton on* or *newton off
on*. Otherwise the data files will not have the same 3\'-to-5\' polarity
as the initial data file. This limitation does not apply to binary
restart files produced with [write_restart](write_restart).
::::

Example input and data files for DNA and RNA duplexes can be found in
examples/PACKAGES/cgdna/examples/oxDNA/ , /oxDNA2/ and /oxRNA2/. A
simple python setup tool which creates single straight or helical DNA
strands, DNA/RNA duplexes or arrays of DNA/RNA duplexes can be found in
examples/PACKAGES/cgdna/util/.

Please cite [(Henrich)](Henrich0) in any publication that uses this
implementation. An updated documentation that contains general
information on the model, its implementation and performance as well as
the structure of the data and input file can be found
[here](PDF/CG-DNA.pdf)\_.

Please cite also the relevant oxDNA/oxRNA publications. These are
[(Ouldridge)](Ouldridge0) and [(Ouldridge-DPhil)](Ouldridge-DPhil0) for
oxDNA, [(Snodin)](Snodin0) for oxDNA2, [(Sulc1)](Sulc01) for oxRNA2 and
for sequence-specific hydrogen-bonding and stacking interactions
[(Sulc2)](Sulc02).

------------------------------------------------------------------------

## Restrictions

This bond style can only be used if LAMMPS was built with the CG-DNA
package and the MOLECULE and ASPHERE package. See the [Build
package](Build_package) page for more info.

## Related commands

[pair_style oxdna/excv](pair_oxdna), [pair_style
oxdna2/excv](pair_oxdna2), [pair_style oxrna2/excv](pair_oxrna2),
[bond_coeff](bond_coeff), [atom_style oxdna](atom_style), [fix
nve/dotc/langevin](fix_nve_dotc_langevin)

## Default

none

------------------------------------------------------------------------

::: {#Henrich0}
**(Henrich)** O. Henrich, Y. A. Gutierrez-Fosado, T. Curk, T. E.
Ouldridge, Eur. Phys. J. E 41, 57 (2018).
:::

::: {#Ouldridge-DPhil0}
**(Ouldridge-DPhil)** T.E. Ouldridge, Coarse-grained modelling of DNA
and DNA self-assembly, DPhil. University of Oxford (2011).
:::

::: {#Ouldridge0}
**(Ouldridge)** T.E. Ouldridge, A.A. Louis, J.P.K. Doye, J. Chem. Phys.
134, 085101 (2011).
:::

::: {#Snodin0}
**(Snodin)** B.E. Snodin, F. Randisi, M. Mosayebi, et al., J. Chem.
Phys. 142, 234901 (2015).
:::

::: {#Sulc01}
**(Sulc1)** P. Sulc, F. Romano, T. E. Ouldridge, et al., J. Chem. Phys.
140, 235102 (2014).
:::

::: {#Sulc02}
**(Sulc2)** P. Sulc, F. Romano, T.E. Ouldridge, L. Rovigatti, J.P.K.
Doye, A.A. Louis, J. Chem. Phys. 137, 135101 (2012).
:::
