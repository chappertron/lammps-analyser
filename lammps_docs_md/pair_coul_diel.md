# pair_style coul/diel command

Accelerator Variants: *coul/diel/omp*

## Syntax

``` LAMMPS
pair_style coul/diel cutoff
```

cutoff = global cutoff (distance units)

## Examples

``` LAMMPS
pair_style coul/diel 3.5
pair_coeff 1 4 78. 1.375 0.112
```

## Description

Style *coul/diel* computes a Coulomb correction for implicit solvent ion
interactions in which the dielectric permittivity is distance dependent.
The dielectric permittivity epsilon_D(r) connects to limiting regimes:
One limit is defined by a small dielectric permittivity (close to
vacuum) at or close to contact separation between the ions. At larger
separations the dielectric permittivity reaches a bulk value used in the
regular Coulomb interaction coul/long or coul/cut. The transition is
modeled by a hyperbolic function which is incorporated in the Coulomb
correction term for small ion separations as follows

$$\begin{aligned}
E  = & \frac{Cq_iq_j}{\epsilon r} \left( \frac{\epsilon}{\epsilon_D(r)}-1\right)                       \qquad r < r_c \\
\epsilon_D(r)  = & \frac{5.2+\epsilon}{2} +  \frac{\epsilon-5.2}{2}\tanh\left(\frac{r-r_{me}}{\sigma_e}\right)
\end{aligned}$$

where $r_{me}$ is the inflection point of $\epsilon_D(r)$ and $\sigma_e$
is a slope defining length scale. C is the same Coulomb conversion
factor as in the pair_styles coul/cut, coul/long, and coul/debye. In
this way the Coulomb interaction between ions is corrected at small
distances r. The lower limit of epsilon_D(r-\>0)=5.2 due to dielectric
saturation [(Stiles)](Stiles) while the Coulomb interaction reaches its
bulk limit by setting $\epsilon_D(r \to \infty) = \epsilon$, the bulk
value of the solvent which is 78 for water at 298K.

Examples of the use of this type of Coulomb interaction include implicit
solvent simulations of salt ions [(Lenart)](Lenart1) and of ionic
surfactants [(Jusufi)](Jusufi1). Note that this potential is only
reasonable for implicit solvent simulations and in combination with
coul/cut or coul/long. It is also usually combined with gauss/cut, see
[(Lenart)](Lenart1) or [(Jusufi)](Jusufi1).

The following coefficients must be defined for each pair of atom types
via the [pair_coeff](pair_coeff) command as in the example above, or in
the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   $\epsilon$ (no units)
-   $r_{me}$ (distance units)
-   $\sigma_e$ (distance units)

The global cutoff ($r_c$) specified in the pair_style command is used.

------------------------------------------------------------------------

## Mixing, shift, table, tail correction, restart, rRESPA info

This pair style does not support parameter mixing. Coefficients must be
given explicitly for each type of particle pairs.

This pair style supports the [pair_modify](pair_modify) shift option for
the energy of the Gauss-potential portion of the pair interaction.

The [pair_modify](pair_modify) table option is not relevant for this
pair style.

This pair style does not support the [pair_modify](pair_modify) tail
option for adding long-range tail corrections to energy and pressure.

This pair style can only be used via the *pair* keyword of the
[run_style respa](run_style) command. It does not support the *inner*,
*middle*, *outer* keywords.

## Restrictions

This style is part of the EXTRA-PAIR package. It is only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

## Related commands

[pair_coeff](pair_coeff) [pair_style gauss/cut](pair_gauss)

## Default

none

------------------------------------------------------------------------

::: {#Stiles}
**(Stiles)** Stiles , Hubbard, and Kayser, J Chem Phys, 77, 6189 (1982).
:::

::: {#Lenart1}
**(Lenart)** Lenart , Jusufi, and Panagiotopoulos, J Chem Phys, 126,
044509 (2007).
:::

::: {#Jusufi1}
**(Jusufi)** Jusufi, Hynninen, and Panagiotopoulos, J Phys Chem B, 112,
13783 (2008).
:::
