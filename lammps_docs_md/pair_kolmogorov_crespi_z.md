# pair_style kolmogorov/crespi/z command

## Syntax

``` LAMMPS
pair_style [hybrid/overlay ...] kolmogorov/crespi/z cutoff
```

## Examples

``` LAMMPS
pair_style hybrid/overlay kolmogorov/crespi/z 20.0
pair_coeff * * none
pair_coeff 1 2 kolmogorov/crespi/z  CC.KC   C C

pair_style hybrid/overlay rebo kolmogorov/crespi/z 14.0
pair_coeff * * rebo                 CH.rebo    C C
pair_coeff 1 2 kolmogorov/crespi/z  CC.KC      C C
```

## Description

The *kolmogorov/crespi/z* style computes the Kolmogorov-Crespi
interaction potential as described in [(Kolmogorov)](KC05). An important
simplification is made, which is to take all normals along the z-axis.

$$\begin{aligned}
E  = & \frac{1}{2} \sum_i \sum_{j \neq i} V_{ij} \\
V_{ij}  = & e^{-\lambda(r_{ij} -z_0)} \left[ C + f(\rho_{ij}) + f(\rho_{ji}) \right] - A \left( \frac{r_{ij}}{z_0}\right)^{-6} + A \left( \frac{\textrm{cutoff}}{z_0}\right)^{-6} \\
\rho_{ij}^2 = & \rho_{ji}^2  =  x_{ij}^2 + y_{ij}^2 \qquad \qquad (\mathbf{n}_i \equiv \mathbf{\hat{z}}) \\
f(\rho)  = &  e^{-(\rho/\delta)^2} \sum_{n=0}^2 C_{2n} \left( \rho/\delta \right)^{2n}
\end{aligned}$$

It is important to have a sufficiently large cutoff to ensure smooth
forces. Energies are shifted so that they go continuously to zero at the
cutoff assuming that the exponential part of $V_{ij}$ (first term)
decays sufficiently fast. This shift is achieved by the last term in the
equation for $V_{ij}$ above.

This potential is intended for interactions between two layers of
graphene. Therefore, to avoid interaction between layers in
multi-layered materials, each layer should have a separate atom type and
interactions should only be computed between atom types of neighboring
layers.

The parameter file (e.g. CC.KC), is intended for use with metal
[units](units), with energies in meV. An additional parameter, *S*, is
available to facilitate scaling of energies in accordance with
[(vanWijk)](vanWijk).

This potential must be used in combination with hybrid/overlay. Other
interactions can be set to zero using pair_style *none*.

## Restrictions

This fix is part of the INTERLAYER package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

## Related commands

[pair_coeff](pair_coeff), [pair_none](pair_none), [pair_style
hybrid/overlay](pair_hybrid), [pair_style drip](pair_drip), [pair_style
ilp/graphene/hbn](pair_ilp_graphene_hbn). [pair_style
kolmogorov/crespi/full](pair_kolmogorov_crespi_full), [pair_style
lebedeva/z](pair_lebedeva_z)

## Default

none

------------------------------------------------------------------------

::: {#KC05}
**(Kolmogorov)** A. N. Kolmogorov, V. H. Crespi, Phys. Rev. B 71, 235415
(2005)
:::

::: {#vanWijk}
**(vanWijk)** M. M. van Wijk, A. Schuring, M. I. Katsnelson, and A.
Fasolino, Physical Review Letters, 113, 135504 (2014)
:::
