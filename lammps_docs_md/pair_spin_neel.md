# pair_style spin/neel command

## Syntax

``` LAMMPS
pair_style spin/neel cutoff
```

-   cutoff = global cutoff pair (distance in metal units)

## Examples

``` LAMMPS
pair_style spin/neel 4.0
pair_coeff * * neel 4.0 0.0048 0.234 1.168 2.6905 0.705 0.652
pair_coeff 1 2 neel 4.0 0.0048 0.234 1.168 0.0 0.0 1.0
```

## Description

Style *spin/neel* computes the Neel pair anisotropy model between pairs
of magnetic spins:

$$\mathcal{H}_{N\acute{e}el}=-\sum_{{ i,j=1,i\neq j}}^N g_1(r_{ij})\left(({\mathbf{e}}_{ij}\cdot {\mathbf{s}}_{i})({\mathbf{e}}_{ij}
\cdot {\mathbf{s}}_{j})-\frac{{\mathbf{s}}_{i}\cdot{\mathbf{s}}_{j}}{3} \right)
+q_1(r_{ij})\left( ({\mathbf{e}}_{ij}\cdot {\mathbf{s}}_{i})^2 -\frac{{\mathbf{s}}_{i}\cdot{\mathbf{s}}_{j}}{3}\right)
\left( ({\mathbf{e}}_{ij}\cdot {\mathbf{s}}_{i})^2 -\frac{{\mathbf{s}}_{i}\cdot{\mathbf{s}}_{j}}{3} \right)
+ q_2(r_{ij}) \Big( ({\mathbf{e}}_{ij}\cdot {\mathbf{s}}_{i}) ({\mathbf{e}}_{ij}\cdot {\mathbf{s}}_{j})^3 + ({\mathbf{e}}_{ij}\cdot
{\mathbf{s}}_{j}) ({\mathbf{e}}_{ij}\cdot {\mathbf{s}}_{i})^3\Big)$$

where $\mathbf{s}_i$ and $\mathbf{s}_j$ are two neighboring magnetic
spins of two particles,
$r_{ij} = \vert \mathbf{r}_i - \mathbf{r}_j \vert$ is the inter-atomic
distance between the two particles,
$\mathbf{e}_{ij} = \frac{\mathbf{r}_i - \mathbf{r}_j}{\vert \mathbf{r}_i - \mathbf{r}_j\vert}$
is their normalized separation vector and $g_1$, $q_1$ and $q_2$ are
three functions defining the intensity of the dipolar and quadrupolar
contributions, with:

$$\begin{aligned}
g_1(r_{ij}) &= g(r_{ij}) + \frac{12}{35} q(r_{ij}) \\
q_1(r_{ij}) &= \frac{9}{5} q(r_{ij}) \\
q_2(r_{ij}) &= - \frac{2}{5} q(r_{ij})
\end{aligned}$$

With the functions $g(r_{ij})$ and $q(r_{ij})$ defined and fitted
according to the same Bethe-Slater function used to fit the exchange
interaction:

$${J}\left( r_{ij} \right) = 4 a \left( \frac{r_{ij}}{d}  \right)^2 \left( 1 - b \left( \frac{r_{ij}}{d}  \right)^2 \right) e^{-\left( \frac{r_{ij}}{d} \right)^2 }\Theta (R_c - r_{ij})$$

where $a$, $b$ and $d$ are the three constant coefficients defined in
the associated \"pair_coeff\" command.

The coefficients $a$, $b$, and $d$ need to be fitted so that the
function above matches with the values of the magneto-elastic constant
of the materials at stake.

Examples and more explanations about this function and its
parameterization are reported in [(Tranchida)](Tranchida6). More
examples of parameterization will be provided in future work.

From this DM interaction, each spin $i$ will be submitted to a magnetic
torque $\mathbf{\omega}$ and its associated atom to a force $\mathbf{F}$
(for spin-lattice calculations only).

More details about the derivation of these torques/forces are reported
in [(Tranchida)](Tranchida6).

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

::: {#Tranchida6}
**(Tranchida)** Tranchida, Plimpton, Thibaudeau and Thompson, Journal of
Computational Physics, 372, 406-425, (2018).
:::
