# pair_style smatb command

# pair_style smatb/single command

## Syntax

``` LAMMPS
pair_style style args
```

-   style = *smatb* or *smatb/single*
-   args = none

## Examples

``` LAMMPS
pair_style smatb
pair_coeff 1 1 2.88 10.35 4.178 0.210 1.818 4.07293506 4.9883063257983666

pair_style smatb/single
pair_coeff 1 1 2.88 10.35 4.178 0.210 1.818 4.07293506 4.9883063257983666
```

## Description

::: versionadded
4May2022
:::

The *smatb* and *smatb/single* styles compute the Second Moment
Approximation to the Tight Binding [(Cyrot)](Cyrot), [(Gupta)](Gupta),
[(Rosato)](Rosato), given by

$$E_{i}  = \sum_{j,R_{ij}\leq R_{c}} \alpha(R_{ij}) - \sqrt{\sum_{j,R_{ij}\leq R_{c}}\Xi^2(R_{ij})}$$

$R_{ij}$ is the distance between the atom $i$ and $j$. And the two
functions $\alpha\left(r\right)$ and $\Xi\left(r\right)$ are:

$$\begin{aligned}
\alpha\left(r\right)=\left\lbrace\begin{array}{ll}
A e^{-p \left(\frac{r}{R_{0}}-1\right)} & r < R_{sc}\\
a_3\left(r-R_{c}\right)^3+a_4\left(r-R_{c}\right)^4
+a_5\left(r-R_{c}\right)^5& R_{sc} < r < R_{c}
\end{array}
\right.
\end{aligned}$$

$$\begin{aligned}
\Xi\left(r\right)=\left\lbrace\begin{array}{ll}
\xi e^{-q \left(\frac{r}{R_{0}}-1\right)} & r < R_{sc}\\
x_3\left(r-R_{c}\right)^3+x_4\left(r-R_{c}\right)^4
+x_5\left(r-R_{c}\right)^5& R_{sc} < r < R_{c}
\end{array}
\right.
\end{aligned}$$

The polynomial coefficients $a_3$, $a_4$, $a_5$, $x_3$, $x_4$, $x_5$ are
computed by LAMMPS: the two exponential terms and their first and second
derivatives are smoothly reduced to zero, from the inner cutoff $R_{sc}$
to the outer cutoff $R_{c}$.

The *smatb/single* style is an optimization when using only a single
atom type.

## Coefficients

The following coefficients must be defined for each pair of atoms types
via the [pair_coeff](pair_coeff) command as in the examples above, or in
the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands, or by mixing as described below:

-   $R_{0}$ (distance units)
-   $p$ (dimensionless)
-   $q$ (dimensionless)
-   $A$ (energy units)
-   $\xi$ (energy units)
-   $R_{cs}$ (distance units)
-   $R_{c}$ (distance units)

Note that: $R_{0}$ is the nearest neighbor distance, usually coincides
with the diameter of the atoms

See the [run_style](run_style) command for details.

------------------------------------------------------------------------

## Mixing info

For atom type pairs I,J and I != J the coefficients are not
automatically mixed.

------------------------------------------------------------------------

## Restrictions

These pair styles are part of the SMTBQ package and are only enabled if
LAMMPS is built with that package. See the [Build
package](Build_package) page for more info.

These pair styles require the [newton](newton) setting to be \"on\" for
pair interactions.

## Related commands

-   [pair_coeff](pair_coeff)

## Default

none

------------------------------------------------------------------------

::: {#Cyrot}
**(Cyrot)** Cyrot-Lackmann and Ducastelle, Phys Rev. B, 4, 2406-2412
(1971).
:::

::: {#Gupta}
**(Gupta)** Gupta ,Phys Rev. B, 23, 6265-6270 (1981).
:::

::: {#Rosato}
**(Rosato)** Rosato and Guillope and Legrand, Philosophical Magazine A,
59.2, 321-336 (1989).
:::
