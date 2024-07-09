# pair_style lj/mdf command

# pair_style buck/mdf command

# pair_style lennard/mdf command

## Syntax

``` LAMMPS
pair_style style args
```

-   style = *lj/mdf* or *buck/mdf* or *lennard/mdf*

-   args = list of arguments for a particular style

        *lj/mdf* args = cutoff1 cutoff2
          cutoff1 = inner cutoff for the start of the tapering function
          cutoff1 = out cutoff for the end of the tapering function
        *buck/mdf* args = cutoff1 cutoff2
          cutoff1 = inner cutoff for the start of the tapering function
          cutoff1 = out cutoff for the end of the tapering function
        *lennard/mdf* args = cutoff1 cutoff2
          cutoff1 = inner cutoff for the start of the tapering function
          cutoff1 = out cutoff for the end of the tapering function

## Examples

``` LAMMPS
pair_style lj/mdf 2.5 3.0
pair_coeff * * 1.0 1.0
pair_coeff 1 1 1.1 2.8 3.0 3.2

pair_style buck/mdf 2.5 3.0
pair_coeff * * 100.0 1.5 200.0
pair_coeff * * 100.0 1.5 200.0 3.0 3.5

pair_style lennard/mdf 2.5 3.0
pair_coeff * * 1.0 1.0
pair_coeff 1 1 1021760.3664 2120.317338 3.0 3.2
```

## Description

The *lj/mdf*, *buck/mdf* and *lennard/mdf* compute the standard 12-6
Lennard-Jones and Buckingham potential with the addition of a taper
function that ramps the energy and force smoothly to zero between an
inner and outer cutoff.

$$E_{smooth}(r) = E(r)*f(r)$$

The tapering, *f(r)*, is done by using the Mei, Davenport, Fernando
function [(Mei)](Mei).

$$\begin{aligned}
f(r) & = 1.0  \qquad \qquad \mathrm{for} \qquad r < r_m \\
f(r) & = (1 - x)^3*(1+3x+6x^2) \quad \mathrm{for} \qquad r_m < r < r_{cut} \\
f(r) & = 0.0  \qquad \qquad \mathrm{for} \qquad  r >= r_{cut} \\
\end{aligned}$$

where

$$x = \frac{(r-r_m)}{(r_{cut}-r_m)}$$

Here $r_m$ is the inner cutoff radius and $r_{cut}$ is the outer cutoff
radius.

------------------------------------------------------------------------

For the *lj/mdf* pair_style, the potential energy, *E(r)*, is the
standard 12-6 Lennard-Jones written in the epsilon/sigma form:

$$E(r) = 4 \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} -
                         \left(\frac{\sigma}{r}\right)^6 \right]$$

Either the first two or all of the following coefficients must be
defined for each pair of atoms types via the pair_coeff command as in
the examples above, or in the data file read by the
[read_data](read_data). The two cutoffs default to the global values and
$\epsilon$ and $\sigma$ can also be determined by mixing as described
below:

-   $\epsilon$ (energy units)
-   $\sigma$ (distance units)
-   $r_m$ (distance units)
-   $r_{cut}$ (distance units)

------------------------------------------------------------------------

For the *buck/mdf* pair_style, the potential energy, *E(r)*, is the
standard Buckingham potential with three required coefficients. The two
cutoffs can be omitted and default to the corresponding global values:

$$E(r) = A e^{(-r/\rho)} -\frac{C}{r^6}$$

-   *A* (energy units)
-   $\rho$ (distance units)
-   *C* (energy-distance\^6 units)
-   $r_m$ (distance units)
-   $r_{cut}$ (distance units)

------------------------------------------------------------------------

For the *lennard/mdf* pair_style, the potential energy, *E(r)*, is the
standard 12-6 Lennard-Jones written in the A/B form:

$$E(r) = \frac{A}{r^{12}} - \frac{B}{r^{6}}$$

The following coefficients must be defined for each pair of atoms types
via the pair_coeff command as in the examples above, or in the data file
read by the read_data commands, or by mixing as described below. The two
cutoffs default to their global values and must be either both given or
both left out:

-   *A* (energy-distance\^12 units)
-   *B* (energy-distance\^6 units)
-   $r_m$ (distance units)
-   $r_{cut}$ (distance units)

------------------------------------------------------------------------

## Mixing, shift, table, tail correction, restart, rRESPA info

For atom type pairs I,J and I != J, the $\epsilon$ and $\sigma$
coefficients and cutoff distances for the lj/mdf pair style can be
mixed. The default mix value is *geometric*. See the \"pair_modify\"
command for details. The other two pair styles buck/mdf and lennard/mdf
do not support mixing, so all I,J pairs of coefficients must be
specified explicitly.

None of the lj/mdf, buck/mdf, or lennard/mdf pair styles supports the
[pair_modify](pair_modify) shift option or long-range tail corrections
to pressure and energy.

These styles write their information to [binary restart files](restart),
so pair_style and pair_coeff commands do not need to be specified in an
input script that reads a restart file.

These styles can only be used via the *pair* keyword of the [run_style
respa](run_style) command. They do not support the *inner*, *middle*,
*outer* keywords.

------------------------------------------------------------------------

## Restrictions

These pair styles can only be used if LAMMPS was built with the
EXTRA-PAIR package. See the [Build package](Build_package) doc page for
more info.

## Related commands

[pair_coeff](pair_coeff)

## Default

none

------------------------------------------------------------------------

::: {#Mei}
**(Mei)** Mei, Davenport, Fernando, Phys Rev B, 43 4653 (1991)
:::
