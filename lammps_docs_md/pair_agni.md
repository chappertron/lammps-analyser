# pair_style agni command

Accelerator Variants: *agni/omp*

## Syntax

``` LAMMPS
pair_style agni
```

## Examples

``` LAMMPS
pair_style      agni
pair_coeff      * * Al.agni Al
```

## Description

Style *agni* style computes the many-body vectorial force components for
an atom as

$$\begin{aligned}
F_i^u                  & = \sum_t^{N_t}\alpha_t \cdot \exp\left[-\frac{\left(d_{i,t}^u\right)^2}{2l^2}\right] \\
d_{i,t}^u              & = \left|\left| V_i^u(\eta) - V_t^u(\eta) \right|\right| \\
V_i^u(\eta)            & = \sum_{j \neq i}\frac{r^u_{ij}}{r_{ij}} \cdot e^{-\left(\frac{r_{ij}}{\eta} \right)^2} \cdot f_d\left(r_{ij}\right) \\
f_d\left(r_{ij}\right) & = \frac{1}{2} \left[\cos\left(\frac{\pi r_{ij}}{R_c}\right) + 1 \right]
\end{aligned}$$

$u$ labels the individual components, i.e. $x$, $y$ or $z$, and $V$ is
the corresponding atomic fingerprint. $d$ is the Euclidean distance
between any two atomic fingerprints. A total of $N_t$ reference atomic
environments are considered to construct the force field file.
$\alpha_t$ and $l$ are the weight coefficients and length scale
parameter of the non-linear regression model.

The method implements the recently proposed machine learning access to
atomic forces as discussed extensively in the following publications
-[(Botu1)](Botu2015adaptive) and [(Botu2)](Botu2015learning). The
premise of the method is to map the atomic environment numerically into
a fingerprint, and use machine learning methods to create a mapping to
the vectorial atomic forces.

Only a single pair_coeff command is used with the *agni* style which
specifies an AGNI potential file containing the parameters of the force
field for the needed elements. These are mapped to LAMMPS atom types by
specifying $N$ additional arguments after the filename in the pair_coeff
command, where $N$ is the number of LAMMPS atom types:

-   filename
-   $N$ element names = mapping of AGNI elements to atom types

See the [pair_coeff](pair_coeff) page for alternate ways to specify the
path for the force field file.

An AGNI force field is fully specified by the filename which contains
the parameters of the force field, i.e., the reference training
environments used to construct the machine learning force field. Example
force field and input files are provided in the examples/PACKAGES/agni
directory.

------------------------------------------------------------------------

Styles with a *gpu*, *intel*, *kk*, *omp*, or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the [Accelerator packages](Speed_packages)
page. The accelerated styles take the same arguments and should produce
the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, INTEL, KOKKOS, OPENMP, and
OPT packages, respectively. They are only enabled if LAMMPS was built
with those packages. See the [Build package](Build_package) page for
more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the [-suffix command-line
switch](Run_options) when you invoke LAMMPS, or you can use the
[suffix](suffix) command in your input script.

See the [Accelerator packages](Speed_packages) page for more
instructions on how to use the accelerated styles effectively.

------------------------------------------------------------------------

## Mixing, shift, table, tail correction, restart, rRESPA info

This pair style does not support the [pair_modify](pair_modify) shift,
table, and tail options.

This pair style does not write its information to [binary restart
files](restart), since it is stored in potential files. Thus, you need
to re-specify the pair_style and pair_coeff commands in an input script
that reads a restart file.

This pair style can only be used via the *pair* keyword of the
[run_style respa](run_style) command. It does not support the *inner*,
*middle*, *outer* keywords.

------------------------------------------------------------------------

## Restrictions

Currently, only elemental systems are implemented. Also, the method only
provides access to the forces and not energies or stresses. The lack of
potential energy data makes this pair style incompatible with several of
the [minimizer algorthms](min_style) like *cg* or *sd*. It should work
with damped dynamics based minimizers like *fire* or *quickmin*.
However, one can access the energy via thermodynamic integration of the
forces as discussed in [(Botu3)](Botu2016construct). This pair style is
part of the MISC package. It is only enabled if LAMMPS was built with
that package. See the [Build package](Build_package) page for more info.

The AGNI force field files provided with LAMMPS (see the potentials
directory) are parameterized for metal [units](units). You can use the
AGNI potential with any LAMMPS units, but you would need to create your
own AGNI potential file with coefficients listed in the appropriate
units if your simulation does not use \"metal\" units.

## Related commands

[pair_coeff](pair_coeff)

## Default

none

------------------------------------------------------------------------

::: {#Botu2015adaptive}
**(Botu1)** V. Botu and R. Ramprasad, Int. J. Quant. Chem., 115(16),
1074 (2015).
:::

::: {#Botu2015learning}
**(Botu2)** V. Botu and R. Ramprasad, Phys. Rev. B, 92(9), 094306
(2015).
:::

::: {#Botu2016construct}
**(Botu3)** V. Botu, R. Batra, J. Chapman and R. Ramprasad,
<https://arxiv.org/abs/1610.02098> (2016).
:::
