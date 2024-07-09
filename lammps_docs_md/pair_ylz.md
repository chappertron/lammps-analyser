# pair_style ylz command

## Syntax

``` LAMMPS
pair_style ylz cutoff
```

-   cutoff = global cutoff for interactions (distance units)

## Examples

``` LAMMPS
pair_style   ylz  2.6
pair_coeff   *  *  1.0  1.0  4  3  0.0  2.6
```

## Description

::: versionadded
3Nov2022
:::

The *ylz* (Yuan-Li-Zhang) style computes an anisotropic interaction
between pairs of coarse-grained particles considering the relative
particle orientations. This potential was originally developed as a
particle-based solvent-free model for biological membranes
[(Yuan2010a)](Yuan). Unlike [pair_style gayberne](pair_gayberne), whose
orientation dependence is strictly derived from the closest distance
between two ellipsoidal rigid bodies, the orientation-dependence of this
pair style is mathematically defined such that the particles can
self-assemble into one-particle-thick fluid membranes. The potential of
this pair style is described by:

$$\begin{aligned}
U ( \mathbf{r}_{ij}, \mathbf{n}_i, \mathbf{n}_j ) =\left\{\begin{matrix} {u}_R(r)+\left [ 1-\phi (\mathbf{\hat{r}}_{ij}, \mathbf{n}_i, \mathbf{n}_j ) \right ]\epsilon, ~~ r<{r}_{min} \\ {u}_A(r)\phi (\mathbf{\hat{r}}_{ij}, \mathbf{n}_i, \mathbf{n}_j ),~~  {r}_{min}<r<{r}_{c} \\ \end{matrix}\right.\\\\ \phi (\mathbf{\hat{r}}_{ij}, \mathbf{n}_i, \mathbf{n}_j )=1+\left [  \mu (a(\mathbf{\hat{r}}_{ij}, \mathbf{n}_i, \mathbf{n}_j )-1) \right ] \\\\a(\mathbf{\hat{r}}_{ij}, \mathbf{n}_i, \mathbf{n}_j )=(\mathbf{n}_i\times\mathbf{\hat{r}}_{ij} )\cdot (\mathbf{n}_j\times\mathbf{\hat{r}}_{ij} )+{\beta}(\mathbf{n}_i-\mathbf{n}_j)\cdot \mathbf{\hat{r}}_{ij}-\beta^{2}\\\\  {u}_R(r)=\epsilon \left [ \left ( \frac{{r}_{min}}{r} \right )^{4}-2\left ( \frac{{r}_{min}}{r}\right )^{2} \right ] \\\\ {u}_A(r)=-\epsilon\;cos^{2\zeta }\left [ \frac{\pi}{2}\frac{\left ( {r}-{r}_{min} \right )}{\left ( {r}_{c}-{r}_{min} \right )} \right ]\\
\end{aligned}$$

where $\mathbf{r}_{i}$ and $\mathbf{r}_{j}$ are the center position
vectors of particles i and j, respectively,
$\mathbf{r}_{ij}=\mathbf{r}_{i}-\mathbf{r}_{j}$ is the inter-particle
distance vector, $r=\left|\mathbf{r}_{ij} \right|$ and
${\hat{\mathbf{r}}}_{ij}=\mathbf{r}_{ij}/r$. The unit vectors
$\mathbf{n}_{i}$ and $\mathbf{n}_{j}$ represent the axes of symmetry of
particles i and j, respectively, $u_R$ and $u_A$ are the repulsive and
attractive potentials, $\phi$ is an angular function which depends on
the relative orientation between pair particles, $\mu$ is the parameter
related to the bending rigidity of the membrane, $\beta$ is the
parameter related to the spontaneous curvature, and $\epsilon$ is the
energy unit, respectively. The $\zeta$ controls the slope of the
attractive branch and hence the diffusivity of the particles in the
in-plane direction of the membrane. ${r}_{c}$ is the cutoff radius,
$r_{min}$ is the distance which minimizes the potential energy
$u_{A}(r)$ and $r_{min}=2^{1/6}\sigma$, where $\sigma$ is the length
unit.

This pair style is suited for solvent-free coarse-grained simulations of
biological systems involving lipid bilayer membranes, such as vesicle
shape transformations [(Yuan2010b)](Yuan), nanoparticle endocytosis
[(Huang)](Huang), modeling of red blood cell membranes [(Fu)](Fu),
[(Appshaw)](Appshaw), and modeling of cell elasticity
[(Becton)](Becton).

Use of this pair style requires the NVE, NVT, or NPT fixes with the
*asphere* extension (e.g. [fix nve/asphere](fix_nve_asphere)) in order
to integrate particle rotation. Additionally, [atom_style
ellipsoid](atom_style) should be used since it defines the rotational
state of each particle.

The following coefficients must be defined for each pair of atoms types
via the [pair_coeff](pair_coeff) command as in the examples above, or in
the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands, or by mixing as described below:

-   $\epsilon$ = well depth (energy units)
-   $\sigma$ = minimum effective particle radii (distance units)
-   $\zeta$ = tuning parameter for the slope of the attractive branch
-   $\mu$ = parameter related to bending rigidity
-   $\beta$ = parameter related to the spontaneous curvature
-   cutoff (distance units)

The last coefficient is optional. If not specified, the global cutoff
specified in the pair_style command is used.

------------------------------------------------------------------------

## Mixing, shift, table, tail correction, restart, rRESPA info

For atom type pairs I,J and I != J, the epsilon and sigma coefficients
and cutoff distance for this pair style can be mixed. The default mix
value is *geometric*. See the \"pair_modify\" command for details.

The [pair_modify](pair_modify) table option is not relevant for this
pair style.

This pair style does not support the [pair_modify](pair_modify) tail
option for adding long-range tail corrections to energy and pressure.

This pair style writes its information to [binary restart
files](restart), so pair_style and pair_coeff commands do not need to be
specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
[run_style respa](run_style) command. It does not support the *inner*,
*middle*, *outer* keywords.

------------------------------------------------------------------------

## Restrictions

The *ylz* style is part of the ASPHERE package. It is only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

This pair style requires that atoms store torque and a quaternion to
represent their orientation, as defined by the [atom_style](atom_style).
It also requires they store a per-atom [shape](set). The particles
cannot store a per-particle diameter. To avoid being mistakenly
considered as point particles, the shape parameters ought to be
non-spherical, like \[1 0.99 0.99\]. Unlike the
[resquared](pair_resquared) pair style for which the shape directly
determines the mathematical expressions of the potential, the shape
parameters for this pair style is only involved in the computation of
the moment of inertia and thus only influences the rotational dynamics
of individual particles.

This pair style requires that **all** atoms are ellipsoids as defined by
the [atom_style ellipsoid](atom_style) command.

## Related commands

[pair_coeff](pair_coeff), [fix nve/asphere](fix_nve_asphere), [compute
temp/asphere](compute_temp_asphere), [pair_style
resquared](pair_resquared), [pair_style gayberne](pair_gayberne)

## Default

none

------------------------------------------------------------------------

::: {#Yuan}
**(Yuan2010a)** Yuan, Huang, Li, Lykotrafitis, Zhang, Phys. Rev. E, 82,
011905(2010).
:::

**(Yuan2010b)** Yuan, Huang, Zhang, Soft. Matter, 6, 4571(2010).

::: {#Huang}
**(Huang)** Huang, Zhang, Yuan, Gao, Zhang, Nano Lett. 13, 4546(2013).
:::

::: {#Fu}
**(Fu)** Fu, Peng, Yuan, Kfoury, Young, Comput. Phys. Commun, 210,
193-203(2017).
:::

::: {#Appshaw}
**(Appshaw)** Appshaw, Seddon, Hanna, Soft. Matter,18, 1747(2022).
:::

::: {#Becton}
**(Becton)** Becton, Averett, Wang, Biomech. Model. Mechanobiology, 18,
425-433(2019).
:::
