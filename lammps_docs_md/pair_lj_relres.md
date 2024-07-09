# pair_style lj/relres command

Accelerator Variants: *lj/relres/omp*

## Syntax

``` LAMMPS
pair_style lj/relres Rsi Rso Rci Rco
```

-   Rsi = inner switching cutoff between the fine-grained and
    coarse-grained potentials (distance units)
-   Rso = outer switching cutoff between the fine-grained and
    coarse-grained potentials (distance units)
-   Rci = inner cutoff beyond which the force smoothing for all
    interactions is applied (distance units)
-   Rco = outer cutoff for all interactions (distance units)

## Examples

``` LAMMPS
pair_style lj/relres 4.0 5.0 8.0 10.0
pair_coeff 1 1 0.5 1.0 1.5 1.1
pair_coeff 2 2 0.5 1.0 0.0 0.0 3.0 3.5 6.0 7.0
```

## Description

Pair style *lj/relres* computes a LJ interaction using the Relative
Resolution (RelRes) framework which applies a fine-grained (FG)
potential between near neighbors and a coarse-grained (CG) potential
between far neighbors [(Chaimovich1)](Chaimovich1). This approach can
improve the computational efficiency by almost an order of magnitude,
while maintaining the correct static and dynamic behavior of a reference
system [(Chaimovich2)](Chaimovich2).

$$\begin{aligned}
E = \left\{\begin{array}{lr}
     4 \epsilon^{\scriptscriptstyle FG} \left[ \left(\frac{\sigma^{FG}}{r}\right)^{12} - \left(\frac{\sigma^{FG}}{r}\right)^6 \right]-\Gamma_{si}, & \quad\mathrm{if}\quad  r< r_{si}, \\
     \sum_{m=0}^{4} \gamma_{sm}\left(r-r_{si}\right)^m-\Gamma_{so} ,   & \quad\mathrm{if}\quad  r_{si}\leq r< r_{so}, \\
     4 \epsilon^{\scriptscriptstyle CG} \left[ \left(\frac{\sigma^{CG}}{r}\right)^{12} -     \left(\frac{\sigma^{CG}}{r}\right)^6 \right]-\Gamma_c, &  \quad\mathrm{if}\quad  r_{so}\leq r<r_{ci}, \\
     \sum_{m=0}^{4} \gamma_{cm}\left(r-r_{ci}\right)^m -\Gamma_c, & \quad\mathrm{if}\quad  r_{ci}\leq r< r_{co}, \\
     0, & \quad\mathrm{if}\quad  r\geq r_{co}.\end{array}\right.
\end{aligned}$$

The FG parameters of the LJ potential ($\epsilon^{FG}$ and
$\sigma^{FG}$) are applied up to the inner switching cutoff, $r_{si}$,
while the CG parameters of the LJ potential ($\epsilon^{CG}$ and
$\sigma^{CG}$) are applied beyond the outer switching cutoff, $r_{so}$.
Between $r_{si}$ and $r_{so}$ a polynomial smoothing function is applied
so that the force and its derivative are continuous between the FG and
CG potentials. An analogous smoothing function is applied between the
inner and outer cutoffs ($r_{ci}$ and $r_{co}$). The offsets
$\Gamma_{si}$, $\Gamma_{so}$ and $\Gamma_{c}$ ensure the continuity of
the energy over the entire domain. The corresponding polynomial
coefficients $\gamma_{sm}$ and $\gamma_{cm}$, as well as the offsets are
automatically computed by LAMMPS.

:::: note
::: title
Note
:::

Energy and force resulting from this methodology can be plotted via the
[pair_write](pair_write) command.
::::

The following coefficients must be defined for each pair of atom types
via the [pair_coeff](pair_coeff) command as in the examples above, or in
the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands, or by mixing as will be described
below:

-   $\epsilon^{FG}$ (energy units)
-   $\sigma^{FG}$ (distance units)
-   $\epsilon^{CG}$ (energy units)
-   $\sigma^{CG}$ (distance units)

Additional parameters can be defined to specify different $r_{si}$,
$r_{so}$, $r_{ci}$, $r_{co}$ for a particular set of atom types:

-   $r_{si}$ (distance units)
-   $r_{so}$ (distance units)
-   $r_{ci}$ (distance units)
-   $r_{co}$ (distance units)

These parameters are optional, and they are used to override the global
cutoffs as defined in the pair_style command. If not specified, the
global values for $r_{si}$, $r_{so}$, $r_{ci}$, and $r_{co}$ are used.
If this override option is employed, all four arguments must be
specified.

------------------------------------------------------------------------

Here are some guidelines for using the pair_style *lj/relres* command.

In general, RelRes focuses on the speedup of pairwise interactions
between all LJ sites. Importantly, it works with any settings and flags
(e.g., [special_bonds](special_bonds) settings and [newton](newton)
flags) that can be used in a molecular simulation with the conventional
LJ potential. In particular, all intramolecular topology with its
energetics (i.e., bonds, angles, etc.) remains unaltered.

At the most basic level in the RelRes framework, all sites are mapped
into clusters. Each cluster is just a collection of sites bonded
together (the bonds themselves are not part of the cluster). In general,
a molecule may be comprised of several clusters, and preferably, no two
sites in a cluster are separated by more than two bonds. There are two
categories of sites in RelRes: \"hybrid\" sites embody both FG and CG
models, while \"ordinary\" sites embody just FG characteristics with no
CG features. A given cluster has a single hybrid site (typically its
central site) and several ordinary sites (typically its peripheral
sites). Notice that while clusters are necessary for the RelRes
parameterization (discussed below), they are not actually defined in
LAMMPS. Besides, the total number of sites in the cluster are called the
\"mapping ratio\", and this substantially impacts the computational
efficiency of RelRes: For a mapping ratio of 3, the efficiency factor is
around 4, and for a mapping ratio of 5, the efficiency factor is around
5 [(Chaimovich2)](Chaimovich2).

The flexibility of LAMMPS allows placing any values for the LJ
parameters in the input script. However, here are the optimal
recommendations for the RelRes parameters, which yield the correct
structural and thermal behavior in a system of interest
[(Chaimovich1)](Chaimovich1). One must first assign a complete set of
parameters for the FG interactions that are applicable to all atom
types. Regarding the parameters for the CG interactions, the rules rely
on the site category (if it is a hybrid or an ordinary site). For atom
types of ordinary sites, $\epsilon^{CG}$ must be set to 0 (zero) while
the specific value of $\sigma^{CG}$ is irrelevant. For atom types of
hybrid sites, the CG parameters should be generally calculated using the
following equations:

$$\sigma_I^{CG}=\frac{\left((\sum_{\alpha\in A}\sqrt{\epsilon_\alpha^{FG}\left(\sigma_\alpha^{FG}\right)^{12}}\right)^{1/2}}{\left((\sum_{\alpha\in A}\sqrt{\epsilon_\alpha^{FG}\left(\sigma_\alpha^{FG}\right)^6}\right)^{1/3}}
\quad\mathrm{and}\quad
\epsilon_I^{CG}=\frac{\left((\sum_{\alpha\in A}\sqrt{\epsilon_\alpha^{FG}\left(\sigma_\alpha^{FG}\right)^6}\right)^4}{\left((\sum_{\alpha\in A}\sqrt{\epsilon_\alpha^{FG}\left(\sigma_\alpha^{FG}\right)^{12}}\right)^2}$$

where $I$ is an atom type of a hybrid site of a particular cluster $A$,
and corresponding with this cluster, the summation proceeds over all of
its sites $\alpha$. These equations are derived from the monopole term
in the underlying Taylor series, and they are indeed relevant only if
geometric mixing is applicable for the FG model; if this is not the
case, Ref. [(Chaimovich2)](Chaimovich2) discusses the alternative
formula, and in such a situation, the pair_coeff command should be
explicitly used for all combinations of atom types $I\;!=J$.

The switching distance (the midpoint between inner and outer switching
cutoffs) is another crucial factor in RelRes: decreasing it improves the
computational efficiency, yet if it is too small, the molecular
simulations may not capture the system behavior correctly. As a rule of
thumb, the switching distance should be approximately
$\,\sim\! 1.5\sigma$ [(Chaimovich1)](Chaimovich1); recommendations can
be found in Ref. [(Chaimovich2)](Chaimovich2). Regarding the switching
smoothing zone, $\,\sim\!0.1\sigma$ is recommended; if desired,
smoothing can be eliminated by setting the inner switching cutoff,
$r_{si}$, equal to the outer switching cutoff, $r_{so}$ (the same is
true for the other cutoffs $r_{ci}$ and $r_{co}$).

------------------------------------------------------------------------

As an example, imagine that in your system, a molecule is comprised just
of one cluster such that one atom type (#1) is associated with its
hybrid site, and another atom type (#2) is associated with its ordinary
sites (in total, there are 2 atom types). If geometric mixing is
applicable, the following commands should be used:

``` LAMMPS
pair_style lj/relres Rsi Rso Rci Rco
pair_coeff 1 1 epsilon_FG1 sigma_FG1 epsilon_CG1 sigma_CG1
pair_coeff 2 2 epsilon_FG2 sigma_FG2 0.0         0.0
pair_modify shift yes
```

In a more complex situation, there may be two distinct clusters in a
system (these two clusters may be on same molecule or on different
molecules), each with its own switching cutoffs. If there are still two
atom types in each cluster as in the earlier example, the commands
should be:

``` LAMMPS
pair_style lj/relres Rsi Rso Rci Rco
pair_coeff 1 1 epsilon_FG1 sigma_FG1 epsilon_CG1 sigma_CG1 Rsi1 Rso1 Rci Rco
pair_coeff 2 2 epsilon_FG2 sigma_FG2 0.0         0.0       Rsi1 Rso1 Rci Rco
pair_coeff 3 3 epsilon_FG3 sigma_FG3 epsilon_CG3 sigma_CG3
pair_coeff 4 4 epsilon_FG4 sigma_FG4 0.0         0.0
pair_modify shift yes
```

In this example, the switching cutoffs for the first cluster (atom types
1 and 2) is defined explicitly in the pair_coeff command which overrides
the global values, while the second cluster (atom types 3 and 4) uses
the global definition from the pair_style command. The emphasis here is
that the atom types that belong to a specific cluster should have the
same switching/cutoff arguments.

In the case that geometric mixing is not applicable, for simulating the
system from the previous example, we recommend using the following
commands:

``` LAMMPS
pair_style lj/relres Rsi Rso Rci Rco
pair_coeff 1 1 epsilon_FG1  sigma_FG1  epsilon_CG1  sigma_CG1  Rsi1  Rso1  Rci Rco
pair_coeff 1 2 epsilon_FG12 sigma_FG12 0.0          0.0        Rsi1  Rso1  Rci Rco
pair_coeff 1 3 epsilon_FG13 sigma_FG13 epsilon_CG13 sigma_CG13 Rsi13 Rso13 Rci Rco
pair_coeff 1 4 epsilon_FG14 sigma_FG14 0.0          0.0        Rsi13 Rso13 Rci Rco
pair_coeff 2 2 epsilon_FG2  sigma_FG2  0.0          0.0        Rsi1  Rso1  Rci Rco
pair_coeff 2 3 epsilon_FG23 sigma_FG23 0.0          0.0        Rsi13 Rso13 Rci Rco
pair_coeff 2 4 epsilon_FG24 sigma_FG24 0.0          0.0        Rsi13 Rso13 Rci Rco
pair_coeff 3 3 epsilon_FG3  sigma_FG3  epsilon_CG3  sigma_CG3
pair_coeff 3 4 epsilon_FG34 sigma_FG34 0.0          0.0
pair_coeff 4 4 epsilon_FG4  sigma_FG4  0.0          0.0
pair_modify shift yes
```

Notice that the CG parameters are mixed only for interactions between
atom types associated with hybrid sites, and that the cutoffs are mixed
on the cluster basis.

More examples can be found in the *examples/relres* folder.

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

For atom type pairs $I,\:J$ with $I\;!=J$, the $\epsilon^{FG}$,
$\sigma^{FG}$, $\epsilon^{CG}$, $\sigma^{CG}$, $r_{si}$, $r_{so}$,
$r_{ci}$, and $r_{co}$ parameters for this pair style can be mixed, if
not defined explicitly. All parameters are mixed according to the
pair_modify mix option. The default mix value is *geometric*, and it is
recommended to use with this *lj/relres* style. See the \"pair_modify\"
command for details.

This pair style supports the [pair_modify](pair_modify) shift option for
the energy of the pair interaction. It is recommended to set this option
to *yes*. Otherwise, the offset $\Gamma_{c}$ is set to zero. Constants
$\Gamma_{si}$ and $\Gamma_{so}$ are not impacted by this option.

The [pair_modify](pair_modify) table option is not relevant for this
pair style.

This pair style does not support the [pair_modify](pair_modify) tail
option for adding long-range tail corrections to energy and pressure,
since the energy of the pair interaction is smoothed to 0.0 at the
cutoff.

This pair style writes its information to [binary restart
files](restart), so pair_style and pair_coeff commands do not need to be
specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
[run_style respa](run_style) command. It does not support the *inner*,
*middle*, *outer* keywords.

------------------------------------------------------------------------

## Restrictions

This pair style is part of the EXTRA-PAIR package. It is only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

## Related commands

[pair_coeff](pair_coeff)

## Default

none

------------------------------------------------------------------------

::: {#Chaimovich1}
**(Chaimovich1)** A. Chaimovich, C. Peter and K. Kremer, J. Chem. Phys.
143, 243107 (2015).
:::

::: {#Chaimovich2}
**(Chaimovich2)** M. Chaimovich and A. Chaimovich, J. Chem. Theory
Comput. 17, 1045-1059 (2021).
:::
