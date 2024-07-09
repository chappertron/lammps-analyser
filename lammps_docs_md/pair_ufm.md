# pair_style ufm command

Accelerator Variants: *ufm/gpu*, *ufm/omp*, *ufm/opt*

## Syntax

``` LAMMPS
pair_style ufm cutoff
```

-   cutoff = global cutoff for *ufm* interactions (distance units)

## Examples

``` LAMMPS
pair_style ufm 4.0
pair_coeff 1 1 100.0 1.0 2.5
pair_coeff * * 100.0 1.0

pair_style ufm 4.0
pair_coeff * * 10.0 1.0
variable prefactor equal ramp(10,100)
fix 1 all adapt 1 pair ufm epsilon * * v_prefactor
```

## Description

Style *ufm* computes pairwise interactions using the Uhlenbeck-Ford
model (UFM) potential [(Paula Leite2016)](PL2) which is given by

$$\begin{aligned}
E & = -\varepsilon\, \ln{\left[1-\exp{\left(-r^{2}/\sigma^{2}\right)}\right]} \qquad  r < r_c \\
\varepsilon & = p\,k_B\,T
\end{aligned}$$

where $r_c$ is the cutoff, $\sigma$ is a distance-scale and $\epsilon$
is an energy-scale, i.e., a product of Boltzmann constant $k_B$,
temperature $T$ and the Uhlenbeck-Ford p-parameter which is responsible
to control the softness of the interactions [(Paula Leite2017)](PL1).
This model is useful as a reference system for fluid-phase free-energy
calculations [(Paula Leite2016)](PL2).

The following coefficients must be defined for each pair of atom types
via the [pair_coeff](pair_coeff) command as in the examples above, or in
the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands, or by mixing as described below:

-   $\epsilon$ (energy units)
-   $\sigma$ (distance units)
-   cutoff (distance units)

The last coefficient is optional. If not specified, the global *ufm*
cutoff is used.

The [fix adapt](fix_adapt) command can be used to vary epsilon and sigma
for this pair style over the course of a simulation, in which case
pair_coeff settings for epsilon and sigma must still be specified, but
will be overridden. For example these commands will vary the prefactor
epsilon for all pairwise interactions from 10.0 at the beginning to
100.0 at the end of a run:

``` LAMMPS
variable prefactor equal ramp(10,100)
fix 1 all adapt 1 pair ufm epsilon * * v_prefactor
```

:::: note
::: title
Note
:::

The thermodynamic integration procedure can be performed with this
potential using [fix adapt](fix_adapt). This command will rescale the
force on each atom by varying a scale variable, which always starts with
value 1.0. The syntax is the same described above, however, changing
epsilon to scale. A detailed explanation of how to use this command and
perform nonequilibrium thermodynamic integration in LAMMPS is given in
the paper by [(Freitas)](Freitas2).
::::

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

For atom type pairs I,J and I != J, the epsilon and sigma coefficients
and cutoff distance for this pair style can be mixed. The default mix
value is *geometric*. See the \"pair_modify\" command for details.

This pair style support the [pair_modify](pair_modify) shift option for
the energy of the pair interaction.

The [pair_modify](pair_modify) table and tail are not relevant for this
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

This pair style is part of the EXTRA-PAIR package. It is only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

## Related commands

[pair_coeff](pair_coeff), [fix adapt](fix_adapt)

## Default

none

::: {#PL1}
**(Paula Leite2017)** Paula Leite, Santos-Florez, and de Koning, Phys
Rev E, 96, 32115 (2017).
:::

::: {#PL2}
**(Paula Leite2016)** Paula Leite , Freitas, Azevedo, and de Koning, J
Chem Phys, 126, 044509 (2016).
:::

::: {#Freitas2}
**(Freitas)** Freitas, Asta, and de Koning, Computational Materials
Science, 112, 333 (2016).
:::
