# fix qeq/reaxff command

Accelerator Variants: *qeq/reaxff/kk*, *qeq/reaxff/omp*

## Syntax

    fix ID group-ID qeq/reaxff Nevery cutlo cuthi tolerance params args

-   ID, group-ID are documented in [fix](fix) command

-   qeq/reaxff = style name of this fix command

-   Nevery = perform QEq every this many steps

-   cutlo,cuthi = lo and hi cutoff for Taper radius

-   tolerance = precision to which charges will be equilibrated

-   params = reaxff or a filename

-   one or more keywords or keyword/value pairs may be appended

        keyword = *dual* or *maxiter* or *nowarn*
          *dual* = process S and T matrix in parallel (only for qeq/reaxff/omp)
          *maxiter* N = limit the number of iterations to *N*
          *nowarn* = do not print a warning message if the maximum number of iterations was reached

## Examples

``` LAMMPS
fix 1 all qeq/reaxff 1 0.0 10.0 1.0e-6 reaxff
fix 1 all qeq/reaxff 1 0.0 10.0 1.0e-6 param.qeq maxiter 500
```

## Description

Perform the charge equilibration (QEq) method as described in [(Rappe
and Goddard)](Rappe2) and formulated in [(Nakano)](Nakano2). It is
typically used in conjunction with the ReaxFF force field model as
implemented in the [pair_style reaxff](pair_reaxff) command, but it can
be used with any potential in LAMMPS, so long as it defines and uses
charges on each atom. The [fix qeq/comb](fix_qeq_comb) command should be
used to perform charge equilibration with the [COMB
potential](pair_comb). For more technical details about the charge
equilibration performed by fix qeq/reaxff, see the
[(Aktulga)](qeq-Aktulga) paper.

The QEq method minimizes the electrostatic energy of the system by
adjusting the partial charge on individual atoms based on interactions
with their neighbors. It requires some parameters for each atom type. If
the *params* setting above is the word \"reaxff\", then these are
extracted from the [pair_style reaxff](pair_reaxff) command and the
ReaxFF force field file it reads in. If a file name is specified for
*params*, then the parameters are taken from the specified file and the
file must contain one line for each atom type. The latter form must be
used when performing QeQ with a non-ReaxFF potential. Each line should
be formatted as follows:

    itype chi eta gamma

where *itype* is the atom type from 1 to Ntypes, *chi* denotes the
electronegativity in eV, *eta* denotes the self-Coulomb potential in eV,
and *gamma* denotes the valence orbital exponent. Note that these 3
quantities are also in the ReaxFF potential file, except that eta is
defined here as twice the eta value in the ReaxFF file. Note that unlike
the rest of LAMMPS, the units of this fix are hard-coded to be A, eV,
and electronic charge.

The optional *dual* keyword allows to perform the optimization of the S
and T matrices in parallel. This is only supported for the
*qeq/reaxff/omp* style. Otherwise they are processed separately. The
*qeq/reaxff/kk* style always solves the S and T matrices in parallel.

The optional *maxiter* keyword allows changing the max number of
iterations in the linear solver. The default value is 200.

The optional *nowarn* keyword silences the warning message printed when
the maximum number of iterations was reached. This can be useful for
comparing serial and parallel results where having the same fixed number
of QEq iterations is desired, which can be achieved by using a very
small tolerance and setting *maxiter* to the desired number of
iterations.

:::: note
::: title
Note
:::

In order to solve the self-consistent equations for electronegativity
equalization, LAMMPS imposes the additional constraint that all the
charges in the fix group must add up to zero. The initial charge
assignments should also satisfy this constraint. LAMMPS will print a
warning if that is not the case.
::::

## Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to [binary restart
files](restart). This fix computes a global scalar (the number of
iterations) for access by various [output commands](Howto_output). No
parameter of this fix can be used with the *start/stop* keywords of the
[run](run) command.

This fix is invoked during [energy minimization](minimize).

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

## Restrictions

This fix is part of the REAXFF package. It is only enabled if LAMMPS was
built with that package. See the [Build package](Build_package) page for
more info.

This fix does not correctly handle interactions involving multiple
periodic images of the same atom. Hence, it should not be used for
periodic cell dimensions less than 10 Angstroms.

This fix may be used in combination with [fix efield](fix_efield) and
will apply the external electric field during charge equilibration, but
there may be only one fix efield instance used and the electric field
vector may only have components in non-periodic directions. Equal-style
variables can be used for electric field vector components without any
further settings. Atom-style variables can be used for spatially-varying
electric field vector components, but the resulting electric potential
must be specified as an atom-style variable using the *potential*
keyword for [fix efield]{.title-ref}.

## Related commands

[pair_style reaxff](pair_reaxff), [fix qeq/shielded](fix_qeq)

## Default

maxiter 200

------------------------------------------------------------------------

::: {#Rappe2}
**(Rappe)** Rappe and Goddard III, Journal of Physical Chemistry, 95,
3358-3363 (1991).
:::

::: {#Nakano2}
**(Nakano)** Nakano, Computer Physics Communications, 104, 59-69 (1997).
:::

::: {#qeq-Aktulga}
**(Aktulga)** Aktulga, Fogarty, Pandit, Grama, Parallel Computing, 38,
245-259 (2012).
:::
