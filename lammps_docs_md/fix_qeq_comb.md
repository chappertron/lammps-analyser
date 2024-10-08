# fix qeq/comb command

Accelerator Variants: *qeq/comb/omp*

## Syntax

    fix ID group-ID qeq/comb Nevery precision keyword value ...

-   ID, group-ID are documented in [fix](fix) command

-   qeq/comb = style name of this fix command

-   Nevery = perform charge equilibration every this many steps

-   precision = convergence criterion for charge equilibration

-   zero or more keyword/value pairs may be appended

-   keyword = *file*

        *file* value = filename
          filename = name of file to write QEQ equilibration info to

## Examples

``` LAMMPS
fix 1 surface qeq/comb 10 0.0001
```

## Description

Perform charge equilibration (QeQ) in conjunction with the COMB
(Charge-Optimized Many-Body) potential as described in
[(COMB_1)](COMB_1) and [(COMB_2)](COMB_2). It performs the charge
equilibration portion of the calculation using the so-called QEq method,
whereby the charge on each atom is adjusted to minimize the energy of
the system. This fix can only be used with the COMB potential; see the
[fix qeq/reaxff](fix_qeq_reaxff) command for a QeQ calculation that can
be used with any potential.

Only charges on the atoms in the specified group are equilibrated. The
fix relies on the pair style (COMB in this case) to calculate the
per-atom electronegativity (effective force on the charges). An
electronegativity equalization calculation (or QEq) is performed in an
iterative fashion, which in parallel requires communication at each
iteration for processors to exchange charge information about nearby
atoms with each other. See [Rappe_and_Goddard](Rappe_and_Goddard) and
[Rick_and_Stuart](Rick_and_Stuart) for details.

During a run, charge equilibration is performed every *Nevery* time
steps. Charge equilibration is also always enforced on the first step of
each run. The *precision* argument controls the tolerance for the
difference in electronegativity for all atoms during charge
equilibration. *Precision* is a trade-off between the cost of performing
charge equilibration (more iterations) and accuracy.

If the *file* keyword is used, then information about each equilibration
calculation is written to the specified file.

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

## Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to [binary restart
files](restart).

The [fix_modify](fix_modify) *respa* option is supported by this fix.
This allows to set at which level of the [r-RESPA](run_style) integrator
the fix is performing charge equilibration. Default is the outermost
level.

This fix produces a per-atom vector which can be accessed by various
[output commands](Howto_output). The vector stores the gradient of the
charge on each atom. The per-atom values be accessed on any timestep.

No parameter of this fix can be used with the *start/stop* keywords of
the [run](run) command.

This fix can be invoked during [energy minimization](minimize).

## Restrictions

This fix command currently only supports [pair style
\*comb\*](pair_comb).

## Related commands

[pair_style comb](pair_comb)

## Default

No file output is performed.

------------------------------------------------------------------------

::: {#COMB_1}
**(COMB_1)** J. Yu, S. B. Sinnott, S. R. Phillpot, Phys Rev B, 75,
085311 (2007),
:::

::: {#COMB_2}
**(COMB_2)** T.-R. Shan, B. D. Devine, T. W. Kemper, S. B. Sinnott, S.
R. Phillpot, Phys Rev B, 81, 125328 (2010).
:::

::: {#Rappe_and_Goddard}
**(Rappe_and_Goddard)** A. K. Rappe, W. A. Goddard, J Phys Chem 95, 3358
(1991).
:::

::: {#Rick_and_Stuart}
**(Rick_and_Stuart)** S. W. Rick, S. J. Stuart, B. J. Berne, J Chem Phys
101, 16141 (1994).
:::
