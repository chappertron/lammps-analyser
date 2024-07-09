# pair_style soft command

Accelerator Variants: *soft/gpu*, *soft/omp*

## Syntax

``` LAMMPS
pair_style soft cutoff
```

-   cutoff = global cutoff for soft interactions (distance units)

## Examples

``` LAMMPS
pair_style soft 1.0
pair_coeff * * 10.0
pair_coeff 1 1 10.0 3.0

pair_style soft 1.0
pair_coeff * * 0.0
variable prefactor equal ramp(0,30)
fix 1 all adapt 1 pair soft a * * v_prefactor
```

## Description

Style *soft* computes pairwise interactions with the formula

$$E = A \left[ 1 + \cos\left(\frac{\pi r}{r_c}\right) \right]
\qquad r < r_c$$

It is useful for pushing apart overlapping atoms, since it does not blow
up as r goes to 0. A is a prefactor that can be made to vary in time
from the start to the end of the run (see discussion below), e.g. to
start with a very soft potential and slowly harden the interactions over
time. Rc is the cutoff. See the [fix nve/limit](fix_nve_limit) command
for another way to push apart overlapping atoms.

The following coefficients must be defined for each pair of atom types
via the [pair_coeff](pair_coeff) command as in the examples above, or in
the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands, or by mixing as described below:

-   A (energy units)
-   cutoff (distance units)

The last coefficient is optional. If not specified, the global soft
cutoff is used.

:::: note
::: title
Note
:::

The syntax for [pair_coeff](pair_coeff) with a single A coeff is
different in the current version of LAMMPS than in older versions which
took two values, Astart and Astop, to ramp between them. This
functionality is now available in a more general form through the [fix
adapt](fix_adapt) command, as explained below. Note that if you use an
old input script and specify Astart and Astop without a cutoff, then
LAMMPS will interpret that as A and a cutoff, which is probably not what
you want.
::::

The [fix adapt](fix_adapt) command can be used to vary A for one or more
pair types over the course of a simulation, in which case pair_coeff
settings for A must still be specified, but will be overridden. For
example these commands will vary the prefactor A for all pairwise
interactions from 0.0 at the beginning to 30.0 at the end of a run:

``` LAMMPS
variable prefactor equal ramp(0,30)
fix 1 all adapt 1 pair soft a * * v_prefactor
```

Note that a formula defined by an [equal-style variable](variable) can
use the current timestep, elapsed time in the current run, elapsed time
since the beginning of a series of runs, as well as access other
variables.

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

For atom type pairs I,J and I != J, the A coefficient and cutoff
distance for this pair style can be mixed. A is always mixed via a
*geometric* rule. The cutoff is mixed according to the pair_modify mix
value. The default mix value is *geometric*. See the \"pair_modify\"
command for details.

This pair style does not support the [pair_modify](pair_modify) shift
option, since the pair interaction goes to 0.0 at the cutoff.

The [pair_modify](pair_modify) table and tail options are not relevant
for this pair style.

This pair style writes its information to [binary restart
files](restart), so pair_style and pair_coeff commands do not need to be
specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
[run_style respa](run_style) command. It does not support the *inner*,
*middle*, *outer* keywords.

------------------------------------------------------------------------

## Restrictions

> none

## Related commands

[pair_coeff](pair_coeff), [fix nve/limit](fix_nve_limit), [fix
adapt](fix_adapt)

## Default

none
