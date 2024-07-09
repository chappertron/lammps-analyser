# fix wall/lj93 command

Accelerator Variants: *wall/lj93/kk*

# fix wall/lj126 command

# fix wall/lj1043 command

# fix wall/colloid command

# fix wall/harmonic command

# fix wall/lepton command

# fix wall/morse command

# fix wall/table command

## Syntax

    fix ID group-ID style [tabstyle] [N] face args ... keyword value ...

-   ID, group-ID are documented in [fix](fix) command
-   style = *wall/lj93* or *wall/lj126* or *wall/lj1043* or
    *wall/colloid* or *wall/harmonic* or *wall/lepton* or *wall/morse*
    or *wall/table*
-   tabstyle = *linear* or *spline* = method of table interpolation
    (only applies to *wall/table*)
-   N = use N values in *linear* or *spline* interpolation (only applies
    to *wall/table*)
-   one or more face/arg pairs may be appended
-   face = *xlo* or *xhi* or *ylo* or *yhi* or *zlo* or *zhi*

<!-- -->

-   args for styles *lj93* or *lj126* or *lj1043* or *colloid* or
    *harmonic*

        args = coord epsilon sigma cutoff
        coord = position of wall = EDGE or constant or variable
          EDGE = current lo or hi edge of simulation box
          constant = number like 0.0 or -30.0 (distance units)
          variable = `equal-style variable <variable>`__ like v_x or v_wiggle
        epsilon = strength factor for wall-particle interaction (energy or energy/distance\^2 units)
          epsilon can be a variable (see below)
        sigma = size factor for wall-particle interaction (distance units)
          sigma can be a variable (see below)
        cutoff = distance from wall at which wall-particle interactions are cut off (distance units)

-   args for style *lepton*

        args = coord expression cutoff
        coord = position of wall = EDGE or constant or variable
          EDGE = current lo or hi edge of simulation box
          constant = number like 0.0 or -30.0 (distance units)
          variable = `equal-style variable <variable>`__ like v_x or v_wiggle
        expression = Lepton expression for the potential  (energy units)
        cutoff = distance from wall at which wall-particle interactions are cut off (distance units)

-   args for style *morse*

        args = coord D_0 alpha r_0 cutoff
        coord = position of wall = EDGE or constant or variable
          EDGE = current lo or hi edge of simulation box
          constant = number like 0.0 or -30.0 (distance units)
          variable = `equal-style variable <variable>`__ like v_x or v_wiggle
        D_0 = depth of the potential (energy units)
          D_0 can be a variable (see below)
        alpha = width factor for wall-particle interaction (1/distance units)
          alpha can be a variable (see below)
        r_0 = distance of the potential minimum from the face of region (distance units)
          r_0 can be a variable (see below)
        cutoff = distance from wall at which wall-particle interactions are cut off (distance units)

-   args for style *table*

        args = coord filename keyword cutoff
        coord = position of wall = EDGE or constant or variable
          EDGE = current lo or hi edge of simulation box
          constant = number like 0.0 or -30.0 (distance units)
          variable = `equal-style variable <variable>`__ like v_x or v_wiggle
        filename = file containing tabulated energy and force values
        keyword = section identifier to select a specific table in table file
        cutoff = distance from wall at which wall-particle interactions are cut off (distance units)

-   zero or more keyword/value pairs may be appended

-   keyword = *units* or *fld* or *pbc*

        *units* value = *lattice* or *box*
          *lattice* = the wall position is defined in lattice units
          *box* = the wall position is defined in simulation box units
        *fld* value = *yes* or *no*
          *yes* = invoke the wall constraint to be compatible with implicit FLD
          *no* = invoke the wall constraint in the normal way
        *pbc* value = *yes* or *no*
          *yes* = allow periodic boundary in a wall dimension
          *no* = require non-perioidic boundaries in any wall dimension

## Examples

``` LAMMPS
fix wallhi all wall/lj93 xlo -1.0 1.0 1.0 2.5 units box
fix wallhi all wall/lj93 xhi EDGE 1.0 1.0 2.5
fix wallhi all wall/harmonic xhi EDGE 100.0 0.0 4.0 units box
fix wallhi all wall/morse xhi EDGE 1.0 1.0 1.0 2.5 units box
fix wallhi all wall/lj126 v_wiggle 23.2 1.0 1.0 2.5
fix zwalls all wall/colloid zlo 0.0 1.0 1.0 0.858 zhi 40.0 1.0 1.0 0.858
fix xwall mobile wall/table spline 200 EDGE -5.0 walltab.dat HARMONIC 4.0
fix xwalls mobile wall/lepton xlo -5.0 "k*(r-rc)^2;k=100.0" 4.0 xhi 5.0 "k*(r-rc)^2;k=100.0" 4.0
```

## Description

Bound the simulation domain on one or more of its faces with a flat wall
that interacts with the atoms in the group by generating a force on the
atom in a direction perpendicular to the wall. The energy of
wall-particle interactions depends on the style.

For style *wall/lj93*, the energy E is given by the 9-3 Lennard-Jones
potential:

$$E = \epsilon \left[ \frac{2}{15} \left(\frac{\sigma}{r}\right)^{9} -
                      \left(\frac{\sigma}{r}\right)^3 \right]
                      \qquad r < r_c$$

For style *wall/lj126*, the energy E is given by the 12-6 Lennard-Jones
potential:

$$E = 4 \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} -
                      \left(\frac{\sigma}{r}\right)^6 \right]
                      \qquad r < r_c$$

For style *wall/lj1043*, the energy E is given by the 10-4-3
Lennard-Jones potential:

$$E = 2 \pi \epsilon \left[ \frac{2}{5} \left(\frac{\sigma}{r}\right)^{10} -
                      \left(\frac{\sigma}{r}\right)^4 -
                      \frac{\sqrt(2)\sigma^3}{3\left(r+\left(0.61/\sqrt(2)\right)\sigma\right)^3}\right]
                      \qquad r < r_c$$

For style *wall/colloid*, the energy E is given by an integrated form of
the [pair_style colloid](pair_colloid) potential:

$$\begin{aligned}
E = & \epsilon \left[ \frac{\sigma^{6}}{7560}
\left(\frac{6R-D}{D^{7}} + \frac{D+8R}{(D+2R)^{7}} \right) \right. \\
 & \left. - \frac{1}{6} \left(\frac{2R(D+R) + D(D+2R)
 \left[ \ln D - \ln (D+2R) \right]}{D(D+2R)} \right) \right] \qquad r < r_c
\end{aligned}$$

For style *wall/harmonic*, the energy E is given by a repulsive-only
harmonic spring potential:

$$E = \epsilon \quad (r - r_c)^2 \qquad r < r_c$$

For style *wall/morse*, the energy E is given by a Morse potential:

$$E = D_0 \left[ e^{- 2 \alpha (r - r_0)} - 2 e^{- \alpha (r - r_0)} \right]
    \qquad r < r_c$$

::: versionadded
28Mar2023
:::

For style *wall/lepton*, the energy E is provided as an Lepton
expression string using \"r\" as the distance variable. The [Lepton
library](https://simtk.org/projects/lepton)\_, that the *wall/lepton*
style interfaces with, evaluates this expression string at run time to
compute the wall-particle energy. It also creates an analytical
representation of the first derivative of this expression with respect
to \"r\" and then uses that to compute the force between the wall and
atoms in the fix group. The Lepton expression must be either enclosed in
quotes or must not contain any whitespace so that LAMMPS recognizes it
as a single keyword.

Optionally, the expression may use \"rc\" to refer to the cutoff
distance for the given wall. Further constants in the expression can be
defined in the same string as additional expressions separated by
semicolons. The expression \"k\*(r-rc)\^2;k=100.0\" represents a
repulsive-only harmonic spring as in fix *wall/harmonic* with a force
constant *K* (same as $\epsilon$ above) of 100 energy units. More
details on the Lepton expression strings are given below.

::: versionadded
28Mar2023
:::

For style *wall/table*, the energy E and forces are determined from
interpolation tables listed in one or more files as a function of
distance. The interpolation tables are used to evaluate energy and
forces between particles and the wall similar to how analytic formulas
are used for the other wall styles.

The interpolation tables are created as a pre-computation by fitting
cubic splines to the file values and interpolating energy and force
values at each of *N* distances. During a simulation, the tables are
used to interpolate energy and force values as needed for each wall and
particle separated by a distance *R*. The interpolation is done in one
of two styles: *linear* or *spline*.

For the *linear* style, the distance *R* is used to find the 2
surrounding table values from which an energy or force is computed by
linear interpolation.

For the *spline* style, cubic spline coefficients are computed and
stored for each of the *N* values in the table, one set of splines for
energy, another for force. Note that these splines are different than
the ones used to pre-compute the *N* values. Those splines were fit to
the *Nfile* values in the tabulated file, where often *Nfile* \< *N*.
The distance *R* is used to find the appropriate set of spline
coefficients which are used to evaluate a cubic polynomial which
computes the energy or force.

For each wall a filename and a keyword must be provided as in the
examples above. The filename specifies a file containing tabulated
energy and force values. The keyword specifies a section of the file.
The format of this file is described below.

In all cases, *r* is the distance from the particle to the wall at
position *coord*, and $r_c$ is the *cutoff* distance at which the
particle and wall no longer interact. The energy of the wall potential
is shifted so that the wall-particle interaction energy is 0.0 at the
cutoff distance.

Up to 6 walls or faces can be specified in a single command: *xlo*,
*xhi*, *ylo*, *yhi*, *zlo*, *zhi*. A *lo* face interacts with particles
near the lower side of the simulation box in that dimension. A *hi* face
interacts with particles near the upper side of the simulation box in
that dimension.

The position of each wall can be specified in one of 3 ways: as the EDGE
of the simulation box, as a constant value, or as a variable. If EDGE is
used, then the corresponding boundary of the current simulation box is
used. If a numeric constant is specified then the wall is placed at that
position in the appropriate dimension (x, y, or z). In both the EDGE and
constant cases, the wall will never move. If the wall position is a
variable, it should be specified as v_name, where name is an
[equal-style variable](variable) name. In this case the variable is
evaluated each timestep and the result becomes the current position of
the reflecting wall. Equal-style variables can specify formulas with
various mathematical functions, and include [thermo_style](thermo_style)
command keywords for the simulation box parameters and timestep and
elapsed time. Thus it is easy to specify a time-dependent wall position.
See examples below.

For the *wall/lj93* and *wall/lj126* and *wall/lj1043* styles,
$\epsilon$ and $\sigma$ are the usual Lennard-Jones parameters, which
determine the strength and size of the particle as it interacts with the
wall. Epsilon has energy units. Note that this $\epsilon$ and $\sigma$
may be different than any $\epsilon$ or $\sigma$ values defined for a
pair style that computes particle-particle interactions.

The *wall/lj93* interaction is derived by integrating over a 3d
half-lattice of Lennard-Jones 12/6 particles. The *wall/lj126*
interaction is effectively a harder, more repulsive wall interaction.
The *wall/lj1043* interaction is yet a different form of wall
interaction, described in Magda et al in [(Magda)](Magda).

For the *wall/colloid* style, *R* is the radius of the colloid particle,
*D* is the distance from the surface of the colloid particle to the wall
(r-R), and $\sigma$ is the size of a constituent LJ particle inside the
colloid particle and wall. Note that the cutoff distance Rc in this case
is the distance from the colloid particle center to the wall. The
prefactor $\epsilon$ can be thought of as an effective Hamaker constant
with energy units for the strength of the colloid-wall interaction. More
specifically, the $\epsilon$ prefactor is
$4\pi^2 \rho_{wall} \rho_{colloid} \epsilon \sigma^6$, where $\epsilon$
and $\sigma$ are the LJ parameters for the constituent LJ particles.
$\rho_{wall}$ and $\rho_{colloid}$ are the number density of the
constituent particles, in the wall and colloid respectively, in units of
1/volume.

The *wall/colloid* interaction is derived by integrating over
constituent LJ particles of size $\sigma$ within the colloid particle
and a 3d half-lattice of Lennard-Jones 12/6 particles of size $\sigma$
in the wall. As mentioned in the preceding paragraph, the density of
particles in the wall and colloid can be different, as specified by the
$\epsilon$ prefactor.

For the *wall/harmonic* style, $\epsilon$ is effectively the spring
constant K, and has units (energy/distance\^2). The input parameter
$\sigma$ is ignored. The minimum energy position of the harmonic spring
is at the *cutoff*. This is a repulsive-only spring since the
interaction is truncated at the *cutoff*

For the *wall/morse* style, the three parameters are in this order:
$D_0$ the depth of the potential, $\alpha$ the width parameter, and
$r_0$ the location of the minimum. $D_0$ has energy units, $\alpha$
inverse distance units, and $r_0$ distance units.

For any wall that supports them, the $\epsilon$ and/or $\sigma$ and/or
$\alpha$ parameter can be specified as an [equal-style
variable](variable), in which case it should be specified as v_name,
where name is the variable name. As with a variable wall position, the
variable is evaluated each timestep and the result becomes the current
epsilon or sigma of the wall. Equal-style variables can specify formulas
with various mathematical functions, and include
[thermo_style](thermo_style) command keywords for the simulation box
parameters and timestep and elapsed time. Thus it is easy to specify a
time-dependent wall interaction.

:::: note
::: title
Note
:::

For all of the styles, you must ensure that r is always \> 0 for all
particles in the group, or LAMMPS will generate an error. This means you
cannot start your simulation with particles at the wall position *coord*
(r = 0) or with particles on the wrong side of the wall (r \< 0). For
the *wall/lj93* and *wall/lj126* styles, the energy of the wall/particle
interaction (and hence the force on the particle) blows up as r -\> 0.
The *wall/colloid* style is even more restrictive, since the energy
blows up as D = r-R -\> 0. This means the finite-size particles of
radius R must be a distance larger than R from the wall position
*coord*. The *harmonic* style is a softer potential and does not blow up
as r -\> 0, but you must use a large enough $\epsilon$ that particles
always reamin on the correct side of the wall (r \> 0).
::::

The *units* keyword determines the meaning of the distance units used to
define a wall position, but only when a numeric constant or variable is
used. It is not relevant when EDGE is used to specify a face position.
In the variable case, the variable is assumed to produce a value
compatible with the *units* setting you specify.

A *box* value selects standard distance units as defined by the
[units](units) command, e.g. Angstroms for units = real or metal. A
*lattice* value means the distance units are in lattice spacings. The
[lattice](lattice) command must have been previously used to define the
lattice spacings.

The *fld* keyword can be used with a *yes* setting to invoke the wall
constraint before pairwise interactions are computed. This allows an
implicit FLD model using [pair_style lubricateU](pair_lubricateU) to
include the wall force in its calculations. If the setting is *no*, wall
forces are imposed after pairwise interactions, in the usual manner.

The *pbc* keyword can be used with a *yes* setting to allow walls to be
specified in a periodic dimension. See the [boundary](boundary) command
for options on simulation box boundaries. The default for *pbc* is *no*,
which means the system must be non-periodic when using a wall. But you
may wish to use a periodic box. E.g. to allow some particles to interact
with the wall via the fix group-ID, and others to pass through it and
wrap around a periodic box. In this case you should ensure that the wall
is sufficiently far enough away from the box boundary. If you do not,
then particles may interact with both the wall and with periodic images
on the other side of the box, which is probably not what you want.

------------------------------------------------------------------------

Here are examples of variable definitions that move the wall position in
a time-dependent fashion using equal-style [variables](variable). The
wall interaction parameters (epsilon, sigma) could be varied with
additional variable definitions.

``` LAMMPS
variable ramp equal ramp(0,10)
fix 1 all wall xlo v_ramp 1.0 1.0 2.5

variable linear equal vdisplace(0,20)
fix 1 all wall xlo v_linear 1.0 1.0 2.5

variable wiggle equal swiggle(0.0,5.0,3.0)
fix 1 all wall xlo v_wiggle 1.0 1.0 2.5

variable wiggle equal cwiggle(0.0,5.0,3.0)
fix 1 all wall xlo v_wiggle 1.0 1.0 2.5
```

The *ramp(lo,hi)* function adjusts the wall position linearly from *lo*
to *hi* over the course of a run. The *vdisplace(c0,velocity)* function
does something similar using the equation *position = c0 +
velocity\*delta*, where *delta* is the elapsed time.

The *swiggle(c0,A,period)* function causes the wall position to
oscillate sinusoidally according to this equation, where *omega = 2 PI /
period*:

    position = c0 + A sin(omega\*delta)

The *cwiggle(c0,A,period)* function causes the wall position to
oscillate sinusoidally according to this equation, which will have an
initial wall velocity of 0.0, and thus may impose a gentler perturbation
on the particles:

    position = c0 + A (1 - cos(omega\*delta))

------------------------------------------------------------------------

Lepton expression syntax and features
\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"

Lepton supports the following operators in expressions:

The following mathematical functions are available:

Numbers may be given in either decimal or exponential form. All of the
following are valid numbers: [5]{.title-ref}, [-3.1]{.title-ref},
[1e6]{.title-ref}, and [3.12e-2]{.title-ref}.

As an extension to the standard Lepton syntax, it is also possible to
use LAMMPS [variables](variable) in the format \"v_name\". Before
evaluating the expression, \"v_name\" will be replaced with the value of
the variable \"name\". This is compatible with all kinds of scalar
variables, but not with vectors, arrays, local, or per-atom variables.
If necessary, a custom scalar variable needs to be defined that can
access the desired (single) item from a non-scalar variable. As an
example, the following lines will instruct LAMMPS to ramp the force
constant for a harmonic bond from 100.0 to 200.0 during the next run:

``` LAMMPS
variable fconst equal ramp(100.0, 200)
bond_style lepton
bond_coeff 1 1.5 "v_fconst * (r^2)"
```

An expression may be followed by definitions for intermediate values
that appear in the expression. A semicolon \";\" is used as a delimiter
between value definitions. For example, the expression:

``` C
a^2+a*b+b^2; a=a1+a2; b=b1+b2
```

is exactly equivalent to

``` C
(a1+a2)^2+(a1+a2)*(b1+b2)+(b1+b2)^2
```

The definition of an intermediate value may itself involve other
intermediate values. Whitespace and quotation characters (\'\'\' and
\'\"\') are ignored. All uses of a value must appear *before* that
value\'s definition. For efficiency reasons, the expression string is
parsed, optimized, and then stored in an internal, pre-parsed
representation for evaluation.

Evaluating a Lepton expression is typically between 2.5 and 5 times
slower than the corresponding compiled and optimized C++ code. If
additional speed or GPU acceleration (via GPU or KOKKOS) is required,
the interaction can be represented as a table. Suitable table files can
be created either internally using the [pair_write](pair_write) or
[bond_write](bond_write) command or through the Python scripts in the
[tools/tabulate](tabulate) folder.

------------------------------------------------------------------------

## Table file format

Suitable tables for use with fix *wall/table* can be created by the
Python code in the `tools/tabulate` folder of the LAMMPS source code
distribution.

The format of a tabulated file is as follows (without the parenthesized
comments):

    # Tabulated wall potential UNITS: real

    HARMONIC                       (keyword is the first text on a line)
    N 100 FP 200 200
                                   (blank line)
    1  0.04    1568.16    792.00   (index, distance to wall, energy, force)
    2  0.08    1536.64    784.00
    3  0.12    1505.44    776.00
    ...
    99  3.96      0.16      8.00
    100  4.00     0         0

A section begins with a non-blank line whose first character is not a
\"#\"; blank lines or lines starting with \"#\" can be used as comments
between sections. The first line begins with a keyword which identifies
the section. The line can contain additional text, but the initial text
must match the argument specified in the fix *wall/table* command. The
next line lists (in any order) one or more parameters for the table.
Each parameter is a keyword followed by one or more numeric values.

The parameter \"N\" is required and its value is the number of table
entries that follow. Note that this may be different than the *N*
specified in the fix *wall/table* command. Let Ntable = *N* in the fix
command, and Nfile = \"N\" in the tabulated file. What LAMMPS does is a
preliminary interpolation by creating splines using the Nfile tabulated
values as nodal points. It uses these to interpolate as needed to
generate energy and force values at Ntable different points. The
resulting tables of length Ntable are then used as described above, when
computing energy and force for wall-particle interactions. This means
that if you want the interpolation tables of length Ntable to match
exactly what is in the tabulated file (with effectively no preliminary
interpolation), you should set Ntable = Nfile.

------------------------------------------------------------------------

## Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to [binary restart
files](restart).

The [fix_modify](fix_modify) *energy* option is supported by this fix to
add the energy of interaction between atoms and all the specified walls
to the global potential energy of the system as part of [thermodynamic
output](thermo_style). The default setting for this fix is [fix_modify
energy no](fix_modify).

The [fix_modify](fix_modify) *virial* option is supported by this fix to
add the contribution due to the interaction between atoms and all the
specified walls to both the global pressure and per-atom stress of the
system via the [compute pressure](compute_pressure) and [compute
stress/atom](compute_stress_atom) commands. The former can be accessed
by [thermodynamic output](thermo_style). The default setting for this
fix is [fix_modify virial no](fix_modify).

The [fix_modify](fix_modify) *respa* option is supported by this fix.
This allows to set at which level of the [r-RESPA](run_style) integrator
the fix is adding its forces. Default is the outermost level.

This fix computes a global scalar energy and a global vector of forces,
which can be accessed by various [output commands](Howto_output). Note
that the scalar energy is the sum of interactions with all defined
walls. If you want the energy on a per-wall basis, you need to use
multiple fix wall commands. The length of the vector is equal to the
number of walls defined by the fix. Each vector value is the normal
force on a specific wall. Note that an outward force on a wall will be a
negative value for *lo* walls and a positive value for *hi* walls. The
scalar and vector values calculated by this fix are \"extensive\".

No parameter of this fix can be used with the *start/stop* keywords of
the [run](run) command.

The forces due to this fix are imposed during an energy minimization,
invoked by the [minimize](minimize) command.

:::: note
::: title
Note
:::

If you want the atom/wall interaction energy to be included in the total
potential energy of the system (the quantity being minimized), you MUST
enable the [fix_modify](fix_modify) *energy* option for this fix.
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

## Restrictions

Fix *wall/lepton* is part of the LEPTON package and only enabled if
LAMMPS was built with this package. See the [Build
package](Build_package) page for more info.

## Related commands

[fix wall/reflect](fix_wall_reflect), [fix wall/gran](fix_wall_gran),
[fix wall/region](fix_wall_region)

## Default

The option defaults units = lattice, fld = no, and pbc = no.

------------------------------------------------------------------------

::: {#Magda}
**(Magda)** Magda, Tirrell, Davis, J Chem Phys, 83, 1888-1901 (1985);
erratum in JCP 84, 2901 (1986).
:::
