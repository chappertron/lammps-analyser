# timestep command

## Syntax

    timestep dt

-   dt = timestep size (time units)

## Examples

``` LAMMPS
timestep 2.0
timestep 0.003
```

## Description

Set the timestep size for subsequent molecular dynamics simulations. See
the [units](units) command for the time units associated with each
choice of units that LAMMPS supports.

The default value for the timestep size also depends on the choice of
units for the simulation; see the default values below.

When the [run style](run_style) is *respa*, dt is the timestep for the
outer loop (largest) timestep.

## Restrictions

> none

## Related commands

[fix dt/reset](fix_dt_reset), [run](run), [run_style](run_style) respa,
[units](units)

## Default

  ------------------------------- --------------- -----------------------
  choice of [units](units) \|                     
  time units \| default timestep                  
  size \|                                         

  lj                              $\tau$          0.005 $\tau$

  real                            fs              1.0 fs

  metal                           ps              0.001 ps

  si                              s               1.0e-8 s (10 ns)

  cgs                             s               1.0e-8 s (10 ns)

  electron                        fs              0.001 fs

  micro                           $\mu$s          2.0 $\mu$s

  nano                            ns              0.00045 ns
  ------------------------------- --------------- -----------------------
