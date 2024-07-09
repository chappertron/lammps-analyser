# Output from LAMMPS (thermo, dumps, computes, fixes, variables)

There are four basic kinds of LAMMPS output:

-   [Thermodynamic output](thermo_style), which is a list of quantities
    printed every few timesteps to the screen and logfile.
-   [Dump files](dump), which contain snapshots of atoms and various
    per-atom values and are written at a specified frequency.
-   Certain fixes can output user-specified quantities to files: [fix
    ave/time](fix_ave_time) for time averaging, [fix
    ave/chunk](fix_ave_chunk) for spatial or other averaging, and [fix
    print](fix_print) for single-line output of [variables](variable).
    Fix print can also output to the screen.
-   [Restart files](restart).

A simulation prints one set of thermodynamic output and (optionally)
restart files. It can generate any number of dump files and fix output
files, depending on what [dump](dump) and [fix](fix) commands you
specify.

As discussed below, LAMMPS gives you a variety of ways to determine what
quantities are computed and printed when the thermodynamics, dump, or
fix commands listed above perform output. Throughout this discussion,
note that users can also [add their own computes and fixes to
LAMMPS](Modify) which can then generate values that can then be output
with these commands.

The following subsections discuss different LAMMPS commands related to
output and the kind of data they operate on and produce:

-   [Global/per-atom/local/per-grid data](global)
-   [Scalar/vector/array data](scalar)
-   [Per-grid data](grid)
-   [Disambiguation](disambiguation)
-   [Thermodynamic output](thermo)
-   [Dump file output](dump)
-   [Fixes that write output files](fixoutput)
-   [Computes that process output quantities](computeoutput)
-   [Fixes that process output quantities](fixprocoutput)
-   [Computes that generate values to output](compute)
-   [Fixes that generate values to output](fix)
-   [Variables that generate values to output](variable)
-   [Summary table of output options and data flow between
    commands](table)

## Global/per-atom/local/per-grid data {#global}

Various output-related commands work with four different styles of data:
global, per-atom, local, and per-grid. A global datum is one or more
system-wide values, e.g. the temperature of the system. A per-atom datum
is one or more values per atom, e.g. the kinetic energy of each atom.
Local datums are calculated by each processor based on the atoms it
owns, but there may be zero or more per atom, e.g. a list of bond
distances.

A per-grid datum is one or more values per grid cell, for a grid which
overlays the simulation domain. The grid cells and the data they store
are distributed across processors; each processor owns the grid cells
whose center point falls within its subdomain.

## Scalar/vector/array data {#scalar}

Global, per-atom, and local datums can come in three kinds: a single
scalar value, a vector of values, or a 2d array of values. The doc page
for a \"compute\" or \"fix\" or \"variable\" that generates data will
specify both the style and kind of data it produces, e.g. a per-atom
vector.

When a quantity is accessed, as in many of the output commands discussed
below, it can be referenced via the following bracket notation, where ID
in this case is the ID of a compute. The leading \"c\_\" would be
replaced by \"f\_\" for a fix, or \"v\_\" for a variable:

  ---------------- --------------------------------------------
  c_ID             entire scalar, vector, or array

  c_ID\[I\]        one element of vector, one column of array

  c_ID\[I\]\[J\]   one element of array
  ---------------- --------------------------------------------

In other words, using one bracket reduces the dimension of the data once
(vector -\> scalar, array -\> vector). Using two brackets reduces the
dimension twice (array -\> scalar). Thus a command that uses scalar
values as input can typically also process elements of a vector or
array.

## Per-grid data {#grid}

Per-grid data can come in two kinds: a vector of values (one per grid
cekk), or a 2d array of values (multiple values per grid ckk). The doc
page for a \"compute\" or \"fix\" that generates data will specify names
for both the grid(s) and datum(s) it produces, e.g. per-grid vectors or
arrays, which can be referenced by other commands. See the [Howto
grid](Howto_grid) doc page for more details.

## Disambiguation

Some computes and fixes produce data in multiple styles, e.g. a global
scalar and a per-atom vector. Usually the context in which the input
script references the data determines which style is meant. Example: if
a compute provides both a global scalar and a per-atom vector, the
former will be accessed by using `c_ID` in an equal-style variable,
while the latter will be accessed by using `c_ID` in an atom-style
variable. Note that atom-style variable formulas can also access global
scalars, but in this case it is not possible to do this directly because
of the ambiguity. Instead, an equal-style variable can be defined which
accesses the global scalar, and that variable can be used in the
atom-style variable formula in place of `c_ID`.

## Thermodynamic output {#thermo}

The frequency and format of thermodynamic output is set by the
[thermo](thermo), [thermo_style](thermo_style), and
[thermo_modify](thermo_modify) commands. The
[thermo_style](thermo_style) command also specifies what values are
calculated and written out. Pre-defined keywords can be specified (e.g.
press, etotal, etc). Three additional kinds of keywords can also be
specified (c_ID, f_ID, v_name), where a [compute](compute) or [fix](fix)
or [variable](variable) provides the value to be output. In each case,
the compute, fix, or variable must generate global values for input to
the [thermo_style custom](dump) command.

Note that thermodynamic output values can be \"extensive\" or
\"intensive\". The former scale with the number of atoms in the system
(e.g. total energy), the latter do not (e.g. temperature). The setting
for [thermo_modify norm](thermo_modify) determines whether extensive
quantities are normalized or not. Computes and fixes produce either
extensive or intensive values; see their individual doc pages for
details. [Equal-style variables](variable) produce only intensive
values; you can include a division by \"natoms\" in the formula if
desired, to make an extensive calculation produce an intensive result.

## Dump file output {#dump}

Dump file output is specified by the [dump](dump) and
[dump_modify](dump_modify) commands. There are several pre-defined
formats (dump atom, dump xtc, etc).

There is also a [dump custom](dump) format where the user specifies what
values are output with each atom. Pre-defined atom attributes can be
specified (id, x, fx, etc). Three additional kinds of keywords can also
be specified (c_ID, f_ID, v_name), where a [compute](compute) or
[fix](fix) or [variable](variable) provides the values to be output. In
each case, the compute, fix, or variable must generate per-atom values
for input to the [dump custom](dump) command.

There is also a [dump local](dump) format where the user specifies what
local values to output. A pre-defined index keyword can be specified to
enumerate the local values. Two additional kinds of keywords can also be
specified (c_ID, f_ID), where a [compute](compute) or [fix](fix) or
[variable](variable) provides the values to be output. In each case, the
compute or fix must generate local values for input to the [dump
local](dump) command.

There is also a [dump grid](dump) format where the user specifies what
per-grid values to output from computes or fixes that generate per-grid
data.

## Fixes that write output files {#fixoutput}

Several fixes take various quantities as input and can write output
files: [fix ave/time](fix_ave_time), [fix ave/chunk](fix_ave_chunk),
[fix ave/histo](fix_ave_histo), [fix ave/correlate](fix_ave_correlate),
and [fix print](fix_print).

The [fix ave/time](fix_ave_time) command enables direct output to a file
and/or time-averaging of global scalars or vectors. The user specifies
one or more quantities as input. These can be global [compute](compute)
values, global [fix](fix) values, or [variables](variable) of any style
except the atom style which produces per-atom values. Since a variable
can refer to keywords used by the [thermo_style custom](thermo_style)
command (like temp or press) and individual per-atom values, a wide
variety of quantities can be time averaged and/or output in this way. If
the inputs are one or more scalar values, then the fix generate a global
scalar or vector of output. If the inputs are one or more vector values,
then the fix generates a global vector or array of output. The
time-averaged output of this fix can also be used as input to other
output commands.

The [fix ave/chunk](fix_ave_chunk) command enables direct output to a
file of chunk-averaged per-atom quantities like those output in dump
files. Chunks can represent spatial bins or other collections of atoms,
e.g. individual molecules. The per-atom quantities can be atom density
(mass or number) or atom attributes such as position, velocity, force.
They can also be per-atom quantities calculated by a [compute](compute),
by a [fix](fix), or by an atom-style [variable](variable). The
chunk-averaged output of this fix is global and can also be used as
input to other output commands.

Note that the [fix ave/grid](fix_ave_grid) command can also average the
same per-atom quantities within spatial bins, but it does this for a
distributed grid whose grid cells are owned by different processors. It
outputs per-grid data, not global data, so it is more efficient for
large numbers of averaging bins.

The [fix ave/histo](fix_ave_histo) command enables direct output to a
file of histogrammed quantities, which can be global or per-atom or
local quantities. The histogram output of this fix can also be used as
input to other output commands.

The [fix ave/correlate](fix_ave_correlate) command enables direct output
to a file of time-correlated quantities, which can be global values. The
correlation matrix output of this fix can also be used as input to other
output commands.

The [fix print](fix_print) command can generate a line of output written
to the screen and log file or to a separate file, periodically during a
running simulation. The line can contain one or more
[variable](variable) values for any style variable except the vector or
atom styles). As explained above, variables themselves can contain
references to global values generated by [thermodynamic
keywords](thermo_style), [computes](compute), [fixes](fix), or other
[variables](variable), or to per-atom values for a specific atom. Thus
the [fix print](fix_print) command is a means to output a wide variety
of quantities separate from normal thermodynamic or dump file output.

## Computes that process output quantities {#computeoutput}

The [compute reduce](compute_reduce) and [compute
reduce/region](compute_reduce) commands take one or more per-atom or
local vector quantities as inputs and \"reduce\" them (sum, min, max,
ave) to scalar quantities. These are produced as output values which can
be used as input to other output commands.

The [compute slice](compute_slice) command take one or more global
vector or array quantities as inputs and extracts a subset of their
values to create a new vector or array. These are produced as output
values which can be used as input to other output commands.

The [compute property/atom](compute_property_atom) command takes a list
of one or more pre-defined atom attributes (id, x, fx, etc) and stores
the values in a per-atom vector or array. These are produced as output
values which can be used as input to other output commands. The list of
atom attributes is the same as for the [dump custom](dump) command.

The [compute property/local](compute_property_local) command takes a
list of one or more pre-defined local attributes (bond info, angle info,
etc) and stores the values in a local vector or array. These are
produced as output values which can be used as input to other output
commands.

The [compute property/grid](compute_property_grid) command takes a list
of one or more pre-defined per-grid attributes (id, grid cell coords,
etc) and stores the values in a per-grid vector or array. These are
produced as output values which can be used as input to the [dump
grid](dump) command.

The [compute property/chunk](compute_property_chunk) command takes a
list of one or more pre-defined chunk attributes (id, count, coords for
spatial bins) and stores the values in a global vector or array. These
are produced as output values which can be used as input to other output
commands.

## Fixes that process output quantities {#fixprocoutput}

The [fix vector](fix_vector) command can create global vectors as output
from global scalars as input, accumulating them one element at a time.

The [fix ave/atom](fix_ave_atom) command performs time-averaging of
per-atom vectors. The per-atom quantities can be atom attributes such as
position, velocity, force. They can also be per-atom quantities
calculated by a [compute](compute), by a [fix](fix), or by an atom-style
[variable](variable). The time-averaged per-atom output of this fix can
be used as input to other output commands.

The [fix store/state](fix_store_state) command can archive one or more
per-atom attributes at a particular time, so that the old values can be
used in a future calculation or output. The list of atom attributes is
the same as for the [dump custom](dump) command, including per-atom
quantities calculated by a [compute](compute), by a [fix](fix), or by an
atom-style [variable](variable). The output of this fix can be used as
input to other output commands.

The [fix ave/grid](fix_ave_grid) command performs time-averaging of
either per-atom or per-grid data.

For per-atom data it performs averaging for the atoms within each grid
cell, similar to the [fix ave/chunk](fix_ave_chunk) command when its
chunks are defined as regular 2d or 3d bins. The per-atom quantities can
be atom density (mass or number) or atom attributes such as position,
velocity, force. They can also be per-atom quantities calculated by a
[compute](compute), by a [fix](fix), or by an atom-style
[variable](variable).

The chief difference between the [fix ave/grid](fix_ave_grid) and [fix
ave/chunk](fix_ave_chunk) commands when used in this context is that the
former uses a distributed grid, while the latter uses a global grid.
Distributed means that each processor owns the subset of grid cells
within its subdomain. Global means that each processor owns a copy of
the entire grid. The [fix ave/grid](fix_ave_grid) command is thus more
efficient for large grids.

For per-grid data, the [fix ave/grid](fix_ave_grid) command takes inputs
for grid data produced by other computes or fixes and averages the
values for each grid point over time.

## Computes that generate values to output {#compute}

Every [compute](compute) in LAMMPS produces either global or per-atom or
local or per-grid values. The values can be scalars or vectors or arrays
of data. These values can be output using the other commands described
in this section. The page for each compute command describes what it
produces. Computes that produce per-atom or local or per-grid values
have the word \"atom\" or \"local\" or \"grid as the last word in their
style name. Computes without the word \"atom\" or \"local\" or \"grid\"
produce global values.

## Fixes that generate values to output {#fix}

Some [fixes](fix) in LAMMPS produces either global or per-atom or local
or per-grid values which can be accessed by other commands. The values
can be scalars or vectors or arrays of data. These values can be output
using the other commands described in this section. The page for each
fix command tells whether it produces any output quantities and
describes them.

## Variables that generate values to output {#variable}

[Variables](variable) defined in an input script can store one or more
strings. But equal-style, vector-style, and atom-style or atomfile-style
variables generate a global scalar value, global vector or values, or a
per-atom vector, respectively, when accessed. The formulas used to
define these variables can contain references to the thermodynamic
keywords and to global and per-atom data generated by computes, fixes,
and other variables. The values generated by variables can be used as
input to and thus output by the other commands described in this
section.

Per-grid variables have not (yet) been implemented.

## Summary table of output options and data flow between commands {#table}

This table summarizes the various commands that can be used for
generating output from LAMMPS. Each command produces output data of some
kind and/or writes data to a file. Most of the commands can take data
from other commands as input. Thus you can link many of these commands
together in pipeline form, where data produced by one command is used as
input to another command and eventually written to the screen or to a
file. Note that to hook two commands together the output and input data
types must match, e.g. global/per-atom/local data and
scalar/vector/array data.

Also note that, as described above, when a command takes a scalar as
input, that could be an element of a vector or array. Likewise a vector
input could be a column of an array.

  ----------------------------------------- --------------------- -----------------------
  Command                                   Input                 Output

  [thermo_style custom](thermo_style) \|                          
  global scalars \| screen, log file \|                           

  [dump custom](dump) \| per-atom vectors                         
  \| dump file \|                                                 

  [dump local](dump) \| local vectors \|                          
  dump file \|                                                    

  [dump grid](dump) \| per-grid vectors \|                        
  dump file \|                                                    

  [fix print](fix_print) \| global scalar                         
  from variable \| screen, file \|                                

  [print](print) \| global scalar from                            
  variable \| screen \|                                           

  [computes](compute) \| N/A \|                                   
  global/per-atom/local/per-grid                                  
  scalar/vector/array \|                                          

  [fixes](fix) \| N/A \|                                          
  global/per-atom/local/per-grid                                  
  scalar/vector/array \|                                          

  [variables](variable) \| global scalars                         
  and vectors, per-atom vectors \| global                         
  scalar and vector, per-atom vector \|                           

  [compute reduce](compute_reduce) \|                             
  per-atom/local vectors \| global                                
  scalar/vector \|                                                

  [compute slice](compute_slice) \| global                        
  vectors/arrays \| global vector/array \|                        

  [compute                                                        
  property/atom](compute_property_atom) \|                        
  N/A \| per-atom vector/array \|                                 

  [compute                                                        
  property/local](compute_property_local)                         
  \| N/A \| local vector/array \|                                 

  [compute                                                        
  property/grid](compute_property_grid) \|                        
  N/A \| per-grid vector/array \|                                 

  [compute                                                        
  property/chunk](compute_property_chunk)                         
  \| N/A \| global vector/array \|                                

  [fix vector](fix_vector) \| global                              
  scalars \| global vector \|                                     

  [fix ave/atom](fix_ave_atom) \| per-atom                        
  vectors \| per-atom vector/array \|                             

  [fix ave/time](fix_ave_time) \| global                          
  scalars/vectors \| global                                       
  scalar/vector/array, file \|                                    

  [fix ave/chunk](fix_ave_chunk) \|                               
  per-atom vectors \| global array, file \|                       

  [fix ave/grid](fix_ave_grid) \| per-atom                        
  vectors or per-grid vectors \| per-grid                         
  vector/array \|                                                 

  [fix ave/histo](fix_ave_histo) \|                               
  global/per-atom/local scalars and vectors                       
  \| global array, file \|                                        

  [fix ave/correlate](fix_ave_correlate) \|                       
  global scalars \| global array, file \|                         

  [fix store/state](fix_store_state) \|                           
  per-atom vectors \| per-atom vector/array                       
  \|                                                              
  ----------------------------------------- --------------------- -----------------------
