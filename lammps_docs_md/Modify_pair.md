# Pair styles

Classes that compute pairwise non-bonded interactions are derived from
the `Pair` class. In LAMMPS, pairwise force calculations include
many-body potentials such as EAM, Tersoff, or ReaxFF where particles
interact without an explicit bond topology but include interactions
beyond pairwise non-bonded contributions. New styles can be created to
add support for additional pair potentials to LAMMPS. When the
modifications are small, sometimes it is more effective to derive from
an existing pair style class. This latter approach is also used by
[Accelerator packages](Speed_packages) where the accelerated style names
differ from their base classes by an appended suffix.

The file `src/pair_lj_cut.cpp` is an example of a Pair class with a very
simple potential function. It includes several optional methods to
enable its use with [run_style respa](run_style) and [compute
group/group](compute_group_group). [Developer_write_pair]{.title-ref}
contains a detailed discussion of writing new pair styles from scratch,
and how simple and more complex pair styles can be implemented with
examples from existing pair styles.

Here is a brief list of some the class methods in the Pair class that
*must* be or *may* be overridden in a derived class for a new pair
style.

  -----------------------------------------------------------------------
  Required                \"pure\" methods that *must* be overridden in a
                          derived class
  ----------------------- -----------------------------------------------
  compute                 workhorse routine that computes pairwise
                          interactions

  settings                processes the arguments to the pair_style
                          command

  coeff                   set coefficients for one i,j type pair, called
                          from pair_coeff
  -----------------------------------------------------------------------

  -----------------------------------------------------------------------------
  Optional                     methods that have a default or dummy
                               implementation
  ---------------------------- ------------------------------------------------
  init_one                     perform initialization for one i,j type pair

  init_style                   style initialization: request neighbor list(s),
                               error checks

  init_list                    Neighbor class callback function to pass
                               neighbor list to pair style

  single                       force/r and energy of a single pairwise
                               interaction between 2 atoms

  compute_inner/middle/outer   versions of compute used by rRESPA

  memory_usage                 return estimated amount of memory used by the
                               pair style

  modify_params                process arguments to pair_modify command

  extract                      provide access to internal scalar or per-type
                               data like cutoffs

  extract_peratom              provide access to internal per-atom data

  setup                        initialization at the beginning of a run

  finish                       called at the end of a run, e.g. to print

  write & read_restart         write/read i,j pair coeffs to restart files

  write &                      write/read global settings to restart files
  read_restart_settings        

  write_data                   write Pair Coeffs section to data file

  write_data_all               write PairIJ Coeffs section to data file

  pack & unpack_forward_comm   copy data to and from buffer if style uses
                               forward communication

  pack & unpack_reverse_comm   copy data to and from buffer if style uses
                               reverse communication

  reinit                       reset all type-based parameters, called by fix
                               adapt for example

  reset_dt                     called when the time step is changed by timestep
                               or fix reset/dt
  -----------------------------------------------------------------------------

Here is a list of flags or settings that should be set in the
constructor of the derived pair class when they differ from the default
setting.

  -----------------------------------------------------------------------------
  Name of flag              Description                               default
  ------------------------- ----------------------------------------- ---------
  single_enable             1 if single() method is implemented, 0 if 1
                            missing                                   

  respa_enable              1 if pair style has                       0
                            compute_inner/middle/outer()              

  restartinfo               1 if pair style writes its settings to a  1
                            restart                                   

  one_coeff                 1 if only a pair_coeff \* \* command is   0
                            allowed                                   

  manybody_flag             1 if pair style is a manybody potential   0

  unit_convert_flag         value != 0 indicates support for unit     0
                            conversion                                

  no_virial_fdotr_compute   1 if pair style does not call             0
                            virial_fdotr_compute()                    

  writedata                 1 if write_data() and write_data_all()    0
                            are implemented                           

  comm_forward              size of buffer (in doubles) for forward   0
                            communication                             

  comm_reverse              size of buffer (in doubles) for reverse   0
                            communication                             

  ghostneigh                1 if cutghost is set and style uses       0
                            neighbors of ghosts                       

  finitecutflag             1 if cutoff depends on diameter of atoms  0

  reinitflag                1 if style has reinit() and is compatible 0
                            with fix adapt                            

  ewaldflag                 1 if compatible with kspace_style ewald   0

  pppmflag                  1 if compatible with kspace_style pppm    0

  msmflag                   1 if compatible with kspace_style msm     0

  dispersionflag            1 if compatible with ewald/disp or        0
                            pppm/disp                                 

  tip4pflag                 1 if compatible with kspace_style         0
                            pppm/tip4p                                

  dipoleflag                1 if compatible with dipole kspace_style  0

  spinflag                  1 if compatible with spin kspace_style    0
  -----------------------------------------------------------------------------
