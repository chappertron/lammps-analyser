# Platform abstraction functions

The `platform` sub-namespace inside the `LAMMPS_NS` namespace provides a
collection of wrapper and convenience functions and utilities that
perform common tasks for which platform specific code would be required
or for which a more high-level abstraction would be convenient and
reduce duplicated code. This reduces redundant implementations and
encourages consistent behavior and thus has some overlap with the
[\"utils\" sub-namespace](Developer_utils).

## Time functions

::: {.doxygenfunction project="progguide"}
cputime
:::

::: {.doxygenfunction project="progguide"}
walltime
:::

::: {.doxygenfunction project="progguide"}
usleep
:::

## Platform information functions

::: {.doxygenfunction project="progguide"}
os_info
:::

::: {.doxygenfunction project="progguide"}
compiler_info
:::

::: {.doxygenfunction project="progguide"}
cxx_standard
:::

::: {.doxygenfunction project="progguide"}
openmp_standard
:::

::: {.doxygenfunction project="progguide"}
mpi_vendor
:::

::: {.doxygenfunction project="progguide"}
mpi_info
:::

::: {.doxygenfunction project="progguide"}
compress_info
:::

## File and path functions and global constants

::: {.doxygenvariable project="progguide"}
filepathsep
:::

::: {.doxygenvariable project="progguide"}
pathvarsep
:::

::: {.doxygenfunction project="progguide"}
guesspath
:::

::: {.doxygenfunction project="progguide"}
path_basename
:::

::: {.doxygenfunction project="progguide"}
path_join
:::

::: {.doxygenfunction project="progguide"}
file_is_readable
:::

::: {.doxygenfunction project="progguide"}
is_console
:::

::: {.doxygenfunction project="progguide"}
path_is_directory
:::

::: {.doxygenfunction project="progguide"}
current_directory
:::

::: {.doxygenfunction project="progguide"}
list_directory
:::

::: {.doxygenfunction project="progguide"}
chdir
:::

::: {.doxygenfunction project="progguide"}
mkdir
:::

::: {.doxygenfunction project="progguide"}
rmdir
:::

::: {.doxygenfunction project="progguide"}
unlink
:::

## Standard I/O function wrappers

::: {.doxygenvariable project="progguide"}
END_OF_FILE
:::

::: {.doxygenfunction project="progguide"}
ftell
:::

::: {.doxygenfunction project="progguide"}
fseek
:::

::: {.doxygenfunction project="progguide"}
ftruncate
:::

::: {.doxygenfunction project="progguide"}
popen
:::

::: {.doxygenfunction project="progguide"}
pclose
:::

## Environment variable functions

::: {.doxygenfunction project="progguide"}
putenv
:::

::: {.doxygenfunction project="progguide"}
unsetenv
:::

::: {.doxygenfunction project="progguide"}
list_pathenv
:::

::: {.doxygenfunction project="progguide"}
find_exe_path
:::

## Dynamically loaded object or library functions

::: {.doxygenfunction project="progguide"}
dlopen
:::

::: {.doxygenfunction project="progguide"}
dlclose
:::

::: {.doxygenfunction project="progguide"}
dlsym
:::

::: {.doxygenfunction project="progguide"}
dlerror
:::

## Compressed file I/O functions

::: {.doxygenfunction project="progguide"}
has_compress_extension
:::

::: {.doxygenfunction project="progguide"}
compressed_read
:::

::: {.doxygenfunction project="progguide"}
compressed_write
:::
