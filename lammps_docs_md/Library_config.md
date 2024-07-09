# Configuration information

This section documents the following functions:

-   `lammps_version`{.interpreted-text role="cpp:func"}
-   `lammps_get_os_info`{.interpreted-text role="cpp:func"}
-   `lammps_config_has_mpi_support`{.interpreted-text role="cpp:func"}
-   `lammps_config_has_gzip_support`{.interpreted-text role="cpp:func"}
-   `lammps_config_has_png_support`{.interpreted-text role="cpp:func"}
-   `lammps_config_has_jpeg_support`{.interpreted-text role="cpp:func"}
-   `lammps_config_has_ffmpeg_support`{.interpreted-text
    role="cpp:func"}
-   `lammps_config_has_exceptions`{.interpreted-text role="cpp:func"}
-   `lammps_config_has_package`{.interpreted-text role="cpp:func"}
-   `lammps_config_package_count`{.interpreted-text role="cpp:func"}
-   `lammps_config_package_name`{.interpreted-text role="cpp:func"}
-   `lammps_config_accelerator`{.interpreted-text role="cpp:func"}
-   `lammps_has_gpu_device`{.interpreted-text role="cpp:func"}
-   `lammps_gpu_device_info`{.interpreted-text role="cpp:func"}
-   `lammps_has_style`{.interpreted-text role="cpp:func"}
-   `lammps_style_count`{.interpreted-text role="cpp:func"}
-   `lammps_style_name`{.interpreted-text role="cpp:func"}
-   `lammps_has_id`{.interpreted-text role="cpp:func"}
-   `lammps_id_count`{.interpreted-text role="cpp:func"}
-   `lammps_id_name`{.interpreted-text role="cpp:func"}

------------------------------------------------------------------------

These library functions can be used to query the LAMMPS library for
compile time settings and included packages and styles. This enables
programs that use the library interface to determine whether the linked
LAMMPS library is compatible with the requirements of the application
without crashing during the LAMMPS functions (e.g. due to missing pair
styles from packages) or to choose between different options (e.g.
whether to use `lj/cut`, `lj/cut/opt`, `lj/cut/omp` or `lj/cut/intel`).
Most of the functions can be called directly without first creating a
LAMMPS instance. While crashes within LAMMPS may be recovered from by
enabling [exceptions](exceptions), avoiding them proactively is a safer
approach.

``` {.c caption="Example for using configuration settings functions"}
#include "library.h"
#include <stdio.h>

int main(int argc, char **argv)
{
  void *handle;

  handle = lammps_open_no_mpi(0, NULL, NULL);
  lammps_file(handle, "in.missing");
  if (lammps_has_error(handle)) {
    char errmsg[256];
    int errtype;
    errtype = lammps_get_last_error_message(handle, errmsg, 256);
    fprintf(stderr, "LAMMPS failed with error: %s\n", errmsg);
    return 1;
  }
  /* write compressed dump file depending on available of options */
  if (lammps_has_style(handle, "dump", "atom/zstd")) {
    lammps_command(handle, "dump d1 all atom/zstd 100 dump.zst");
  } else if (lammps_has_style(handle, "dump", "atom/gz")) {
    lammps_command(handle, "dump d1 all atom/gz 100 dump.gz");
  } else if (lammps_config_has_gzip_support()) {
    lammps_command(handle, "dump d1 all atom 100 dump.gz");
  } else {
    lammps_command(handle, "dump d1 all atom 100 dump");
  }
  lammps_close(handle);
  return 0;
}
```

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_version
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_get_os_info
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_config_has_mpi_support
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_config_has_gzip_support
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_config_has_png_support
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_config_has_jpeg_support
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_config_has_ffmpeg_support
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_config_has_exceptions
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_config_has_package
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_config_package_count
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_config_package_name
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_config_accelerator
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_has_gpu_device
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_get_gpu_device_info
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_has_style
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_style_count
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_style_name
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_has_id
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_id_count
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_id_name
:::
