# Utility functions

To simplify some tasks, the library interface contains these utility
functions. They do not directly call the LAMMPS library.

-   `lammps_encode_image_flags`{.interpreted-text role="cpp:func"}
-   `lammps_decode_image_flags`{.interpreted-text role="cpp:func"}
-   `lammps_set_fix_external_callback`{.interpreted-text
    role="cpp:func"}
-   `lammps_fix_external_get_force`{.interpreted-text role="cpp:func"}
-   `lammps_fix_external_set_energy_global`{.interpreted-text
    role="cpp:func"}
-   `lammps_fix_external_set_energy_peratom`{.interpreted-text
    role="cpp:func"}
-   `lammps_fix_external_set_virial_global`{.interpreted-text
    role="cpp:func"}
-   `lammps_fix_external_set_virial_peratom`{.interpreted-text
    role="cpp:func"}
-   `lammps_fix_external_set_vector_length`{.interpreted-text
    role="cpp:func"}
-   `lammps_fix_external_set_vector`{.interpreted-text role="cpp:func"}
-   `lammps_flush_buffers`{.interpreted-text role="cpp:func"}
-   `lammps_free`{.interpreted-text role="cpp:func"}
-   `lammps_is_running`{.interpreted-text role="cpp:func"}
-   `lammps_force_timeout`{.interpreted-text role="cpp:func"}
-   `lammps_has_error`{.interpreted-text role="cpp:func"}
-   `lammps_get_last_error_message`{.interpreted-text role="cpp:func"}
-   `lammps_python_api_version`{.interpreted-text role="cpp:func"}

The `lammps_free`{.interpreted-text role="cpp:func"} function is a
clean-up function to free memory that the library had allocated
previously via other function calls. Look for notes in the descriptions
of the individual commands where such memory buffers were allocated that
require the use of `lammps_free`{.interpreted-text role="cpp:func"}.

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_encode_image_flags
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_decode_image_flags(int image, int \*flags)
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_set_fix_external_callback(void *, const char*, FixExternalFnPtr,
void\*)
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_fix_external_get_force
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_fix_external_set_energy_global
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_fix_external_set_energy_peratom
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_fix_external_set_virial_global
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_fix_external_set_virial_peratom
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_fix_external_set_vector_length
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_fix_external_set_vector
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_flush_buffers
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_free
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_is_running
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_force_timeout
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_has_error
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_get_last_error_message
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_python_api_version
:::
