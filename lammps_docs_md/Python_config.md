# Configuration information

The following methods can be used to query the LAMMPS library about
compile time settings and included packages and styles.

``` {.python caption="Example for using configuration settings functions"}
from lammps import lammps

lmp = lammps()

try:
    lmp.file("in.missing")
except Exception as e:
    print("LAMMPS failed with error:", e)

# write compressed dump file depending on available of options

if lmp.has_style("dump", "atom/zstd"):
    lmp.command("dump d1 all atom/zstd 100 dump.zst")
elif lmp.has_style("dump", "atom/gz"):
    lmp.command("dump d1 all atom/gz 100 dump.gz")
elif lmp.has_gzip_support():
    lmp.command("dump d1 all atom 100 dump.gz")
else:
    lmp.command("dump d1 all atom 100 dump")
```

------------------------------------------------------------------------

**Methods:**

-   `lammps.has_mpi_support <lammps.lammps.has_mpi_support>`{.interpreted-text
    role="py:attr"}\_\_
-   `lammps.has_exceptions <lammps.lammps.has_exceptions>`{.interpreted-text
    role="py:attr"}\_\_
-   `lammps.has_gzip_support <lammps.lammps.has_gzip_support>`{.interpreted-text
    role="py:attr"}\_\_
-   `lammps.has_png_support <lammps.lammps.has_png_support>`{.interpreted-text
    role="py:attr"}\_\_
-   `lammps.has_jpeg_support <lammps.lammps.has_jpeg_support>`{.interpreted-text
    role="py:attr"}\_\_
-   `lammps.has_ffmpeg_support <lammps.lammps.has_ffmpeg_support>`{.interpreted-text
    role="py:attr"}\_\_
-   `lammps.installed_packages <lammps.lammps.installed_packages>`{.interpreted-text
    role="py:attr"}\_\_
-   `lammps.get_accelerator_config <lammps.lammps.accelerator_config>`{.interpreted-text
    role="py:meth"}\_\_
-   `lammps.has_style() <lammps.lammps.has_style()>`{.interpreted-text
    role="py:meth"}\_\_
-   `lammps.available_styles() <lammps.lammps.available_styles()>`{.interpreted-text
    role="py:meth"}\_\_
