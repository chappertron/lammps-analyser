# Utility functions

The `utils` sub-namespace inside the `LAMMPS_NS` namespace provides a
collection of convenience functions and utilities that perform common
tasks that are required repeatedly throughout the LAMMPS code like
reading or writing to files with error checking or translation of
strings into specific types of numbers with checking for validity. This
reduces redundant implementations and encourages consistent behavior and
thus has some overlap with the [\"platform\"
sub-namespace](Developer_platform).

## I/O with status check and similar functions

The the first two functions are wrappers around the corresponding C
library calls `fgets()` or `fread()`. They will check if there were
errors on reading or an unexpected end-of-file state was reached. In
that case, the functions will stop with an error message, indicating the
name of the problematic file, if possible unless the *error* argument is
a NULL pointer.

The
`utils::fgets_trunc() <LAMMPS_NS::utils::fgets_trunc>`{.interpreted-text
role="cpp:func"}\_\_ function will work similar for `fgets()` but it
will read in a whole line (i.e. until the end of line or end of file),
but store only as many characters as will fit into the buffer including
a final newline character and the terminating NULL byte. If the line in
the file is longer it will thus be truncated in the buffer. This
function is used by
`utils::read_lines_from_file() <LAMMPS_NS::utils::read_lines_from_file>`{.interpreted-text
role="cpp:func"}\_\_ to read individual lines but make certain they
follow the size constraints.

The
`utils::read_lines_from_file() <LAMMPS_NS::utils::read_lines_from_file>`{.interpreted-text
role="cpp:func"}\_\_ function will read the requested number of lines of
a maximum length into a buffer and will return 0 if successful or 1 if
not. It also guarantees that all lines are terminated with a newline
character and the entire buffer with a NULL character.

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
sfgets
:::

::: {.doxygenfunction project="progguide"}
sfread
:::

::: {.doxygenfunction project="progguide"}
fgets_trunc
:::

::: {.doxygenfunction project="progguide"}
read_lines_from_file
:::

------------------------------------------------------------------------

## String to number conversions with validity check

These functions should be used to convert strings to numbers. They are
are strongly preferred over C library calls like `atoi()` or `atof()`
since they check if the **entire** string is a valid (floating-point or
integer) number, and will error out instead of silently returning the
result of a partial conversion or zero in cases where the string is not
a valid number. This behavior improves detecting typos or issues when
processing input files.

Similarly the
`utils::logical() <LAMMPS_NS::utils::logical>`{.interpreted-text
role="cpp:func"}\_\_ function will convert a string into a boolean and
will only accept certain words.

The *do_abort* flag should be set to `true` in case this function is
called only on a single MPI rank, as that will then trigger the a call
to `Error::one()` for errors instead of `Error::all()` and avoids a
\"hanging\" calculation when run in parallel.

Please also see
`utils::is_integer() <LAMMPS_NS::utils::is_integer>`{.interpreted-text
role="cpp:func"}\_\_ and
`utils::is_double() <LAMMPS_NS::utils::is_double>`{.interpreted-text
role="cpp:func"}\_\_ for testing strings for compliance without
conversion.

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
numeric(const char *file, int line, const std::string &str, bool
do_abort, LAMMPS*lmp)
:::

::: {.doxygenfunction project="progguide"}
numeric(const char *file, int line, const char*str, bool do_abort,
LAMMPS \*lmp)
:::

::: {.doxygenfunction project="progguide"}
inumeric(const char *file, int line, const std::string &str, bool
do_abort, LAMMPS*lmp)
:::

::: {.doxygenfunction project="progguide"}
inumeric(const char *file, int line, const char*str, bool do_abort,
LAMMPS \*lmp)
:::

::: {.doxygenfunction project="progguide"}
bnumeric(const char *file, int line, const std::string &str, bool
do_abort, LAMMPS*lmp)
:::

::: {.doxygenfunction project="progguide"}
bnumeric(const char *file, int line, const char*str, bool do_abort,
LAMMPS \*lmp)
:::

::: {.doxygenfunction project="progguide"}
tnumeric(const char *file, int line, const std::string &str, bool
do_abort, LAMMPS*lmp)
:::

::: {.doxygenfunction project="progguide"}
tnumeric(const char *file, int line, const char*str, bool do_abort,
LAMMPS \*lmp)
:::

::: {.doxygenfunction project="progguide"}
logical(const char *file, int line, const std::string &str, bool
do_abort, LAMMPS*lmp)
:::

::: {.doxygenfunction project="progguide"}
logical(const char *file, int line, const char*str, bool do_abort,
LAMMPS \*lmp)
:::

## String processing

The following are functions to help with processing strings and parsing
files or arguments.

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
strdup
:::

::: {.doxygenfunction project="progguide"}
lowercase
:::

::: {.doxygenfunction project="progguide"}
uppercase
:::

::: {.doxygenfunction project="progguide"}
trim
:::

::: {.doxygenfunction project="progguide"}
trim_comment
:::

::: {.doxygenfunction project="progguide"}
strip_style_suffix
:::

::: {.doxygenfunction project="progguide"}
star_subst
:::

::: {.doxygenfunction project="progguide"}
has_utf8
:::

::: {.doxygenfunction project="progguide"}
utf8_subst
:::

::: {.doxygenfunction project="progguide"}
count_words(const char \*text)
:::

::: {.doxygenfunction project="progguide"}
count_words(const std::string &text)
:::

::: {.doxygenfunction project="progguide"}
count_words(const std::string &text, const std::string &separators)
:::

::: {.doxygenfunction project="progguide"}
trim_and_count_words
:::

::: {.doxygenfunction project="progguide"}
join_words
:::

::: {.doxygenfunction project="progguide"}
split_words
:::

::: {.doxygenfunction project="progguide"}
split_lines
:::

::: {.doxygenfunction project="progguide"}
strmatch
:::

::: {.doxygenfunction project="progguide"}
strfind
:::

::: {.doxygenfunction project="progguide"}
is_integer
:::

::: {.doxygenfunction project="progguide"}
is_double
:::

::: {.doxygenfunction project="progguide"}
is_id
:::

::: {.doxygenfunction project="progguide"}
is_type
:::

## Potential file functions

::: {.doxygenfunction project="progguide"}
get_potential_file_path
:::

::: {.doxygenfunction project="progguide"}
get_potential_date
:::

::: {.doxygenfunction project="progguide"}
get_potential_units
:::

::: {.doxygenfunction project="progguide"}
get_supported_conversions
:::

::: {.doxygenfunction project="progguide"}
get_conversion_factor
:::

::: {.doxygenfunction project="progguide"}
open_potential(const std::string &name, LAMMPS *lmp, int*auto_convert)
:::

## Argument processing

::: {.doxygenfunction project="progguide"}
bounds
:::

::: {.doxygenfunction project="progguide"}
expand_args
:::

::: {.doxygenfunction project="progguide"}
parse_grid_id
:::

::: {.doxygenfunction project="progguide"}
expand_type
:::

## Convenience functions

::: {.doxygenfunction project="progguide"}
logmesg(LAMMPS \*lmp, const std::string &format, Args&&\... args)
:::

::: {.doxygenfunction project="progguide"}
logmesg(LAMMPS \*lmp, const std::string &mesg)
:::

::: {.doxygenfunction project="progguide"}
errorurl
:::

::: {.doxygenfunction project="progguide"}
missing_cmd_args
:::

::: {.doxygenfunction project="progguide"}
flush_buffers(LAMMPS \*lmp)
:::

::: {.doxygenfunction project="progguide"}
getsyserror
:::

::: {.doxygenfunction project="progguide"}
check_packages_for_style
:::

::: {.doxygenfunction project="progguide"}
timespec2seconds
:::

::: {.doxygenfunction project="progguide"}
date2num
:::

::: {.doxygenfunction project="progguide"}
current_date
:::

## Customized standard functions

::: {.doxygenfunction project="progguide"}
binary_search
:::

::: {.doxygenfunction project="progguide"}
merge_sort
:::

------------------------------------------------------------------------

# Special Math functions

The `MathSpecial` namespace implements a selection of custom and
optimized mathematical functions for a variety of applications.

::: {.doxygenfunction project="progguide"}
factorial
:::

::: {.doxygenfunction project="progguide"}
exp2_x86
:::

::: {.doxygenfunction project="progguide"}
fm_exp
:::

::: {.doxygenfunction project="progguide"}
my_erfcx
:::

::: {.doxygenfunction project="progguide"}
expmsq
:::

::: {.doxygenfunction project="progguide"}
square
:::

::: {.doxygenfunction project="progguide"}
cube
:::

::: {.doxygenfunction project="progguide"}
powsign
:::

::: {.doxygenfunction project="progguide"}
powint
:::

::: {.doxygenfunction project="progguide"}
powsinxx
:::

------------------------------------------------------------------------

# Tokenizer classes

The purpose of the tokenizer classes is to simplify the recurring task
of breaking lines of text down into words and/or numbers. Traditionally,
LAMMPS code would be using the `strtok()` function from the C library
for that purpose, but that function has two significant
disadvantages: 1) it cannot be used concurrently from different LAMMPS
instances since it stores its status in a global variable and 2) it
modifies the string that it is processing. These classes were
implemented to avoid both of these issues and also to reduce the amount
of code that needs to be written.

The basic procedure is to create an instance of the tokenizer class with
the string to be processed as an argument and then do a loop until all
available tokens are read. The constructor has a default set of
separator characters, but that can be overridden. The default separators
are all \"whitespace\" characters, i.e. the space character, the
tabulator character, the carriage return character, the linefeed
character, and the form feed character.

``` {.c++ caption="Tokenizer class example listing entries of the PATH environment variable"}
#include "tokenizer.h"
#include <cstdlib>
#include <string>
#include <iostream>

using namespace LAMMPS_NS;

int main(int, char **)
{
    const char *path = getenv("PATH");

    if (path != nullptr) {
        Tokenizer p(path,":");
        while (p.has_next())
            std::cout << "Entry: " << p.next() << "\n";
    }
    return 0;
}
```

Most tokenizer operations cannot fail except for
`LAMMPS_NS::Tokenizer::next`{.interpreted-text role="cpp:func"} (when
used without first checking with
`LAMMPS_NS::Tokenizer::has_next`{.interpreted-text role="cpp:func"}) and
`LAMMPS_NS::Tokenizer::skip`{.interpreted-text role="cpp:func"}. In case
of failure, the class will throw an exception, so you may need to wrap
the code using the tokenizer into a `try` / `catch` block to handle
errors. The `LAMMPS_NS::ValueTokenizer`{.interpreted-text
role="cpp:class"} class may also throw an exception when a (type of)
number is requested as next token that is not compatible with the string
representing the next word.

``` {.c++ caption="ValueTokenizer class example with exception handling"}
#include "tokenizer.h"
#include <cstdlib>
#include <string>
#include <iostream>

using namespace LAMMPS_NS;

int main(int, char **)
{
    const char *text = "1 2 3 4 5 20.0 21 twentytwo 2.3";
    double num1(0),num2(0),num3(0),num4(0);

    ValueTokenizer t(text);
    // read 4 doubles after skipping over 5 numbers
    try {
        t.skip(5);
        num1 = t.next_double();
        num2 = t.next_double();
        num3 = t.next_double();
        num4 = t.next_double();
    } catch (TokenizerException &e) {
        std::cout << "Reading numbers failed: " << e.what() << "\n";
    }
    std::cout << "Values: " << num1 << " " << num2 << " " << num3 << " " << num4 << "\n";
    return 0;
}
```

This code example should produce the following output:

``` 
Reading numbers failed: Not a valid floating-point number: 'twentytwo'
Values: 20 21 0 0
```

------------------------------------------------------------------------

::: {.doxygenclass project="progguide" members=""}
LAMMPS_NS::Tokenizer
:::

::: {.doxygenclass project="progguide" members=""}
LAMMPS_NS::TokenizerException
:::

::: {.doxygenclass project="progguide" members=""}
LAMMPS_NS::ValueTokenizer
:::

::: {.doxygenclass project="progguide" members=""}
LAMMPS_NS::InvalidIntegerException
:::

::: {.doxygenclass project="progguide" members=""}
LAMMPS_NS::InvalidFloatException
:::

------------------------------------------------------------------------

# Argument parsing classes

The purpose of argument parsing classes it to simplify and unify how
arguments of commands in LAMMPS are parsed and to make abstractions of
repetitive tasks.

The `LAMMPS_NS::ArgInfo`{.interpreted-text role="cpp:class"} class
provides an abstraction for parsing references to compute or fix styles,
variables or custom integer or double properties handled by [fix
property/atom](fix_property_atom). These would start with a \"c\_\",
\"f\_\", \"v\_\", \"d\_\", \"d2\_\", \"i\_\", or \"i2\_\" followed by
the ID or name of than instance and may be postfixed with one or two
array indices \"\[\<number\>\]\" with numbers \> 0.

A typical code segment would look like this:

``` {.c++ caption="Usage example for ArgInfo class"}
int nvalues = 0;
for (iarg = 0; iarg < nargnew; iarg++) {
  ArgInfo argi(arg[iarg]);

  which[nvalues] = argi.get_type();
  argindex[nvalues] = argi.get_index1();
  ids[nvalues] = argi.copy_name();

  if ((which[nvalues] == ArgInfo::UNKNOWN)
       || (which[nvalues] == ArgInfo::NONE)
       || (argi.get_dim() > 1))
    error->all(FLERR,"Illegal compute XXX command");

  nvalues++;
}
```

------------------------------------------------------------------------

::: {.doxygenclass project="progguide" members=""}
LAMMPS_NS::ArgInfo
:::

------------------------------------------------------------------------

# File reader classes

The purpose of the file reader classes is to simplify the recurring task
of reading and parsing files. They can use the
`ValueTokenizer <LAMMPS_NS::ValueTokenizer>`{.interpreted-text
role="cpp:class"}\_\_ class to process the read in text. The
`TextFileReader <LAMMPS_NS::TextFileReader>`{.interpreted-text
role="cpp:class"}\_\_ is a more general version while
`PotentialFileReader <LAMMPS_NS::PotentialFileReader>`{.interpreted-text
role="cpp:class"}\_\_ is specialized to implement the behavior expected
for looking up and reading/parsing files with potential parameters in
LAMMPS. The potential file reader class requires a LAMMPS instance,
requires to be run on MPI rank 0 only, will use the
`utils::get_potential_file_path <LAMMPS_NS::utils::get_potential_file_path>`{.interpreted-text
role="cpp:func"}\_\_ function to look up and open the file, and will
call the `LAMMPS_NS::Error`{.interpreted-text role="cpp:class"} class in
case of failures to read or to convert numbers, so that LAMMPS will be
aborted.

``` {.c++ caption="Use of PotentialFileReader class in pair style coul/streitz"}
PotentialFileReader reader(lmp, file, "coul/streitz");
char * line;

while((line = reader.next_line(NPARAMS_PER_LINE))) {
  try {
    ValueTokenizer values(line);
    std::string iname = values.next_string();

    int ielement;
    for (ielement = 0; ielement < nelements; ielement++)
      if (iname == elements[ielement]) break;

    if (nparams == maxparam) {
      maxparam += DELTA;
      params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                          "pair:params");
    }

    params[nparams].ielement = ielement;
    params[nparams].chi = values.next_double();
    params[nparams].eta = values.next_double();
    params[nparams].gamma = values.next_double();
    params[nparams].zeta = values.next_double();
    params[nparams].zcore = values.next_double();

  } catch (TokenizerException & e) {
    error->one(FLERR, e.what());
  }
  nparams++;
}
```

A file that would be parsed by the reader code fragment looks like this:

    # DATE: 2015-02-19 UNITS: metal CONTRIBUTOR: Ray Shan CITATION: Streitz and Mintmire, Phys Rev B, 50, 11996-12003 (1994)
    #
    # X (eV)                J (eV)          gamma (1/\AA)   zeta (1/\AA)    Z (e)

    Al      0.000000        10.328655       0.000000        0.968438        0.763905
    O       5.484763        14.035715       0.000000        2.143957        0.000000

------------------------------------------------------------------------

::: {.doxygenclass project="progguide" members=""}
LAMMPS_NS::TextFileReader
:::

::: {.doxygenclass project="progguide" members=""}
LAMMPS_NS::PotentialFileReader
:::

------------------------------------------------------------------------

# Memory pool classes

The memory pool classes are used for cases where otherwise many small
memory allocations would be needed and where the data would be either
all used or all freed. One example for that is the storage of neighbor
lists. The memory management strategy is based on the assumption that
allocations will be in chunks of similar sizes. The allocation is then
not done per individual call for a reserved chunk of memory, but for a
\"page\" that can hold multiple chunks of data. A parameter for the
maximum chunk size must be provided, as that is used to determine
whether a new page of memory must be used.

The `MyPage <LAMMPS_NS::MyPage>`{.interpreted-text role="cpp:class"}\_\_
class offers two ways to reserve a chunk: 1) with
`MyPage::get() <LAMMPS_NS::MyPage::get>`{.interpreted-text
role="cpp:func"}\_\_ the chunk size needs to be known in advance, 2)
with `MyPage::vget() <LAMMPS_NS::MyPage::vget>`{.interpreted-text
role="cpp:func"}\_\_ a pointer to the next chunk is returned, but its
size is registered later with
`MyPage::vgot() <LAMMPS_NS::MyPage::vgot>`{.interpreted-text
role="cpp:func"}\_\_.

``` {.c++ caption="Example of using :cpp:class:`MyPage <LAMMPS_NS::MyPage>`__"}
#include "my_page.h"
using namespace LAMMPS_NS;

MyPage<double> *dpage = new MyPage<double>;
// max size of chunk: 256, size of page: 10240 doubles (=81920 bytes)
dpage->init(256,10240);

double **build_some_lists(int num)
{
    dpage->reset();
    double **dlist = new double*[num];
    for (int i=0; i < num; ++i) {
        double *dptr = dpage.vget();
        int jnum = 0;
        for (int j=0; j < jmax; ++j) {
            // compute some dvalue for eligible loop index j
            dptr[j] = dvalue;
            ++jnum;
        }
        if (dpage.status() != 0) {
            // handle out of memory or jnum too large errors
        }
        dpage.vgot(jnum);
        dlist[i] = dptr;
    }
    return dlist;
}
```

------------------------------------------------------------------------

::: {.doxygenclass project="progguide" members=""}
LAMMPS_NS::MyPage
:::

::: {.doxygenclass project="progguide" members=""}
LAMMPS_NS::MyPoolChunk
:::

------------------------------------------------------------------------

# Eigensolver functions

The `MathEigen` sub-namespace of the `LAMMPS_NS` namespace contains
functions and classes for eigensolvers. Currently only the
`jacobi3 function <MathEigen::jacobi3>`{.interpreted-text
role="cpp:func"}\_\_ is used in various places in LAMMPS. That function
is built on top of a group of more generic eigensolvers that are
maintained in the `math_eigen_impl.h` header file. This header contains
the implementation of three template classes:

1.  \"Jacobi\" calculates all of the eigenvalues and eigenvectors of a
    dense, symmetric, real matrix.
2.  The \"PEigenDense\" class only calculates the principal eigenvalue
    (ie. the largest or smallest eigenvalue), and its corresponding
    eigenvector. However it is much more efficient than \"Jacobi\" when
    applied to large matrices (larger than 13x13). PEigenDense also can
    understand complex-valued Hermitian matrices.
3.  The \"LambdaLanczos\" class is a generalization of \"PEigenDense\"
    which can be applied to arbitrary sparse matrices.

The \"math_eigen_impl.h\" code is an amalgamation of
[jacobi_pd](https://github.com/jewettaij/jacobi_pd)\_ by Andrew Jewett
at Scripps Research (under CC0-1.0 license) and [Lambda
Lanczos](https://github.com/mrcdr/lambda-lanczos)\_ by Yuya Kurebayashi
at Tohoku University (under MIT license)

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
MathEigen::jacobi3(double const *const*mat, double *eval, double*\*evec)
:::

::: {.doxygenfunction project="progguide"}
MathEigen::jacobi3(double const mat\[3\]\[3\], double \*eval, double
evec\[3\]\[3\])
:::

------------------------------------------------------------------------

# Communication buffer coding with *ubuf* {#communication_buffer_coding_with_ubuf}

LAMMPS uses communication buffers where it collects data from various
class instances and then exchanges the data with neighboring subdomains.
For simplicity those buffers are defined as `double` buffers and used
for doubles and integer numbers. This presents a unique problem when
64-bit integers are used. While the storage needed for a `double` is
also 64-bit, it cannot be used by a simple assignment. To get around
that limitation, LAMMPS uses the
`ubuf <LAMMPS_NS::ubuf>`{.interpreted-text role="cpp:union"}\_\_ union.
It is used in the various \"pack\" and \"unpack\" functions in the
LAMMPS classes to store and retrieve integers that may be 64-bit from
the communication buffers.

------------------------------------------------------------------------

::: {.doxygenunion project="progguide"}
LAMMPS_NS::ubuf
:::
