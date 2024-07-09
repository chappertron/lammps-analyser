# Compute, fixes, variables

This section documents accessing or modifying data from objects like
computes, fixes, or variables in LAMMPS using the
`lammps`{.interpreted-text role="py:mod"} module.

::::: tabs
::: tab
lammps API

For
`lammps.extract_compute() <lammps.lammps.extract_compute()>`{.interpreted-text
role="py:meth"}\_\_ and
`lammps.extract_fix() <lammps.lammps.extract_fix()>`{.interpreted-text
role="py:meth"}\_\_, the global, per-atom, or local data calculated by
the compute or fix can be accessed. What is returned depends on whether
the compute or fix calculates a scalar or vector or array. For a scalar,
a single double value is returned. If the compute or fix calculates a
vector or array, a pointer to the internal LAMMPS data is returned,
which you can use via normal Python subscripting.

The one exception is that for a fix that calculates a global vector or
array, a single double value from the vector or array is returned,
indexed by I (vector) or I and J (array). I,J are zero-based indices.
See the [Howto output](Howto_output) page for a discussion of global,
per-atom, and local data, and of scalar, vector, and array data types.
See the doc pages for individual [computes](compute) and [fixes](fix)
for a description of what they calculate and store.

For
`lammps.extract_variable() <lammps.lammps.extract_variable()>`{.interpreted-text
role="py:meth"}\_\_, an [equal-style or atom-style variable](variable)
is evaluated and its result returned.

For equal-style variables a single `c_double` value is returned and the
group argument is ignored. For atom-style variables, a vector of
`c_double` is returned, one value per atom, which you can use via normal
Python subscripting. The values will be zero for atoms not in the
specified group.

`lammps.numpy.extract_compute() <lammps.numpy_wrapper.numpy_wrapper.extract_compute()>`{.interpreted-text
role="py:meth"}\_\_,
`lammps.numpy.extract_fix() <lammps.numpy_wrapper.numpy_wrapper.extract_fix()>`{.interpreted-text
role="py:meth"}\_\_, and
`lammps.numpy.extract_variable() <lammps.numpy_wrapper.numpy_wrapper.extract_variable()>`{.interpreted-text
role="py:meth"}\_\_ are equivalent NumPy implementations that return
NumPy arrays instead of `ctypes` pointers.

The
`lammps.set_variable() <lammps.lammps.set_variable()>`{.interpreted-text
role="py:meth"}\_\_ method sets an existing string-style variable to a
new string value, so that subsequent LAMMPS commands can access the
variable.

**Methods**:

-   `lammps.extract_compute() <lammps.lammps.extract_compute()>`{.interpreted-text
    role="py:meth"}\_\_: extract value(s) from a compute
-   `lammps.extract_fix() <lammps.lammps.extract_fix()>`{.interpreted-text
    role="py:meth"}\_\_: extract value(s) from a fix
-   `lammps.extract_variable() <lammps.lammps.extract_variable()>`{.interpreted-text
    role="py:meth"}\_\_: extract value(s) from a variable
-   `lammps.set_variable() <lammps.lammps.set_variable()>`{.interpreted-text
    role="py:meth"}\_\_: set existing named string-style variable to
    value

**NumPy Methods**:

-   `lammps.numpy.extract_compute() <lammps.numpy_wrapper.numpy_wrapper.extract_compute()>`{.interpreted-text
    role="py:meth"}\_\_: extract value(s) from a compute, return arrays
    as numpy arrays
-   `lammps.numpy.extract_fix() <lammps.numpy_wrapper.numpy_wrapper.extract_fix()>`{.interpreted-text
    role="py:meth"}\_\_: extract value(s) from a fix, return arrays as
    numpy arrays
-   `lammps.numpy.extract_variable() <lammps.numpy_wrapper.numpy_wrapper.extract_variable()>`{.interpreted-text
    role="py:meth"}\_\_: extract value(s) from a variable, return arrays
    as numpy arrays
:::

::: tab
PyLammps/IPyLammps API

PyLammps and IPyLammps classes currently do not add any additional ways
of retrieving information out of computes and fixes. This information
can still be accessed by using the lammps API:

``` python
L.lmp.extract_compute(...)
L.lmp.extract_fix(...)
# OR
L.lmp.numpy.extract_compute(...)
L.lmp.numpy.extract_fix(...)
```

LAMMPS variables can be both defined and accessed via the
`PyLammps <lammps.PyLammps>`{.interpreted-text role="py:class"}\_\_
interface.

To define a variable you can use the [variable](variable) command:

``` python
L.variable("a index 2")
```

A dictionary of all variables is returned by the
`PyLammps.variables <lammps.PyLammps.variables>`{.interpreted-text
role="py:attr"}\_\_ property:

you can access an individual variable by retrieving a variable object
from the `L.variables` dictionary by name

``` python
a = L.variables['a']
```

The variable value can then be easily read and written by accessing the
value property of this object.

``` python
print(a.value)
a.value = 4
```
:::
:::::
