# dihedral_style spherical command

## Syntax

``` LAMMPS
dihedral_style spherical
```

## Examples

``` LAMMPS
dihedral_coeff 1 1  286.1  1 124  1    1 90.0 0    1 90.0 0
dihedral_coeff 1 3  69.3   1 93.9 1    1 90   0    1 90   0  &
                    49.1   0 0.00 0    1 74.4 1    0 0.00 0  &
                    25.2   0 0.00 0    0 0.00 0    1 48.1 1
```

## Description

The *spherical* dihedral style uses the potential:

![image](JPG/dihedral_spherical_angles.jpg){.align-center}

$$\begin{aligned}
E(\phi,\theta_1,\theta_2) & = \sum_{i=1}^N\nolimits\ C_i\ \Phi_i(\phi)\ \Theta_{1i}(\theta_1)\ \Theta_{2i}(\theta_2) \\
\Phi_{i}(\phi)            & = u_i - \mathrm{cos}((\phi   - a_i)K_i) \\
\Theta_{1i}(\theta_1)     & = v_i - \mathrm{cos}((\theta_1-b_i)L_i) \\
\Theta_{2i}(\theta_2)     & = w_i - \mathrm{cos}((\theta_2-c_i)M_i)
\end{aligned}$$

For this dihedral style, the energy can be any function that combines
the 4-body dihedral-angle ($\phi$) and the two 3-body bond-angles
($\theta_1$, $\theta_2$). For this reason, there is usually no need to
define 3-body \"angle\" forces separately for the atoms participating in
these interactions. It is probably more efficient to incorporate 3-body
angle forces into the dihedral interaction even if it requires adding
additional terms to the expansion (as was done in the second example). A
careful choice of parameters can prevent singularities that occur with
traditional force-fields whenever theta1 or theta2 approach 0 or 180
degrees.

The last example above corresponds to an interaction with a single
energy minima located near $\phi=93.9$, $\theta_1=74.4$, $\theta_2=48.1$
degrees, and it remains numerically stable at all angles ($\phi$,
$\theta_1$, $\theta_2$). In this example, the coefficients 49.1, and
25.2 can be physically interpreted as the harmonic spring constants for
theta1 and theta2 around their minima. The coefficient 69.3 is the
harmonic spring constant for phi after division by sin(74.4)\*sin(48.1)
(the minima positions for theta1 and theta2).

The following coefficients must be defined for each dihedral type via
the [dihedral_coeff](dihedral_coeff) command as in the example above, or
in the Dihedral Coeffs section of a data file read by the
[read_data](read_data) command:

-   $n$ (integer \>= 1)
-   $C_1$ (energy)
-   $K_1$ (typically an integer)
-   $a_1$ (degrees)
-   $u_1$ (typically 0.0 or 1.0)
-   $L_1$ (typically an integer)
-   $b_1$ (degrees, typically 0.0 or 90.0)
-   $v_1$ (typically 0.0 or 1.0)
-   $M_1$ (typically an integer)
-   $c_1$ (degrees, typically 0.0 or 90.0)
-   $w_1$ (typically 0.0 or 1.0)
-   \[\...\]
-   $C_n$ (energy)
-   $K_n$ (typically an integer)
-   $a_n$ (degrees)
-   $u_n$ (typically 0.0 or 1.0)
-   $L_n$ (typically an integer)
-   $b_n$ (degrees, typically 0.0 or 90.0)
-   $v_n$ (typically 0.0 or 1.0)
-   $M_n$ (typically an integer)
-   $c_n$ (degrees, typically 0.0 or 90.0)
-   $w_n$ (typically 0.0 or 1.0)

------------------------------------------------------------------------

## Restrictions

This dihedral style can only be used if LAMMPS was built with the
EXTRA-MOLECULE package. See the [Build package](Build_package) doc page
for more info.

## Related commands

[dihedral_coeff](dihedral_coeff)

## Default

none
