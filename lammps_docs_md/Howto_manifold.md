# Manifolds (surfaces)

**Overview:**

This page is not about a LAMMPS input script command, but about
manifolds, which are generalized surfaces, as defined and used by the
MANIFOLD package, to track particle motion on the manifolds. See the
src/MANIFOLD/README file for more details about the package and its
commands.

Below is a list of currently supported manifolds by the MANIFOLD
package, their parameters and a short description of them. The
parameters listed here are in the same order as they should be passed to
the relevant fixes.

  --------------- -------------- ------------------------------- ------------------------------
  *manifold*      *parameters*   *equation*                      *description*

  cylinder        R              x\^2 + y\^2 - R\^2 = 0          Cylinder along z-axis, axis
                                                                 going through (0,0,0)

  cylinder_dent   R l a          x\^2 + y\^2 - r(z)\^2 = 0, r(x) A cylinder with a dent around
                                 = R if \| z \| \> l, r(z) = R - z = 0
                                 a\*(1 + cos(z/l))/2 otherwise   

  dumbbell        a A B c        -( x\^2 + y\^2 ) + (a\^2 -      A dumbbell
                                 z\^2/c\^2) \* ( 1 +             
                                 (A\*sin(B\*z\^2))\^4) = 0       

  ellipsoid       a b c          (x/a)\^2 + (y/b)\^2 + (z/c)\^2  An ellipsoid
                                 = 0                             

  gaussian_bump   A l rc1 rc2    if( x \< rc1) -z + A \* exp(    A Gaussian bump at x = y = 0,
                                 -x\^2 / (2 l\^2) ); else if( x  smoothly tapered to a flat
                                 \< rc2 ) -z + a + b\*x +        plane z = 0.
                                 c\*x\^2 + d\*x\^3; else z       

  plane           a b c x0 y0 z0 a\*(x-x0) + b\*(y-y0) +         A plane with normal (a,b,c)
                                 c\*(z-z0) = 0                   going through point (x0,y0,z0)

  plane_wiggle    a w            z - a\*sin(w\*x) = 0            A plane with a sinusoidal
                                                                 modulation on z along x.

  sphere          R              x\^2 + y\^2 + z\^2 - R\^2 = 0   A sphere of radius R

  supersphere     R q            \| x \|\^q + \| y \|\^q + \| z  A supersphere of hyperradius R
                                 \|\^q - R\^q = 0                

  spine           a, A, B, B2, c -(x\^2 + y\^2) + (a\^2 -        An approximation to a
                                 z\^2/f(z)\^2)\*(1 +             dendritic spine
                                 (A\*sin(g(z)\*z\^2))\^4), f(z)  
                                 = c if z \> 0, 1 otherwise;     
                                 g(z) = B if z \> 0, B2          
                                 otherwise                       

  spine_two       a, A, B, B2, c -(x\^2 + y\^2) + (a\^2 -        Another approximation to a
                                 z\^2/f(z)\^2)\*(1 +             dendritic spine
                                 (A\*sin(g(z)\*z\^2))\^2), f(z)  
                                 = c if z \> 0, 1 otherwise;     
                                 g(z) = B if z \> 0, B2          
                                 otherwise                       

  thylakoid       wB LB lB       Various, see                    
                                 [(Paquay)](Paquay1) \| A model  
                                 grana thylakoid consisting of   
                                 two block-like compartments     
                                 connected by a bridge of width  
                                 wB, length LB and taper length  
                                 lB \|                           

  torus           R r            (R - sqrt( x\^2 + y\^2 ) )\^2 + A torus with large radius R
                                 z\^2 - r\^2                     and small radius r, centered
                                                                 on (0,0,0)
  --------------- -------------- ------------------------------- ------------------------------

::: {#Paquay1}
**(Paquay)** Paquay and Kusters, Biophys. J., 110, 6, (2016). preprint
available at [arXiv:1411.3019](https://arxiv.org/abs/1411.3019/)\_.
:::
