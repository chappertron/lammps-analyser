# pair_style granular command

## Syntax

``` LAMMPS
pair_style granular cutoff
```

-   cutoff = global cutoff (optional). See discussion below.

## Examples

``` LAMMPS
pair_style granular
pair_coeff * * hooke 1000.0 50.0 tangential linear_nohistory 1.0 0.4 damping mass_velocity

pair_style granular
pair_coeff * * hooke 1000.0 50.0 tangential linear_history 500.0 1.0 0.4 damping mass_velocity

pair_style granular
pair_coeff * * hertz 1000.0 50.0 tangential mindlin 1000.0 1.0 0.4 limit_damping

pair_style granular
pair_coeff * * hertz/material 1e8 0.3 0.3 tangential mindlin_rescale NULL 1.0 0.4 damping tsuji

pair_style granular
pair_coeff 1 * jkr 1000.0 500.0 0.3 10 tangential mindlin 800.0 1.0 0.5 rolling sds 500.0 200.0 0.5 twisting marshall
pair_coeff 2 2 hertz 200.0 100.0 tangential linear_history 300.0 1.0 0.1 rolling sds 200.0 100.0 0.1 twisting marshall

pair_style granular
pair_coeff 1 1 dmt 1000.0 50.0 0.3 0.0 tangential mindlin NULL 0.5 0.5 rolling sds 500.0 200.0 0.5 twisting marshall
pair_coeff 2 2 dmt 1000.0 50.0 0.3 10.0 tangential mindlin NULL 0.5 0.1 rolling sds 500.0 200.0 0.1 twisting marshall

pair_style granular
pair_coeff * * hertz 1000.0 50.0 tangential mindlin 1000.0 1.0 0.4 heat area 0.1
```

## Description

The *granular* styles support a variety of options for the normal,
tangential, rolling and twisting forces resulting from contact between
two granular particles. This expands on the options offered by the [pair
gran/\*](pair_gran) pair styles. The total computed forces and torques
are the sum of various models selected for the normal, tangential,
rolling and twisting modes of motion.

All model choices and parameters are entered in the
[pair_coeff](pair_coeff) command, as described below. Unlike e.g. [pair
gran/hooke](pair_gran), coefficient values are not global, but can be
set to different values for different combinations of particle types, as
determined by the [pair_coeff](pair_coeff) command. If the contact model
choice is the same for two particle types, the mixing for the
cross-coefficients can be carried out automatically. This is shown in
the last example, where model choices are the same for type 1 - type 1
as for type 2 - type2 interactions, but coefficients are different. In
this case, the mixed coefficients for type 1 - type 2 interactions can
be determined from mixing rules discussed below. For additional
flexibility, coefficients as well as model forms can vary between
particle types, as shown in the fourth example: type 1 - type 1
interactions are based on a Johnson-Kendall-Roberts normal contact model
and 2-2 interactions are based on a DMT cohesive model (see below). In
that example, 1-1 and 2-2 interactions have different model forms, in
which case mixing of coefficients cannot be determined, so 1-2
interactions must be explicitly defined via the *pair_coeff 1 \**
command, otherwise an error would result.

------------------------------------------------------------------------

The first required keyword for the *pair_coeff* command is the normal
contact model. Currently supported options for normal contact models and
their required arguments are:

1.  *hooke* : $k_n$, $\eta_{n0}$ (or $e$)
2.  *hertz* : $k_n$, $\eta_{n0}$ (or $e$)
3.  *hertz/material* : E, $\eta_{n0}$ (or $e$), $\nu$
4.  *dmt* : E, $\eta_{n0}$ (or $e$), $\nu$, $\gamma$
5.  *jkr* : E, $\eta_{n0}$ (or $e$), $\nu$, $\gamma$

Here, $k_n$ is spring stiffness (with units that depend on model choice,
see below); $\eta_{n0}$ is a damping prefactor (or, in its place a
coefficient of restitution $e$, depending on the choice of damping mode,
see below); E is Young\'s modulus in units of *force*/\*length\*\^2,
i.e. *pressure*; $\nu$ is Poisson\'s ratio and $\gamma$ is a surface
energy density, in units of *energy*/\*length\*\^2.

For the *hooke* model, the normal, elastic component of force acting on
particle *i* due to contact with particle *j* is given by:

$$\mathbf{F}_{ne, Hooke} = k_n \delta_{ij} \mathbf{n}$$

Where $\delta_{ij} = R_i + R_j - \|\mathbf{r}_{ij}\|$ is the particle
overlap, $R_i, R_j$ are the particle radii,
$\mathbf{r}_{ij} = \mathbf{r}_i - \mathbf{r}_j$ is the vector separating
the two particle centers (note the i-j ordering so that
$\mathbf{F}_{ne}$ is positive for repulsion), and
$\mathbf{n} = \frac{\mathbf{r}_{ij}}{\|\mathbf{r}_{ij}\|}$. Therefore,
for *hooke*, the units of the spring constant $k_n$ are
*force*/\*distance\*, or equivalently *mass*/\*time\^2\*.

For the *hertz* model, the normal component of force is given by:

$$\mathbf{F}_{ne, Hertz} = k_n R_{eff}^{1/2}\delta_{ij}^{3/2} \mathbf{n}$$

Here, $R_{eff} = \frac{R_i R_j}{R_i + R_j}$ is the effective radius,
denoted for simplicity as *R* from here on. For *hertz*, the units of
the spring constant $k_n$ are *force*/\*length\*\^2, or equivalently
*pressure*.

For the *hertz/material* model, the force is given by:

$$\mathbf{F}_{ne, Hertz/material} = \frac{4}{3} E_{eff} R_{eff}^{1/2}\delta_{ij}^{3/2} \mathbf{n}$$

Here,
$E_{eff} = E = \left(\frac{1-\nu_i^2}{E_i} + \frac{1-\nu_j^2}{E_j}\right)^{-1}$
is the effective Young\'s modulus, with $\nu_i, \nu_j$ the Poisson
ratios of the particles of types *i* and *j*. Note that if the elastic
modulus and the shear modulus of the two particles are the same, the
*hertz/material* model is equivalent to the *hertz* model with
$k_n = 4/3 E_{eff}$

The *dmt* model corresponds to the [(Derjaguin-Muller-Toporov)](DMT1975)
cohesive model, where the force is simply Hertz with an additional
attractive cohesion term:

$$\mathbf{F}_{ne, dmt} = \left(\frac{4}{3} E R^{1/2}\delta_{ij}^{3/2} - 4\pi\gamma R\right)\mathbf{n}$$

The *jkr* model is the [(Johnson-Kendall-Roberts)](JKR1971) model, where
the force is computed as:

$$\mathbf{F}_{ne, jkr} = \left(\frac{4Ea^3}{3R} - 2\pi a^2\sqrt{\frac{4\gamma E}{\pi a}}\right)\mathbf{n}$$

Here, $a$ is the radius of the contact zone, related to the overlap
$\delta$ according to:

$$\delta = a^2/R - 2\sqrt{\pi \gamma a/E}$$

LAMMPS internally inverts the equation above to solve for *a* in terms
of $\delta$, then solves for the force in the previous equation.
Additionally, note that the JKR model allows for a tensile force beyond
contact (i.e. for $\delta < 0$), up to a maximum of $3\pi\gamma R$ (also
known as the \'pull-off\' force). Note that this is a hysteretic effect,
where particles that are not contacting initially will not experience
force until they come into contact $\delta \geq 0$; as they move apart
and ($\delta < 0$), they experience a tensile force up to
$3\pi\gamma R$, at which point they lose contact.

------------------------------------------------------------------------

In addition, the normal force is augmented by a damping term of the
following general form:

$$\mathbf{F}_{n,damp} = -\eta_n \mathbf{v}_{n,rel}$$

Here,
$\mathbf{v}_{n,rel} = (\mathbf{v}_j - \mathbf{v}_i) \cdot \mathbf{n}\ \mathbf{n}$
is the component of relative velocity along $\mathbf{n}$.

The optional *damping* keyword to the *pair_coeff* command followed by a
keyword determines the model form of the damping factor $\eta_n$, and
the interpretation of the $\eta_{n0}$ or $e$ coefficients specified as
part of the normal contact model settings. The *damping* keyword and
corresponding model form selection may be appended anywhere in the *pair
coeff* command. Note that the choice of damping model affects both the
normal and tangential damping (and depending on other settings,
potentially also the twisting damping). The options for the damping
model currently supported are:

1.  *velocity*
2.  *mass_velocity*
3.  *viscoelastic*
4.  *tsuji*

If the *damping* keyword is not specified, the *viscoelastic* model is
used by default.

For *damping velocity*, the normal damping is simply equal to the
user-specified damping coefficient in the *normal* model:

$$\eta_n = \eta_{n0}$$

Here, $\eta_{n0}$ is the damping coefficient specified for the normal
contact model, in units of *mass*/\*time\*.

For *damping mass_velocity*, the normal damping is given by:

$$\eta_n = \eta_{n0} m_{eff}$$

Here, $\eta_{n0}$ is the damping coefficient specified for the normal
contact model, in units of 1/\*time\* and
$m_{eff} = m_i m_j/(m_i + m_j)$ is the effective mass. Use *damping
mass_velocity* to reproduce the damping behavior of *pair
gran/hooke/\**.

The *damping viscoelastic* model is based on the viscoelastic treatment
of [(Brilliantov et al)](Brill1996), where the normal damping is given
by:

$$\eta_n = \eta_{n0}\ a m_{eff}$$

Here, *a* is the contact radius, given by $a =\sqrt{R\delta}$ for all
models except *jkr*, for which it is given implicitly according to
$\delta = a^2/R - 2\sqrt{\pi \gamma a/E}$. For *damping viscoelastic*,
$\eta_{n0}$ is in units of 1/(\*time\*\*\*distance\*).

The *tsuji* model is based on the work of [(Tsuji et al)](Tsuji1992).
Here, the damping coefficient specified as part of the normal model is
interpreted as a restitution coefficient $e$. The damping constant
$\eta_n$ is given by:

$$\eta_n = \alpha (m_{eff}k_n)^{1/2}$$

For normal contact models based on material parameters, $k_n = 4/3Ea$.
This damping model is not compatible with cohesive normal models such as
*JKR* or *DMT*. The parameter $\alpha$ is related to the restitution
coefficient *e* according to:

$$\alpha = 1.2728-4.2783e+11.087e^2-22.348e^3+27.467e^4-18.022e^5+4.8218e^6$$

The dimensionless coefficient of restitution $e$ specified as part of
the normal contact model parameters should be between 0 and 1, but no
error check is performed on this.

The total normal force is computed as the sum of the elastic and damping
components:

$$\mathbf{F}_n = \mathbf{F}_{ne} + \mathbf{F}_{n,damp}$$

------------------------------------------------------------------------

The *pair_coeff* command also requires specification of the tangential
contact model. The required keyword *tangential* is expected, followed
by the model choice and associated parameters. Currently supported
tangential model choices and their expected parameters are as follows:

1.  *linear_nohistory* : $x_{\gamma,t}$, $\mu_s$
2.  *linear_history* : $k_t$, $x_{\gamma,t}$, $\mu_s$
3.  *mindlin* : $k_t$ or NULL, $x_{\gamma,t}$, $\mu_s$
4.  *mindlin/force* : $k_t$ or NULL, $x_{\gamma,t}$, $\mu_s$
5.  *mindlin_rescale* : $k_t$ or NULL, $x_{\gamma,t}$, $\mu_s$
6.  *mindlin_rescale/force* : $k_t$ or NULL, $x_{\gamma,t}$, $\mu_s$

Here, $x_{\gamma,t}$ is a dimensionless multiplier for the normal
damping $\eta_n$ that determines the magnitude of the tangential
damping, $\mu_t$ is the tangential (or sliding) friction coefficient,
and $k_t$ is the tangential stiffness coefficient.

For *tangential linear_nohistory*, a simple velocity-dependent Coulomb
friction criterion is used, which mimics the behavior of the *pair
gran/hooke* style. The tangential force $\mathbf{F}_t$ is given by:

$$\mathbf{F}_t =  -\min(\mu_t F_{n0}, \|\mathbf{F}_\mathrm{t,damp}\|) \mathbf{t}$$

The tangential damping force $\mathbf{F}_\mathrm{t,damp}$ is given by:

$$\mathbf{F}_\mathrm{t,damp} = -\eta_t \mathbf{v}_{t,rel}$$

The tangential damping prefactor $\eta_t$ is calculated by scaling the
normal damping $\eta_n$ (see above):

$$\eta_t = -x_{\gamma,t} \eta_n$$

The normal damping prefactor $\eta_n$ is determined by the choice of the
*damping* keyword, as discussed above. Thus, the *damping* keyword also
affects the tangential damping. The parameter $x_{\gamma,t}$ is a
scaling coefficient. Several works in the literature use
$x_{\gamma,t} = 1$ ([Marshall](Marshall2009), [Tsuji et al](Tsuji1992),
[Silbert et al](Silbert2001)). The relative tangential velocity at the
point of contact is given by
$\mathbf{v}_{t, rel} = \mathbf{v}_{t} - (R_i\mathbf{\Omega}_i + R_j\mathbf{\Omega}_j) \times \mathbf{n}$,
where
$\mathbf{v}_{t} = \mathbf{v}_r - \mathbf{v}_r\cdot\mathbf{n}\ \mathbf{n}$,
$\mathbf{v}_r = \mathbf{v}_j - \mathbf{v}_i$ . The direction of the
applied force is
$\mathbf{t} = \mathbf{v_{t,rel}}/\|\mathbf{v_{t,rel}}\|$ .

The normal force value $F_{n0}$ used to compute the critical force
depends on the form of the contact model. For non-cohesive models
(\*hertz\*, *hertz/material*, *hooke*), it is given by the magnitude of
the normal force:

$$F_{n0} = \|\mathbf{F}_n\|$$

For cohesive models such as *jkr* and *dmt*, the critical force is
adjusted so that the critical tangential force approaches
$\mu_t F_{pulloff}$, see [Marshall](Marshall2009), equation 43, and
[Thornton](Thornton1991). For both models, $F_{n0}$ takes the form:

$$F_{n0} = \|\mathbf{F}_{ne} + 2 F_{pulloff}\|$$

Where $F_{pulloff} = 3\pi \gamma R$ for *jkr*, and
$F_{pulloff} = 4\pi \gamma R$ for *dmt*.

The remaining tangential options all use accumulated tangential
displacement (i.e. contact history), except for the options
*mindlin/force* and *mindlin_rescale/force*, that use accumulated
tangential force instead, and are discussed further below. The
accumulated tangential displacement is discussed in details below in the
context of the *linear_history* option. The same treatment of the
accumulated displacement applies to the other options as well.

For *tangential linear_history*, the tangential force is given by:

$$\mathbf{F}_t =  -\min(\mu_t F_{n0}, \|-k_t\mathbf{\xi} + \mathbf{F}_\mathrm{t,damp}\|) \mathbf{t}$$

Here, $\mathbf{\xi}$ is the tangential displacement accumulated during
the entire duration of the contact:

$$\mathbf{\xi} = \int_{t0}^t \mathbf{v}_{t,rel}(\tau) \mathrm{d}\tau$$

This accumulated tangential displacement must be adjusted to account for
changes in the frame of reference of the contacting pair of particles
during contact. This occurs due to the overall motion of the contacting
particles in a rigid-body-like fashion during the duration of the
contact. There are two modes of motion that are relevant: the
\'tumbling\' rotation of the contacting pair, which changes the
orientation of the plane in which tangential displacement occurs; and
\'spinning\' rotation of the contacting pair about the vector connecting
their centers of mass ($\mathbf{n}$). Corrections due to the former mode
of motion are made by rotating the accumulated displacement into the
plane that is tangential to the contact vector at each step, or
equivalently removing any component of the tangential displacement that
lies along $\mathbf{n}$, and rescaling to preserve the magnitude. This
follows the discussion in [Luding](Luding2008), see equation 17 and
relevant discussion in that work:

$$\mathbf{\xi} = \left(\mathbf{\xi'} - (\mathbf{n} \cdot \mathbf{\xi'})\mathbf{n}\right) \frac{\|\mathbf{\xi'}\|}{\|\mathbf{\xi'} - (\mathbf{n}\cdot\mathbf{\xi'})\mathbf{n}\|}$$

Here, $\mathbf{\xi'}$ is the accumulated displacement prior to the
current time step and $\mathbf{\xi}$ is the corrected displacement.
Corrections to the displacement due to the second mode of motion
described above (rotations about $\mathbf{n}$) are not currently
implemented, but are expected to be minor for most simulations.

Furthermore, when the tangential force exceeds the critical force, the
tangential displacement is re-scaled to match the value for the critical
force (see [Luding](Luding2008), equation 20 and related discussion):

$$\mathbf{\xi} = -\frac{1}{k_t}\left(\mu_t F_{n0}\mathbf{t} - \mathbf{F}_{t,damp}\right)$$

The tangential force is added to the total normal force (elastic plus
damping) to produce the total force on the particle. The tangential
force also acts at the contact point (defined as the center of the
overlap region) to induce a torque on each particle according to:

$$\mathbf{\tau}_i = -(R_i - 0.5 \delta) \mathbf{n} \times \mathbf{F}_t$$

$$\mathbf{\tau}_j = -(R_j - 0.5 \delta) \mathbf{n} \times \mathbf{F}_t$$

For *tangential mindlin*, the [Mindlin](Mindlin1949) no-slip solution is
used which differs from the *linear_history* option by an additional
factor of $a$, the radius of the contact region. The tangential force is
given by:

$$\mathbf{F}_t =  -\min(\mu_t F_{n0}, \|-k_t a \mathbf{\xi} + \mathbf{F}_\mathrm{t,damp}\|) \mathbf{t}$$

Here, $a$ is the radius of the contact region, given by
$a =\sqrt{R\delta}$ for all normal contact models, except for *jkr*,
where it is given implicitly by
$\delta = a^2/R - 2\sqrt{\pi \gamma a/E}$, see discussion above. To
match the Mindlin solution, one should set $k_t = 8G_{eff}$, where
$G_{eff}$ is the effective shear modulus given by:

$$G_{eff} = \left(\frac{2-\nu_i}{G_i} + \frac{2-\nu_j}{G_j}\right)^{-1}$$

where $G$ is the shear modulus, related to Young\'s modulus $E$ and
Poisson\'s ratio $\nu$ by $G = E/(2(1+\nu))$. This can also be achieved
by specifying *NULL* for $k_t$, in which case a normal contact model
that specifies material parameters $E$ and $\nu$ is required (e.g.
*hertz/material*, *dmt* or *jkr*). In this case, mixing of the shear
modulus for different particle types *i* and *j* is done according to
the formula above.

:::: note
::: title
Note
:::

The radius of the contact region $a$ depends on the normal overlap. As a
result, the tangential force for *mindlin* can change due to a variation
in normal overlap, even with no change in tangential displacement.
::::

For *tangential mindlin/force*, the accumulated elastic tangential force
characterizes the contact history, instead of the accumulated tangential
displacement. This prevents the dependence of the tangential force on
the normal overlap as noted above. The tangential force is given by:

$$\mathbf{F}_t =  -\min(\mu_t F_{n0}, \|\mathbf{F}_{te} + \mathbf{F}_\mathrm{t,damp}\|) \mathbf{t}$$

The increment of the elastic component of the tangential force
$\mathbf{F}_{te}$ is given by:

$$\mathrm{d}\mathbf{F}_{te} = -k_t a \mathbf{v}_{t,rel} \mathrm{d}\tau$$

The changes in frame of reference of the contacting pair of particles
during contact are accounted for by the same formula as above, replacing
the accumulated tangential displacement $\xi$, by the accumulated
tangential elastic force $F_{te}$. When the tangential force exceeds the
critical force, the tangential force is directly re-scaled to match the
value for the critical force:

$$\mathbf{F}_{te} = - \mu_t F_{n0}\mathbf{t} + \mathbf{F}_{t,damp}$$

The same rules as those described for *mindlin* apply regarding the
tangential stiffness and mixing of the shear modulus for different
particle types.

The *mindlin_rescale* option uses the same form as *mindlin*, but the
magnitude of the tangential displacement is re-scaled as the contact
unloads, i.e. if $a < a_{t_{n-1}}$:

$$\mathbf{\xi} = \mathbf{\xi_{t_{n-1}}} \frac{a}{a_{t_{n-1}}}$$

Here, $t_{n-1}$ indicates the value at the previous time step. This
rescaling accounts for the fact that a decrease in the contact area upon
unloading leads to the contact being unable to support the previous
tangential loading, and spurious energy is created without the rescaling
above ([Walton](WaltonPC) ).

:::: note
::: title
Note
:::

For *mindlin*, a decrease in the tangential force already occurs as the
contact unloads, due to the dependence of the tangential force on the
normal force described above. By re-scaling $\xi$, *mindlin_rescale*
effectively re-scales the tangential force twice, i.e., proportionally
to $a^2$. This peculiar behavior results from use of the accumulated
tangential displacement to characterize the contact history. Although
*mindlin_rescale* remains available for historic reasons and backward
compatibility purposes, it should be avoided in favor of
*mindlin_rescale/force*.
::::

The *mindlin_rescale/force* option uses the same form as
*mindlin/force*, but the magnitude of the tangential elastic force is
re-scaled as the contact unloads, i.e. if $a < a_{t_{n-1}}$:

$$\mathbf{F}_{te} = \mathbf{F}_{te, t_{n-1}} \frac{a}{a_{t_{n-1}}}$$

This approach provides a better approximation of the
[Mindlin-Deresiewicz](Mindlin1953) laws and is more consistent than
*mindlin_rescale*. See discussions in [Thornton et al,
2013](Thornton2013), particularly equation 18(b) of that work and
associated discussion, and [Agnolin and Roux, 2007](AgnolinRoux2007),
particularly Appendix A.

------------------------------------------------------------------------

The optional *rolling* keyword enables rolling friction, which resists
pure rolling motion of particles. The options currently supported are:

1.  *none*
2.  *sds* : $k_{roll}$, $\gamma_{roll}$, $\mu_{roll}$

If the *rolling* keyword is not specified, the model defaults to *none*.

For *rolling sds*, rolling friction is computed via a
spring-dashpot-slider, using a \'pseudo-force\' formulation, as detailed
by [Luding](Luding2008). Unlike the formulation in
[Marshall](Marshall2009), this allows for the required adjustment of
rolling displacement due to changes in the frame of reference of the
contacting pair. The rolling pseudo-force is computed analogously to the
tangential force:

$$\mathbf{F}_{roll,0} =  k_{roll} \mathbf{\xi}_{roll}  - \gamma_{roll} \mathbf{v}_{roll}$$

Here,
$\mathbf{v}_{roll} = -R(\mathbf{\Omega}_i - \mathbf{\Omega}_j) \times \mathbf{n}$
is the relative rolling velocity, as given in [Wang et al](Wang2015) and
[Luding](Luding2008). This differs from the expressions given by [Kuhn
and Bagi](Kuhn2004) and used in [Marshall](Marshall2009); see [Wang et
al](Wang2015) for details. The rolling displacement is given by:

$$\mathbf{\xi}_{roll} = \int_{t_0}^t \mathbf{v}_{roll} (\tau) \mathrm{d} \tau$$

A Coulomb friction criterion truncates the rolling pseudo-force if it
exceeds a critical value:

$$\mathbf{F}_{roll} =  \min(\mu_{roll} F_{n,0}, \|\mathbf{F}_{roll,0}\|)\mathbf{k}$$

Here, $\mathbf{k} = \mathbf{v}_{roll}/\|\mathbf{v}_{roll}\|$ is the
direction of the pseudo-force. As with tangential displacement, the
rolling displacement is rescaled when the critical force is exceeded, so
that the spring length corresponds the critical force. Additionally, the
displacement is adjusted to account for rotations of the frame of
reference of the two contacting particles in a manner analogous to the
tangential displacement.

The rolling pseudo-force does not contribute to the total force on
either particle (hence \'pseudo\'), but acts only to induce an equal and
opposite torque on each particle, according to:

$$\tau_{roll,i} =  R_{eff} \mathbf{n} \times \mathbf{F}_{roll}$$

$$\tau_{roll,j} =  -\tau_{roll,i}$$

------------------------------------------------------------------------

The optional *twisting* keyword enables twisting friction, which resists
rotation of two contacting particles about the vector $\mathbf{n}$ that
connects their centers. The options currently supported are:

1.  *none*
2.  *sds* : $k_{twist}$, $\gamma_{twist}$, $\mu_{twist}$
3.  *marshall*

If the *twisting* keyword is not specified, the model defaults to
*none*.

For both *twisting sds* and *twisting marshall*, a history-dependent
spring-dashpot-slider is used to compute the twisting torque. Because
twisting displacement is a scalar, there is no need to adjust for
changes in the frame of reference due to rotations of the particle pair.
The formulation in [Marshall](Marshall2009) therefore provides the most
straightforward treatment:

$$\tau_{twist,0} = -k_{twist}\xi_{twist} - \gamma_{twist}\Omega_{twist}$$

Here $\xi_{twist} = \int_{t_0}^t \Omega_{twist} (\tau) \mathrm{d}\tau$
is the twisting angular displacement, and
$\Omega_{twist} = (\mathbf{\Omega}_i - \mathbf{\Omega}_j) \cdot \mathbf{n}$
is the relative twisting angular velocity. The torque is then truncated
according to:

$$\tau_{twist} = \min(\mu_{twist} F_{n,0}, \tau_{twist,0})$$

Similar to the sliding and rolling displacement, the angular
displacement is rescaled so that it corresponds to the critical value if
the twisting torque exceeds this critical value:

$$\xi_{twist} = \frac{1}{k_{twist}} (\mu_{twist} F_{n,0}sgn(\Omega_{twist}) - \gamma_{twist}\Omega_{twist})$$

For *twisting sds*, the coefficients $k_{twist}, \gamma_{twist}$ and
$\mu_{twist}$ are simply the user input parameters that follow the
*twisting sds* keywords in the *pair_coeff* command.

For *twisting_marshall*, the coefficients are expressed in terms of
sliding friction coefficients, as discussed in [Marshall](Marshall2009)
(see equations 32 and 33 of that work):

$$k_{twist} = 0.5k_ta^2$$

$$\eta_{twist} = 0.5\eta_ta^2$$

$$\mu_{twist} = \frac{2}{3}a\mu_t$$

Finally, the twisting torque on each particle is given by:

$$\mathbf{\tau}_{twist,i} = \tau_{twist}\mathbf{n}$$

$$\mathbf{\tau}_{twist,j} = -\mathbf{\tau}_{twist,i}$$

------------------------------------------------------------------------

If two particles are moving away from each other while in contact, there
is a possibility that the particles could experience an effective
attractive force due to damping. If the optional *limit_damping* keyword
is used, this option will zero out the normal component of the force if
there is an effective attractive force. This keyword cannot be used with
the JKR or DMT models.

------------------------------------------------------------------------

The optional *heat* keyword enables heat conduction. The options
currently supported are:

1.  *none*
2.  *area* : $k_{s}$

If the *heat* keyword is not specified, the model defaults to *none*.

For *heat* *area*, the heat $Q$ conducted between two particles is given
by

$$Q = k_{s} A \Delta T$$

where $\Delta T$ is the difference in the two particles\' temperature,
$k_{s}$ is a non-negative numeric value for the conductivity, and $A$ is
the area of the contact and depends on the normal force model.

Note that the option *none* must either be used in all or none of of the
*pair_coeff* calls. See [fix heat/flow](fix_heat_flow) and [fix
property/atom](fix_property_atom) for more information on this option.

------------------------------------------------------------------------

The *granular* pair style can reproduce the behavior of the *pair
gran/\** styles with the appropriate settings (some very minor
differences can be expected due to corrections in displacement history
frame-of-reference, and the application of the torque at the center of
the contact rather than at each particle). The first example above is
equivalent to *pair gran/hooke 1000.0 NULL 50.0 50.0 0.4 1*. The second
example is equivalent to *pair gran/hooke/history 1000.0 500.0 50.0 50.0
0.4 1*. The third example is equivalent to *pair gran/hertz/history
1000.0 500.0 50.0 50.0 0.4 1 limit_damping*.

------------------------------------------------------------------------

LAMMPS automatically sets pairwise cutoff values for *pair_style
granular* based on particle radii (and in the case of *jkr* pull-off
distances). In the vast majority of situations, this is adequate.
However, a cutoff value can optionally be appended to the *pair_style
granular* command to specify a global cutoff (i.e. a cutoff for all atom
types). Additionally, the optional *cutoff* keyword can be passed to the
*pair_coeff* command, followed by a cutoff value. This will set a
pairwise cutoff for the atom types in the *pair_coeff* command. These
options may be useful in some rare cases where the automatic cutoff
determination is not sufficient, e.g. if particle diameters are being
modified via the *fix adapt* command. In that case, the global cutoff
specified as part of the *pair_style granular* command is applied to all
atom types, unless it is overridden for a given atom type combination by
the *cutoff* value specified in the *pair coeff* command. If *cutoff* is
only specified in the *pair coeff* command and no global cutoff is
appended to the *pair_style granular* command, then LAMMPS will use that
cutoff for the specified atom type combination, and automatically set
pairwise cutoffs for the remaining atom types.

------------------------------------------------------------------------

Styles with a *gpu*, *intel*, *kk*, *omp*, or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the [Accelerator packages](Speed_packages)
page. The accelerated styles take the same arguments and should produce
the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, INTEL, KOKKOS, OPENMP, and
OPT packages, respectively. They are only enabled if LAMMPS was built
with those packages. See the [Build package](Build_package) page for
more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the [-suffix command-line
switch](Run_options) when you invoke LAMMPS, or you can use the
[suffix](suffix) command in your input script.

See the [Accelerator packages](Speed_packages) page for more
instructions on how to use the accelerated styles effectively.

------------------------------------------------------------------------

## Mixing, shift, table, tail correction, restart, rRESPA info

The [pair_modify](pair_modify) mix, shift, table, and tail options are
not relevant for granular pair styles.

Mixing of coefficients is carried out using geometric averaging for most
quantities, e.g. if friction coefficient for type 1-type 1 interactions
is set to $\mu_1$, and friction coefficient for type 2-type 2
interactions is set to $\mu_2$, the friction coefficient for type1-type2
interactions is computed as $\sqrt{\mu_1\mu_2}$ (unless explicitly
specified to a different value by a *pair_coeff 1 2 \...* command). The
exception to this is elastic modulus, only applicable to
*hertz/material*, *dmt* and *jkr* normal contact models. In that case,
the effective elastic modulus is computed as:

$$E_{eff,ij} = \left(\frac{1-\nu_i^2}{E_i} + \frac{1-\nu_j^2}{E_j}\right)^{-1}$$

If the *i-j* coefficients $E_{ij}$ and $\nu_{ij}$ are explicitly
specified, the effective modulus is computed as:

$$E_{eff,ij} = \left(\frac{1-\nu_{ij}^2}{E_{ij}} + \frac{1-\nu_{ij}^2}{E_{ij}}\right)^{-1}$$

or

$$E_{eff,ij} = \frac{E_{ij}}{2(1-\nu_{ij}^2)}$$

These pair styles write their information to [binary restart
files](restart), so a pair_style command does not need to be specified
in an input script that reads a restart file.

These pair styles can only be used via the *pair* keyword of the
[run_style respa](run_style) command. They do not support the *inner*,
*middle*, *outer* keywords.

The single() function of these pair styles returns 0.0 for the energy of
a pairwise interaction, since energy is not conserved in these
dissipative potentials. It also returns only the normal component of the
pairwise interaction force. However, the single() function also
calculates 12 extra pairwise quantities. The first 3 are the components
of the tangential force between particles I and J, acting on particle I.
The fourth is the magnitude of this tangential force. The next 3 (5-7)
are the components of the rolling torque acting on particle I. The next
entry (8) is the magnitude of the rolling torque. The next entry (9) is
the magnitude of the twisting torque acting about the vector connecting
the two particle centers. The last 3 (10-12) are the components of the
vector connecting the centers of the two particles (x_I - x_J).

These extra quantities can be accessed by the [compute
pair/local](compute_pair_local) command, as *p1*, *p2*, \..., *p12*.

------------------------------------------------------------------------

## Restrictions

This pair style is part of the GRANULAR package. It is only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

This pair style requires that atoms store per-particle radius, torque,
and angular velocity (omega) as defined by the [atom_style
sphere](atom_style).

This pair style requires you to use the [comm_modify vel
yes](comm_modify) command so that velocities are stored by ghost atoms.

This pair style will not restart exactly when using the
[read_restart](read_restart) command, though it should provide
statistically similar results. This is because the forces it computes
depend on atom velocities and the atom velocities have been propagated
half a timestep between the force computation and when the restart is
written, due to using Velocity Verlet time integration. See the
[read_restart](read_restart) command for more details.

Accumulated values for individual contacts are saved to restart files
but are not saved to data files. Therefore, forces may differ
significantly when a system is reloaded using the [read_data](read_data)
command.

## Related commands

[pair_coeff](pair_coeff) [pair gran/\*](pair_gran)

## Default

For the *pair_coeff* settings: *damping viscoelastic*, *rolling none*,
*twisting none*.

## References

::: {#Brill1996}
**(Brilliantov et al, 1996)** Brilliantov, N. V., Spahn, F., Hertzsch,
J. M., & Poschel, T. (1996). Model for collisions in granular gases.
Physical review E, 53(5), 5382.
:::

::: {#Tsuji1992}
**(Tsuji et al, 1992)** Tsuji, Y., Tanaka, T., & Ishida, T. (1992).
Lagrangian numerical simulation of plug flow of cohesionless particles
in a horizontal pipe. Powder technology, 71(3), 239-250.
:::

::: {#JKR1971}
**(Johnson et al, 1971)** Johnson, K. L., Kendall, K., & Roberts, A. D.
(1971). Surface energy and the contact of elastic solids. Proc. R. Soc.
Lond. A, 324(1558), 301-313.
:::

::: {#DMT1975}
**Derjaguin et al, 1975)** Derjaguin, B. V., Muller, V. M., & Toporov,
Y. P. (1975). Effect of contact deformations on the adhesion of
particles. Journal of Colloid and interface science, 53(2), 314-326.
:::

::: {#Luding2008}
**(Luding, 2008)** Luding, S. (2008). Cohesive, frictional powders:
contact models for tension. Granular matter, 10(4), 235.
:::

::: {#Marshall2009}
**(Marshall, 2009)** Marshall, J. S. (2009). Discrete-element modeling
of particulate aerosol flows. Journal of Computational Physics, 228(5),
1541-1561.
:::

::: {#Silbert2001}
**(Silbert, 2001)** Silbert, L. E., Ertas, D., Grest, G. S., Halsey, T.
C., Levine, D., & Plimpton, S. J. (2001). Granular flow down an inclined
plane: Bagnold scaling and rheology. Physical Review E, 64(5), 051302.
:::

::: {#Kuhn2004}
**(Kuhn and Bagi, 2005)** Kuhn, M. R., & Bagi, K. (2004). Contact
rolling and deformation in granular media. International journal of
solids and structures, 41(21), 5793-5820.
:::

::: {#Wang2015}
**(Wang et al, 2015)** Wang, Y., Alonso-Marroquin, F., & Guo, W. W.
(2015). Rolling and sliding in 3-D discrete element models.
Particuology, 23, 49-55.
:::

::: {#Thornton1991}
**(Thornton, 1991)** Thornton, C. (1991). Interparticle sliding in the
presence of adhesion. J. Phys. D: Appl. Phys. 24 1942
:::

::: {#Mindlin1949}
**(Mindlin, 1949)** Mindlin, R. D. (1949). Compliance of elastic bodies
in contact. J. Appl. Mech., ASME 16, 259-268.
:::

::: {#Thornton2013}
**(Thornton et al, 2013)** Thornton, C., Cummins, S. J., & Cleary, P. W.
(2013). An investigation of the comparative behavior of alternative
contact force models during inelastic collisions. Powder Technology,
233, 30-46.
:::

::: {#WaltonPC}
**(Otis R. Walton)** Walton, O.R., Personal Communication
:::

::: {#Mindlin1953}
**(Mindlin and Deresiewicz, 1953)** Mindlin, R.D., & Deresiewicz, H
(1953). Elastic Spheres in Contact under Varying Oblique Force. J. Appl.
Mech., ASME 20, 327-344.
:::

::: {#AgnolinRoux2007}
**(Agnolin and Roux 2007)** Agnolin, I. & Roux, J-N. (2007). Internal
states of model isotropic granular packings. I. Assembling process,
geometry, and contact networks. Phys. Rev. E, 76, 061302.
:::
