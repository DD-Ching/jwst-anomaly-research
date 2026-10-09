# R2-A: Co-Witness Consistency (CWC)

*System A, round 2, group R2-A. Status: **speculative hypothesis**. Nothing here is an observation.
I label each claim as follows: [AX] axiom, [DER] derived from the axioms (the argument is given),
[EST] numerical estimate using standard GR/ΛCDM inputs, [HYP] an extra assumption that the axioms
do not force.*

## 0. Idea in one paragraph

The primitives are records, plus one requirement: records that can ever be compared must agree.
Spacetime is whatever structure makes that requirement satisfiable. The axioms give the theory
**no kinematic scale and no notion of nearest neighbour**. The first fact makes the light cone
exact, with no energy-dependent deviation. The second avoids the finite-valency obstruction. The
one place where the theory departs from GR is a single strong form of the consistency requirement:
**any two records must have a common future witness**, meaning some later record that both of them
can reach. GR allows event horizons, which are pairs of events with no common future. CWC forbids
them. The "connection between distant regions" is therefore a constraint that is enforced where
the regions' futures meet. It is not a path between them. The derived consequence that can be
tested is a closure rule on energy at black-hole remnants, and its size is fixed by GR's own
horizon-absorbed flux.

## 1. Axioms

**A0 (records, witnessing) [AX].** There is a countable set $R$ of records and a strict partial
order $\prec$ on it, read "$r \prec s$: $s$ witnesses $r$", meaning $s$'s content depends on $r$'s.
The order is irreflexive and transitive, so it is acyclic: no record witnesses itself.

**A1 (density, no valency) [AX].** For every $r \prec s$, the interval $I(r,s)=\{x: r\prec x\prec s\}$
is infinite, and no record has an immediate successor. Records are dense in the order, the way
$\mathbb{Q}$ is dense in $\mathbb{R}$. They are not locally finite. There is **no "neighbour"
relation and no valency**.

**A2 (consistency functor) [AX].** Each record has *content* in a fixed space $V$: frequencies,
directions and counts. Each witness relation $r\prec s$ carries a transfer
$T_{sr}\in G=\mathrm{Aut}(V)$ that converts $r$'s content into $s$'s terms. For a chain
$r\prec x\prec s$, the rule is $T_{sr}=T_{sx}T_{xr}$. Two different chains from $r$ to $s$ may
differ only by a *holonomy* $H$. $H$ must itself be recorded by the records enclosed between the
two chains (see A5): no unrecorded inconsistency is allowed.

**A3 (relativity) [AX].** An *observer* is a maximal chain. The order-automorphisms of $(R,\prec)$
act transitively on locally unaccelerated observers, so no observer is preferred. Reciprocity
follows: the frequency ratio $k_{AB}$ that $B$ records for $A$'s ticks equals $k_{BA}$.

**A4 (co-witnessability; this is the consistency requirement) [AX].** "All observers' records are
mutually consistent" can only be checked where both records are present. So it is formalised as:
$(R,\prec)$ is **upward directed**. For all $r,s$ there is a record $w$ with $r\preceq w$ and
$s\preceq w$.

**A5 (carrier conservation) [AX].** A witness relation is realised by *carriers*, and each carrier
is counted once at emission and once at absorption. Carriers are not created or destroyed between
records. Every record that can be counted is the end of a carrier.

**A6 (no kinematic scale) [AX].** The axioms that define $\prec$ and $T$ contain no dimensionful
constant. The only dimensionful number in the theory is one *coupling* $a_*$, an area that converts
carrier-count density into holonomy (§3.4). It appears only on the source side.

## 2. Light cone and Lorentz invariance

**2.1 Cone = boundary of the order [DER].** For an observer $O$ and a record $r$, define $e(r)$ as
the earliest record on $O$ that witnesses $r$. A1 makes this an infimum, not a "next" record. The
*light cone* of $r$ is the order-boundary $\partial J^+(r)$: the records that witness $r$ and are
not interior to $J^+(r)$. This definition uses $\prec$ alone, so every observer's light cone is the
same set. The cone is invariant by construction, not to some approximation.

**2.2 The group is Lorentz [DER].** Radar records give each observer coordinates (Bondi-style).
$O$ sends at its own tick $t_1$, the cone-boundary carrier reaches $r$, and the reply returns at
$t_2$. Then $t=(t_1+t_2)/2$ and $x=(t_2-t_1)/2$, where the factor that converts tick counts to
length is a unit choice, not a constant (A6). For collinear observers A2 forces
$k_{AC}=k_{AB}k_{BC}$, and A3 forces $k_{AB}=k_{BA}$. Together these give the composition law for
velocities, $v=(k^2-1)/(k^2+1)$, which is relativistic. More generally, homogeneity plus isotropy
plus the group property leave a one-parameter family of kinematic groups with parameter $\kappa$:

* $\kappa>0$ is rotation-like and mixes past and future. It does not preserve $\prec$, so it is
  excluded by A3.
* $\kappa=0$ is Galilean. Its cone is a hyperplane: two simultaneous records would witness each
  other, which gives $r\prec s\prec r$. This is excluded by A0 (acyclicity).
* Only $\kappa<0$ survives, which is the Lorentz group.

The invariant speed has no separate status. It is the slope of $\partial J^+$ in radar units.

**2.3 Escaping the finite-valency obstruction [DER].** The no-go theorem applies to locally finite
graphs whose vertices have finitely many neighbours, and boosts must then stretch the edge set. A1
removes the premise: there are no neighbours to stretch. One concrete model is the set of rational
points of Minkowski space with $\prec$ the causal order. Its automorphism group contains
$SO^+(3,1;\mathbb{Q})$, which is **dense** in $SO^+(3,1)$, and it maps the null cone
$t^2=|\mathbf{x}|^2$ exactly onto itself. Any finite-precision measurement cannot tell a dense
subgroup from the full group, so no preferred frame can be detected.

**2.4 Why the bounds hold, with predicted violation exactly zero [DER].** A Lorentz-violating
dispersion relation $\omega^2=k^2(1+\xi(k/k_*)^n)$ needs a scale $k_*$. A6 provides none. The
coupling $a_*$ enters only sources (§3.4) and never $\prec$, so it cannot reach the dispersion of
free carriers. Every carrier that realises boundary links therefore has $\omega=|k|$ in radar units.

* **Gravity.** In CWC, gravity is the holonomy of $T$ (§3.4), and holonomy changes propagate
  along the same $\partial J^+$. So $c_{\rm g}=c_\gamma$ **identically**. GW170817 bounds
  $|c_g-c|/c\lesssim10^{-15}$; CWC predicts zero.
* **Laboratory.** Anisotropy or boost-dependence at $10^{-18}$ would need a preferred frame. A3
  forbids one, and the dense automorphism group leaves no residual direction.

The bounds are satisfied because the predicted violation is zero. It is not small and tuned. The
attack surface is A6: if anything forces $a_*$ into the kinematics, this argument fails (§6).

## 3. Flux, distances, conservation

**3.1 Inverse square [DER].** Take an isotropic emitter in a holonomy-free region. The radar-sphere
of radius $r$ has area $4\pi r^2$, because §2 gives a Minkowski metric locally. A5 says the number
of carriers crossing every sphere equals the number emitted. So the count per unit area falls as
$1/(4\pi r^2)$, and each count arrives with frequency ratio $k$. Flux is
$F=L/(4\pi r^2)\cdot k^{-2}$, with one power of $k$ from energy per carrier and one from arrival rate.

**3.2 The four distances [DER].** Each comes from a different kind of record:

* **Radar (proper) distance**, from round-trip counts. Over cosmological baselines, relays of
  comoving observers compose $k$ multiplicatively (A2), which gives the FLRW form $1+z=a_0/a$ once
  homogeneity and isotropy are imposed.
* **Angular-diameter distance $d_A$**, from direction records.
* **Luminosity distance $d_L$**, from counts and A5.
* **Parallax / proper-motion distance**, from direction records made by two observers.

Reciprocity has a direct reason here. A carrier is a *two-ended record*. The solid angle the
emitter assigns to the receiver and the one the receiver assigns to the emitter describe the same
bundle of witness relations, and A2 converts between them with $k$. Hence
$d_L=(1+z)^2d_A$ **exactly**. CWC predicts no cosmic-opacity signal: lesson 3 is respected, not
exploited.

**3.3 Conservation [DER].** Carrier-count conservation is primitive (A5). Energy-momentum
conservation is Noether's theorem for the automorphisms of $(R,\prec)$: the translations are among
them. Angular momentum follows from rotations.

**3.4 Gravity [DER, with HYP].**

1. *Curvature is holonomy.* In A2, holonomy is the failure of two chains to agree. A2 also requires
   that this failure be recorded. For infinitesimal loops this gives a curvature 2-form.
2. *Divergence-free source.* "Consistency of consistency" means the holonomies of the faces of any
   closed 3-cell compose to the identity. That is the Bianchi identity, so the source tensor must
   be divergence-free.
3. *Einstein equations.* Set the source equal to the carrier energy-momentum, using $a_*$
   [HYP: linear, minimal coupling]. Lovelock-type uniqueness in four dimensions then gives the
   Einstein equations, with a cosmological term entering **only as an integration constant**.
   A6 forbids a fundamental $\Lambda$.
4. *The coupling.* $a_*$ is identified with $8\pi G\hbar/c^3$, which is $8\pi\ell_P^2$.

## 4. Where CWC departs from GR + ΛCDM

**4.1 The structural deviation [DER from A4].** GR allows event horizons: pairs $(r,s)$ with
$J^+(r)\cap J^+(s)=\varnothing$. A4 forbids **all** of them. This is the "non-shortcut connection":
two regions that never exchange a carrier are still bound by the condition that their futures
intersect, and their records must agree there. No path is added and no distance shrinks.

**4.2 Cosmology [DER; no present-day amplitude].** Using ΛCDM with $\Omega_m=0.31$ and
$H_0=67.7$ [EST]:

* The comoving event horizon today is $\chi_{eh}\approx5.08$ Gpc.
* The particle horizon is $\chi_p\approx14.2$ Gpc.
* "Now" events on our worldline and on a comoving worldline at $\chi$ share a future only if
  $\chi\le2\chi_{eh}\approx10.2$ Gpc. That corresponds to $z\gtrsim13$ for sources we already see.

So ΛCDM contains pairs of record-holders, both inside our past light cone, that have no common
witness. CWC instead requires the remaining conformal time to diverge,
$\int_{t_0}^{\infty}c\,dt/a=\infty$. Dark energy therefore cannot be an eternal constant. A6 forbids
it as a fundamental constant anyway. The present phase of acceleration must end, with
$w_{\rm eff}\to\ge-1/3$ asymptotically. **Honest status:** this is an inequality about the future.
It sets no amplitude today, so by lesson 4 it is **not** offered as a prediction.

**4.3 Compact objects: the derived amplitude [DER].**

* *No permanent trapping.* By A4, carriers that cross the would-be event horizon of a merger
  remnant must eventually become co-witnessable with the outside. By A5, all of them must come back
  out: none are lost.
* *Where they come back.* A6 allows only one length to place the turning surface: $\ell_P$.
  So the surface sits at proper distance $\zeta\ell_P$ outside the trapping surface, with
  $\zeta=O(1)$ undetermined.
* *Delay.* The round trip from the light-ring barrier to that surface takes
  $$\Delta t\simeq\frac{4GM}{c^3}\Big(1+\frac{1}{\sqrt{1-\chi^2}}\Big)\ln\frac{GM}{c^2\ell_P}\;(1\pm0.02).$$
  The unknown $\zeta$ appears only inside the logarithm. Changing $\zeta$ by a factor of $e$ shifts
  $\Delta t$ by about 1 %, so the delay is fixed by the remnant's mass $M$ and spin $\chi$ alone.
  Values [EST]:

  | Remnant (detector frame) | $\Delta t$ |
  |---|---|
  | GW150914-like, $M\approx68M_\odot$, $\chi\approx0.67$ | 0.29 s |
  | $20M_\odot$ | 0.085 s |
  | $150M_\odot$ | 0.66 s |

* *Energy closure (the amplitude).* By A5, the total energy re-emitted to infinity, summed over all
  later leakage through the barrier, must equal the energy GR sends **into** the horizon:
  $$\sum_n E_n^{\rm after}=E_{\rm hor}^{\rm GR}(m_1,m_2,\vec\chi_1,\vec\chi_2).$$
  $E_{\rm hor}^{\rm GR}$ is a deterministic, GR-computable number. It can be read from
  numerical-relativity horizon-flux or apparent-horizon mass histories, or from the ingoing part of
  the ringdown in perturbation theory. It is **not** a free parameter and cannot be tuned to zero.
  CWC adds nothing to it. GR predicts that this energy never returns; CWC predicts all of it
  returns.
* *Phase.* A4 and A5 fix the number and energy of returning carriers but not their phase. Infinite
  redshift at the trapping surface makes the transfer $T$ there ill-defined. The robust prediction
  is therefore **energy closure within a delay window, with the spectrum filtered by the barrier**:
  it is concentrated near the remnant's fundamental quasinormal-mode frequency.
* *Coherent case.* Comb-like echoes are the special case where phase survives. CWC does not require
  it.

## 5. Predictions and kill criteria

**P1: the energy-closure afterglow** (meets requirements a, b, c).

* *What is predicted.* After every binary-black-hole ringdown, excess gravitational-wave power
  appears at times $t_{\rm rd}+[\Delta t, \text{a few}\times\Delta t]$. It is band-limited around
  $f_{220}$, and its integrated energy equals $E^{\rm GR}_{\rm hor}$.
* *(a) Not mimicked by ordinary astrophysics.* Environmental matter produces no delay that scales
  as $M\ln(M/\ell_P)$ and no band-limit at $f_{220}$. Instrumental glitches are not coherent across
  a detector network with the astrophysical time offsets and antenna patterns.
* *(b) Data.* Public LIGO/Virgo/KAGRA strain from the GWOSC open-data releases (O1 to O4
  catalogue events), together with public posterior samples for $M_f$ and $\chi_f$. Method: a
  coherent excess-power stack across events, in event-specific windows. **No template parameters
  are fitted.**
* *(c) Kill criterion.* First compute $E^{\rm GR}_{\rm hor}$ per event from numerical-relativity
  surrogates. Then form the predicted stacked SNR. If that SNR is ≥ 5 and the observed stack is
  < 1/3 of the predicted energy at 95 % confidence, CWC is **dead**. If the predicted SNR is < 5,
  the result must be reported as **untested, not survived**.
* *Pre-registered worry.* $E^{\rm GR}_{\rm hor}$ may be only percent-level of the radiated energy.
  If so, only the loudest events (SNR ≳ 50) contribute meaningfully.

**P2: zero for all Lorentz and speed-of-gravity tests.** CWC predicts:

* no energy-dependent photon arrival times in GRB, blazar or FRB data (Fermi, H.E.S.S./MAGIC
  public light curves);
* $c_g=c$ exactly.

Any significant, reproducible detection kills CWC through A6.

**P3: Etherington ratio $\eta\equiv d_L/[(1+z)^2d_A]=1$ exactly** (supernova plus BAO/cluster
data). Any significant deviation kills A5.

**Quantities current searches never measure:**

1. **The energy-closure ratio $\mathcal{R}=\sum E^{\rm after}/E^{\rm GR}_{\rm hor}$.** Published
   echo searches fit free reflectivity and damping to coherent templates. None compares
   integrated late-time power with GR's own horizon-absorbed energy. CWC predicts $\mathcal{R}=1$;
   GR predicts $\mathcal{R}=0$.
2. **The remaining conformal reach $\Delta\eta_\infty=\int_{t_0}^\infty c\,dt/a$.** It can be
   reported from any dark-energy posterior, $w(z)$ or $(w_0,w_a)$ extrapolated. Equivalently: the
   fraction of currently observed sources that are co-witnessable with us. ΛCDM gives
   $\Delta\eta_\infty\approx5.1$ Gpc comoving; CWC requires $\infty$. This is extrapolation, not a
   test, but nobody reports it, and it would make the A4 tension explicit.

## 6. The three weakest points

1. **The ergoregion and spinning remnants.** A spinning compact object that returns everything
   that falls in is known to be prone to superradiant ergoregion instability. If phase-scrambled
   re-emission does not quench it, CWC predicts that rapidly spinning remnants spin down quickly.
   That conflicts with observed high black-hole spins in X-ray binaries and with GW remnant spins
   near 0.7. This could kill P1's premise **before** any data analysis. I have not shown that
   incoherent return prevents the growth.
2. **A6 is asserted, not earned.** The coupling $a_*\sim\ell_P^2$ is a scale. I claimed it enters
   only on the source side. In any quantum version, loops of carriers feed sources back into
   propagation, which is how $a_*$ could reach the dispersion relation and reopen Lorentz
   violation. The exactness argument in §2.4 holds only at the classical level. Also, A1 (dense
   records) does not sit easily with counting: carriers must be discrete while records are dense,
   and I have only stipulated this.
3. **Thin and partly borrowed derivation of gravity and cosmology.** Going from "holonomy must be
   recorded" to the Einstein equations uses a minimal-coupling assumption [HYP] plus a uniqueness
   theorem. Nothing in CWC derives $G$, the value of $\Lambda$ today, or $n_s$. The cosmological
   consequence (§4.2) is only an inequality about the future. The surviving fixed-amplitude
   prediction, P1, sits in strong gravity. In the cosmological sector the theory has no present-day
   amplitude at all.

## Resemblances I noticed

* **Bondi k-calculus and Ignatowski-type derivations of Lorentz kinematics:** used in §2.2.
* **Causal-order-determines-conformal-geometry theorems** (Malament; Hawking–King–McCarthy):
  behind §2.1.
* **Causal set theory:** shares order primitives. A1 deliberately rejects local finiteness, which
  is the opposite choice.
* **Unimodular gravity:** $\Lambda$ as an integration constant (§3.4).
* **Weyl tensor / 2+1-gravity holonomy:** intuitions behind "curvature is holonomy".
* **Black-hole echoes and exotic-compact-object models** (Cardoso, Pani, Abedi et al.): the delay
  formula in §4.3 is theirs. My 0.29 s for GW150914 matches their published estimate. What is new
  here is only the energy-closure amplitude derived from A5, and the claim that phase is not fixed.
* **The firewall/information-loss debates and the "no event horizon" proposals** (e.g. apparent-
  horizon-only black holes): A4 resembles these.
* **Vacuum-energy sequestering and "no eternal de Sitter" conjectures:** both resemble the §4.2
  consequence.
* **Relational and observer-consistency views of quantum mechanics:** resemble the seed.
* **Distance-duality tests:** P3 is the standard test.
