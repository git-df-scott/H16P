# LIT_B — Global classifications of quadratic systems, and which graphics (separatrix cycles) can surround a focus

Slice B of the H(2) literature sweep. Compiled 2026-09-06.
Everything marked **VERBATIM** is copied from the PDF text layer (ligature loss from `pypdf`
extraction has been repaired silently; `fi`/`fl` were dropped by the extractor in the Cambridge
scan of Llibre–Schlomiuk 2004). Everything marked *(paraphrase / secondary)* is NOT verbatim
and must be re-checked against the source before being used in a written argument.

PDF root: `/Users/scottg/Claude_all/papers/`

> **Housekeeping note for the other agents:** a file that was already in `papers/` under the name
> `cambridge_quadratic_weak_focus.pdf` was identified as Zhang Pingguang & Cai Suilin, *Quadratic
> systems with a weak focus*, Bull. Austral. Math. Soc. 44 (1991) 511–526, and **renamed** to
> `zhang_pingguang_cai_suilin_1991_quadratic_systems_weak_focus_BAMS44.pdf`.

---

## 0. Executive summary (the part that matters for the >=5 limit cycle hunt)

Five theorems in this slice, taken together, close off almost every "stack a weak focus on top of
a surrounding graphic" recipe:

1. **Ye Yanqian.** A limit cycle of a quadratic system surrounds a *unique* singular point, and
   that point is a focus. (So all counting is per-focus, and a quadratic system has at most two
   "limit-cycle-carrying" points, since it has at most two anti-saddles that are foci.)
2. **Li Chengzhi (1986).** There are **no** limit cycles surrounding a weak focus of order three.
3. **Zhang Pingguang (2001/2002).** If a quadratic system has limit cycles around **two** foci,
   then around one of them there is **at most one**. (=> distribution is (1,n); (2,2) and worse
   are dead.)
4. **Zhang Pingguang–Cai Suilin (1991).** A quadratic system with a weak focus and a strong focus
   cannot have (2,2) distribution; and with 4 finite singular points (strong focus + weak focus +
   2 saddles) the strong focus carries at most one limit cycle.
5. **Kooij–Zegeling (1998).** With exactly two finite singularities (a weak focus + a strong
   focus) and >=2 singularities at infinity: if the weak focus is of order **2 or 3**, the strong
   focus is surrounded by **at most one** limit cycle, and it is hyperbolic.

Consequence for the classification papers below: in `QW3` the unique graphic that occurs
(portraits W13, W15, W18) surrounds the **strong** focus, never the weak focus of order three,
and its cyclicity is **1**. So the QW3 recipe caps at 3 (Bautin) + 1 (graphic) = **4**.
This is exactly Llibre–Schlomiuk's Theorem 18 / Corollary 19.

The live gaps are (i) the **weak focus of order 1** families (QW1, QW1SN, QW1SN11, QW1SN02,
QW^I), where Coppel's Conjecture W is still open and Artés-school bifurcation diagrams show parts
"with at least two limit cycles"; and (ii) the ~30 DRR graphics whose finite cyclicity is still
unproved — but note that finite cyclicity being *unproved* is not the same as cyclicity being
*large*, and no open DRR graphic has a proven cyclicity >= 2 in the quadratic family.

---

## 1. Llibre & Schlomiuk 2004 — "The geometry of quadratic differential systems with a weak focus of third order"

**Citation.** J. Llibre and D. Schlomiuk, *The geometry of quadratic differential systems with a
weak focus of third order*, Canad. J. Math. **56** (2004), no. 2, 310–343. DOI 10.4153/CJM-2004-015-2.

**PDF: OBTAINED (open "bronze" copy on Cambridge Core).**
`/Users/scottg/Claude_all/papers/llibre_schlomiuk_2004_CJM56_weak_focus_third_order.pdf` (34 pp.)
Direct URL that worked (found via the Semantic Scholar `openAccessPdf` field, not via the article
landing page, which is paywalled):
`https://www.cambridge.org/core/services/aop-cambridge-core/content/view/FCA46BEF322C4B8A02C2C47270177DF5/S0008414X00031758a.pdf/div-class-title-the-geometry-of-quadratic-differential-systems-with-a-weak-focus-of-third-order-div.pdf`

### 1.1 Abstract (VERBATIM)

> "In this article we determine the global geometry of the planar quadratic differential systems
> with a weak focus of third order. This class plays a significant role in the context of Hilbert's
> 16-th problem. Indeed, all examples of quadratic differential systems with at least four limit
> cycles, were obtained by perturbing a system in this family. We use the algebro-geometric
> concepts of divisor and zero-cycle to encode global properties of the systems and to give
> structure to this class. We give a theorem of topological classification of such systems in terms
> of integer-valued affine invariants. According to the possible values taken by them in this
> family we obtain a total of 18 topologically distinct phase portraits. We show that inside the
> class of all quadratic systems with the topology of the coefficients, there exists a neighborhood
> of the family of quadratic systems with a weak focus of third order and which may have graphics
> but no polycycle in the sense of [15] and no limit cycle, such that any quadratic system in this
> neighborhood has at most four limit cycles."

### 1.2 THEOREM 16 (VERBATIM) — the requested list of portraits with graphics

> "**Theorem 16** The class QW3 is partitioned in the following three subclasses:
>
> (I) Systems without a limit cycle and without graphics. We have a total of twelve such phase
> portraits: Wi with i = 1, ..., 12. The systems with Wi, i = 1, ..., 10 are classified by J(S)
> (see Theorem 13). J(S) = (7, 2, 2, 3) for W11 and for W12 and these phase portraits are
> distinguished by m1. m1 = 1 for W11 and m1 = 2 for W12.
>
> (II) Systems with a limit cycle. These have no graphic and the limit cycle is unique. These
> yield three phase portraits which are topologically classified by J(S). More precisely we have:
> (II.1) J(S) = (7, 2, 2, 1) with phase portrait W14; (II.2) J(S) = (7, 2, 2, 3) with phase
> portrait W17; (II.3) J(S) = (6, 2, 2, 2) with phase portrait W16 occurring as a bifurcation from
> (II.2) to (II.1), when two of the three points at infinity collide. (In Figure 3 the region where
> we have limit cycles is delimited by the curves G and a = 0: Σ ∈ R9 ∪ Σ3 ∪ R8^dd).
>
> (III) **Systems with a graphic. These have no limit cycle and the graphic is unique, surrounding
> a strong focus.** We have exactly three phase portraits in this class. These are: W18 with
> J(S) = (7, 2, 2, 3) and W13, W15, both with J(S) = (6, 2, 2, 2). W13 is distinguished from W15
> by m1. For W13, m1 = 3; for W15, m1 = 2.
>
> We note that if J(S) = (7, 2, 2, 3), S could be of any of the types (I), (II), (III)."

**(a) Graphics surrounding a weak focus of order 3: NONE.** Theorem 16(III) is explicit — the
graphic "surround[s] a strong focus". Cross-check, §7.1(viii) VERBATIM:

> "(viii) There are no limit cycles of systems (6) surrounding the weak focus of third order,
> see [24]."
> ([24] = Chengzhi Li, *Non-existence of limit cycles around a weak focus of order three for any
> quadratic system*, Chinese Ann. Math. Ser. B **7** (1986), 174–190.)

**(b) The graphics all go through infinity (VERBATIM, §1):**

> "It is in B that a connected bifurcation curve G, part analytic, part algebraic, occurs on the
> diagram of Figure 3. On G the systems do not have limit cycles. **G is the only subset in the
> parameter space where graphics occur; more precisely, for all points in G the systems have a
> unique graphic (which in two of the three possible cases is a polycycle (cf. [15])) with two
> singularities, both at infinity and with two path curves, one of them part of the line at
> infinity.**"

Structure of the three graphics (VERBATIM, proof of Prop. 17):
- **W13** (parameters on the algebraic curve G1): "elementary (in the sense of [15]), having two
  singular points, a saddle and a saddle-node".
- **W15** (the point G2): "elementary having two singular points, a saddle and a saddle-node."
- **W18** (the analytic curve G3): "an elementary graphic having two saddles."

**(c) Cyclicity of those graphics — Proposition 17 (VERBATIM):**

> "**Proposition 17** The following statement holds.
> (a) At most one limit cycle could arise near the graphic, from a quadratic system corresponding
> to a point on G1, in any perturbation of the system inside the class of all quadratic systems.
> Furthermore, if this limit cycle exists, it is hyperbolic.
> In addition there is numerical evidence for the following affirmation:
> (b) From the graphic corresponding to a quadratic system having its parameters on the curve
> G2 ∪ G3, at most one limit cycle could arise near the graphic, in any perturbation of the system
> inside the class of all quadratic systems. When this limit cycle exists it is hyperbolic."

Proof text (VERBATIM): "The graphics of the systems associated to these parameters (see W13) are
elementary (in the sense of [15]), having two singular points, a saddle and a saddle-node, which
satisfy the assumptions of Theorem 1 of [15]. **Therefore, the cyclicity of such a graphic is one**,
i.e., by a small perturbation, at most one limit cycle could bifurcate near the graphic."
(Reference [15] = Dumortier–Roussarie–Rousseau, JDE 110 (1994); Theorem 1 there is the cyclicity-1
criterion, further developed in [16] = *Elementary graphics of cyclicity 1 and 2*, Nonlinearity 7
(1994), 1001–1043.)

**The 4-limit-cycle ceiling in QW3 — Theorem 18 (VERBATIM):**

> "**Theorem 18** The following statements hold.
> (a) Inside the class of all quadratic systems there exists a neighborhood U1 of the systems in
> QW3 which have neither a limit cycle nor a graphic, i.e., whose phase portrait is Wi for some
> i = 1, ..., 12 such that any quadratic system in U1 has at most three limit cycles.
> (b) Inside the class of all quadratic systems there exists a neighborhood U2 of the family of
> systems in QW3 with phase portrait W13, such that any quadratic system in U2 has at most four
> limit cycles.
> In addition there is numerical evidence for the following two affirmations:
> (c) Inside the class of all quadratic systems, there exists a neighborhood U3 of the systems in
> QW3 with phase portraits Wi for i = 15, 18 such that any quadratic system in U3 has at most four
> limit cycles.
> (d) Inside the class of all quadratic systems, there exists a neighborhood U4 of the systems in
> QW3 with phase portraits Wi for i = 14, 16, 17, such that any quadratic system in U4 has at most
> four limit cycles."

And the stacking arithmetic, VERBATIM from the proof of Theorem 18:

> "A system with portrait W13 has two foci, the weak focus of third order and a strong focus which
> is surrounded by a graphic. This system has no limit cycles. By Proposition 17(a), if the
> perturbation of W13 is sufficiently small within the quadratic family, the only limit cycles which
> could be obtained are those produced by the weak focus plus the unique one produced by the graphic.
> So, the perturbed system can have at most four limit cycles."

> "**Corollary 19** Inside the class of all quadratic systems there exists a neighborhood of the
> class of quadratic systems with a weak focus of third order and without any polycycle or limit
> cycle such that any quadratic system in this neighborhood has at most 4 limit cycles."

### 1.3 The background theorems LS2004 relies on (§7.1, VERBATIM)

> "(iv) If a quadratic system has a limit cycle, then it surrounds a unique singular point, and
> this point is a focus; see [11]."
> "(v) If in a quadratic system the separatrix of an infinite saddle connects with the separatrix
> of the diametrically opposite infinite saddle, then this separatrix is an invariant straight line;
> see [49]." ([49] = Sotomayor–Paterlini 1983.)
> "(vi) If a quadratic system has a center, then it is integrable ..."
> "Finally we state here an important result on general quadratic systems obtained by Zhang
> Pingguang [52], [53]. **(vii) If there are limit cycles surrounding two foci of a quadratic system,
> then around one of the foci there is at most one limit cycle.**"
> "(viii) There are no limit cycles of systems (6) surrounding the weak focus of third order,
> see [24]."

> "[52] Zhang Pingguang, *On the distribution and number of limit cycles for quadratic systems with
> two foci.* (chinese), Acta Math. Sinica **44**(2001), 37–44.
> [53] ——, *On the distribution and number of limit cycles for quadratic systems with two foci.*
> Qualitative Theory [of Dynamical] Systems **3**(2002), 1–28."

### 1.4 Status of the DRR program as reported by LS2004 (VERBATIM, §1)

> "If we consider the family of all quadratic vector fields, it is not even proved that H(2) is a
> finite number. There is a program under way (cf. [15]) to prove this. The program proposes to show
> that all the 121 graphics listed in [15], which intervene in the problem of proving the finiteness
> of H(2), have finite cyclicity (cf. [36]). Because in over one hundred years since Hilbert stated
> the problem, no example was found where we can prove that more than four limit cycles exist, it is
> not only conjectured that H(2) is finite but also that H(2) = 4. All known examples of quadratic
> vector fields having at least four limit cycles were obtained by perturbing a quadratic vector
> field with a weak focus of third order, a good motivation for understanding this family."

And §10 (VERBATIM), on where the program is stuck:

> "It is thus not surprising that the program outlined in [15] and pursued in numerous articles has
> stopped short (so far) of proving the finite cyclicity of graphics present in these more degenerate
> systems, such as for example the system corresponding to the point Σ = [0:1:0] in Figure 1."

Also VERBATIM (§10(I)), the three known 4-limit-cycle mechanisms:

> "The first known examples of quadratic systems for which it was possible to show that they have at
> least four limit cycles are the example of Shi Songling [44] and that of Chen and Wang [9]. Both
> these examples are produced by perturbing systems S in the region of QW3 determined by
> J(S) = (7, 2, 2, 1), i.e., in R9 of Figures 2 and 3. In the early 1980s it became clear from the
> work [3] that there are two other distinct subsets of QW3 from which two new types of phase
> portraits with at least four limit cycles could be obtained by perturbations of systems in the
> class QW3. ... These are the systems in R8^dd and on Σ3."

---

## 2. Artés & Llibre 1997 — the QW3 precursor (open access, already in `papers/`)

**Citation.** J. C. Artés and J. Llibre, *Quadratic vector fields with a weak focus of third order*,
Publicacions Matemàtiques **41** (1997), 7–39.
**PDF:** `/Users/scottg/Claude_all/papers/pubmat41.pdf` (33 pp.) — this is where the labels
**W1 … W20** were introduced (LS2004 later merged W6=W9 and W13=W14 topologically, giving 18).

Abstract (VERBATIM): "We study phase portraits of quadratic vector fields with a weak focus of third
order at the origin. We show numerically the existence of at least 20 different global phase
portraits for such vector fields coming from exactly 16 different local phase portraits available for
these vector fields. Among these 20 phase portraits, 17 have no limit cycles and three have at least
one limit cycle."

Definition of graphic used, and the graphic-death mechanism (VERBATIM, §4):

> "As usual a *graphic* of p(X) ∈ P_n(S²) denotes a closed simple curve which is the union of a
> finite number of critical points and separatrices, and at least in one of the two sides of the
> curve is defined a return or Poincaré map."

> "We know that when we arrive to the curve C13 this limit cycle has disappeared because we cannot
> have a limit cycle around a node. As the point (0,1) does not change from unstable to stable
> anywhere close to where we are now, the limit cycle cannot die at this point and, so, **it has to
> die in an infinite graphic** in the same way that it has appeared."

The two classical constraints, quoted VERBATIM by Artés–Llibre:

> "[Yu]: If a quadratic vector field has a limit cycle, then it surrounds one and only one critical
> point, and this point is a focus."
> "[Li2]: There are no limit cycles around a weak focus of third order."
> "[SP]: If in a quadratic vector field the separatrix of an infinite saddle connects with the
> symmetric one, then they form an invariant straight line."

Conjecture (VERBATIM, §4):

> "So, we conjecture that quadratic vector fields having a weak focus of third order can not have
> limit cycles anywhere except the corresponding ones to the regions R14, C16 and R15, and that those
> only have one limit cycle."

And the standard expectation (VERBATIM, §1):

> "A quadratic vector field has at most 4 limit cycles, and when they exist, three are surrounding one
> focus and the other is surrounding a different focus. All the analytic examples realizing these four
> limit cycles are obtained perturbing a quadratic system having a weak focus of third order."

---

## 3. Artés–Llibre–Schlomiuk 2006 — weak focus of second order (QW2)

**Citation.** J. C. Artés, J. Llibre and D. Schlomiuk, *The geometry of quadratic differential
systems with a weak focus of second order*, Internat. J. Bifur. Chaos **16** (2006), no. 11,
3127–3194. DOI 10.1142/S0218127406016720.

**PDF: PAYWALLED — not obtained.** Exhaustive negative search record:
- UAB preprint server `https://mat.uab.cat/dpt/Publ/prep/` — **the archive only starts at 2008**
  (`p*_08.pdf` … `p*_25.pdf`); there are no `pNN_05.pdf` / `pNN_06.pdf` files at all. I downloaded and
  text-scanned the *entire* 2008–2012 run (167 PDFs, first page each): no weak-focus classification
  preprint is there.
- `ddd.uab.cat` (UAB digital repository): searched "weak focus quadratic", `"weak focus of second
  order"` — the only QW hits are the 1997 Publ. Mat. paper record and unrelated papers. No QW2 file.
- OpenAlex on DOI 10.1142/S0218127406016720: `{'is_oa': False, 'oa_status': 'closed',
  'any_repository_has_fulltext': False}`.
- Semantic Scholar `openAccessPdf`: `{"url": "", "status": "CLOSED"}`.
- `worldscientific.com` returns HTTP 403.
- Artés's own article page for this paper (`http://mat.uab.cat/~artes/articles/qvfwf2o.html`)
  hosts **only the supplementary Mathematica/CorelDraw/JPG slice files**, and says explicitly:
  *"We do not add here the paper due to the copyright, but we place only the extra files."*
  (I did not download the .zip bundles — they are 3-D bifurcation-surface pictures and slice
  drawings, ~hundreds of MB unzipped. They are available if you want the actual portrait plates:
  `3d1.zip`–`3d4.zip`, `qvf2owfmath.zip`, `slicesppcdr.zip`, `slicesnoppcdr.zip`,
  `slicesnoppjpg.zip`, `slicesexcel.zip`, all under `http://mat.uab.cat/~artes/articles/`.)

**What I could get VERBATIM — the authors' own abstract, from their web page**
(`http://mat.uab.cat/~artes/articles/qvfwf2o.html`):

> "In this paper we classify all the quadratic vector fields which are in the closure (within the
> quadratic family) of the family of systems with a weak focus of second order. This family includes
> apart from systems with a weak focus of second order, those with a weak focus of third order, and
> some with a center. The bifurcation diagram for this class, done in the adequate parameter space
> which is the 3-dimensional real projective space, is quite rich in its complexity and yields 373
> subsets with 126 phase portraits for the whole class, 95 of them specifically for quadratic vector
> fields with a weak focus of second order."

**The key question — "does any phase portrait in the closure of QW2 have a graphic (finite loop or
through infinity) surrounding the weak focus of order two?" — is NOT settled by any open source I
could reach.** What I can say from open sources:

- *Indirect but strong evidence for NO, in the "two finite singularities" case.* Kooij–Zegeling
  1998 (§6 below), Theorem 1.3 VERBATIM: "If the weak focus is of second or third order then the
  strong focus is surrounded by at most one limit cycle, and if it exists it is hyperbolic." That
  constrains the *strong* focus, not the weak one. The constraint on the weak focus of order 2
  itself is Bautin's bound (at most 2 small-amplitude cycles from an order-2 weak focus) plus
  Zhang Pingguang's (1,n) theorem.
- A widely-repeated secondary statement — *(paraphrase, UNVERIFIED, do not quote)* — is that a system
  in QW2 has at most one limit cycle surrounding the weak focus, hyperbolic when it exists, and at
  most two limit cycles in the whole plane. **I could not verify this against the paper.** Treat as
  a hypothesis to check; if true it kills the QW2 recipe outright.
- Independent contemporaneous testimony that the authors themselves think QW2 is where the record
  lives — BIRS 2007 workshop report (Rousseau, ed.), VERBATIM:
  > "They now have a complete bifurcation diagrams of quadratic systems with a weak focus of order
  > greater or equal to 2 [1] and they start attacking the case of a weak focus of order 1. **Their
  > study of the quadratic systems with weak order of order 2 has enlightened some corners of the
  > bifurcation diagram which they conjecture will produce the larger number of limit cycles in the
  > family of quadratic vector fields.** For their work they mix algebro-geometric techniques with
  > numerical simulations and Joan Carlos Artés lectured on the numerous phase portraits expected for
  > a system with a weak focus of order 1."
  (`/Users/scottg/Claude_all/papers/birs2007_hilbert16_workshop_report.pdf`, §3.6.)

**Recommendation:** this single paper is the highest-value paywalled item in slice B. Getting
IJBC 16(11) pp. 3127–3194 through an institutional library (or asking Artés directly — his page
invites contact at `artes@mat.uab.es`) would settle the QW2 question directly.

---

## 4. Artés–Llibre–Schlomiuk 2010 — weak focus + invariant straight line (QW^I)

**Citation.** J. C. Artés, J. Llibre and D. Schlomiuk, *The geometry of quadratic polynomial
differential systems with a weak focus and an invariant straight line*, Internat. J. Bifur. Chaos
**20** (2010), no. 11, 3627–3662. DOI 10.1142/S021812741002791X.

**PDF: PAYWALLED.** OpenAlex: `{'is_oa': False, 'oa_status': 'closed'}`. Not on ddd.uab.cat, not on
the UAB preprint server (checked 2008–2012 exhaustively). Author page
`http://mat.uab.cat/~artes/articles/qvfwfsl.html` holds only supplementary zips
(`qvfwfsl/qvfwfslmath.zip`, `qvfwfsl/slicesppcdr.zip`, `qvfwfsl/slicesnoppcdr.zip`).

**Abstract (VERBATIM, from the authors' own page):**

> "In this article we make a global study of the class $QW^I$ of all real quadratic differential
> systems which have a weak focus and invariant straight lines of total multiplicity at least two.
> **As it turns out, these conditions imply that all the systems in $QW^I$, having a weak focus at
> some point necessarily have a weak focus of order one at that point.** In the closure of this
> family we find some systems having a center. The bifurcation diagram for this class, done in the
> adequate parameter space which is the 3-dimensional real projective space, is quite rich in its
> complexity and yields 151 subsets with 99 phase portraits for the whole class, 73 of them
> specifically for quadratic vector fields with a weak focus of first order."

**Relevance to stacking.** An invariant straight line forces the weak focus down to order 1, so the
"weak focus order k >= 2 + invariant line + graphic" recipe is dead by construction inside QW^I.
Note also LS2004 §7.1(v): an infinite-saddle-to-opposite-infinite-saddle separatrix connection in a
quadratic system *is* an invariant straight line — i.e. many "graphics through infinity" force an
invariant line, which then forces the weak focus to order 1.

---

## 5. QW1 — weak focus of first order (the live frontier)

All of these are from the Artés school; abstracts VERBATIM from `http://mat.uab.cat/~artes/articles/`.

### 5.1 Artés & Trullàs 2026 — QW1 + a finite saddle-node
**Citation.** J. C. Artés and C. Trullàs, *Quadratic differential systems with a weak focus of first
order and a finite saddle-node*, Internat. J. Bifur. Chaos **36** (2026), no. 6.
**PDF: PAYWALLED / not yet posted.** Supplement page: `http://mat.uab.cat/~artes/articles/qwf1sn/qwf1sn.html`.
Abstract (VERBATIM, excerpt):

> "In this article we perform a global study (modulo islands) of the class $\overline{Qwf1sn}$ which
> is the closure within real quadratic differential systems, of the family Qwf1sn of all such systems
> which have a weak focus of first order and a finite saddle-node. The bifurcation diagram for this
> class, done in the adequate parameter space which is the 3-dimensional real projective space RP^3,
> is quite rich in its complexity since yields 399 subsets with 192 topologically distinct phase
> portraits for the closure of Qwf1sn, 146 of which have a representative in Qwf1sn. **It can be shown
> that some of these parts have at least two limit cycles.** ... The bifurcation set is formed by an
> algebraic set of bifurcations of singularities, finite or infinite and by a set of bifurcations which
> we suspect to be analytic corresponding to **global separatrices which have connections, or double
> limit cycles**."

### 5.2 Artés & Cairó — QW1 + a (1,1) infinite saddle-node
**Citation.** J. C. Artés and L. Cairó, *Phase portraits of quadratic differential systems with a weak
focus and a $\binom{1}{1}SN$* (Qwf1SN11). **PDF: not posted** (supplement page
`http://mat.uab.cat/~artes/articles/QWF1SN11/QWF1SN11.html`).
Abstract (VERBATIM, excerpt): "... gives 303 subsets with 165 phase portraits for
$\overline{\textbf{Qwf1SN11}}$, 120 for **Qwf1SN11**, **28 having limit cycles** and 13 with center.
... The bifurcation set is formed by an algebraic set of finite or infinite singularity bifurcations
and by a non-algebraic set of points **corresponding to global separators that have connections or
double limit cycles.**"

### 5.3 Artés & Cairó — QW1 + a (0,2) infinite saddle-node
**Citation.** J. C. Artés and L. Cairó, *Phase portraits of quadratic differential systems with a weak
focus and a $\binom{0}{2}SN$* (Qwf1SN02). **PDF: not posted** (`.../QWF1SN02/QWF1SN02.html`).
Abstract (VERBATIM, excerpt): "... gives 295 subsets with 159 phase portraits for
$\overline{\textbf{Qwf1SN02}}$, **25 having limit cycles** and 10 with centers. ... a non-algebraic
set of points corresponding to global separators that have **separatrix connections or double limit
cycles**."

**Reading.** Every one of these QW1 families reports a non-algebraic bifurcation set consisting of
*separatrix connections* (i.e. graphics) and *double limit cycles*. None of the public abstracts
reports more than 2 limit cycles anywhere. There is no public claim of a QW1 portrait with a graphic
surrounding the weak focus itself.

---

## 6. Zhang Pingguang and the Chinese uniqueness/distribution results

### 6.1 Zhang Pingguang & Cai Suilin 1991 — OBTAINED
**Citation.** Zhang Pingguang and Cai Suilin, *Quadratic systems with a weak focus*, Bull. Austral.
Math. Soc. **44** (1991), 511–526. DOI 10.1017/S0004972700030008.
**PDF:** `/Users/scottg/Claude_all/papers/zhang_pingguang_cai_suilin_1991_quadratic_systems_weak_focus_BAMS44.pdf`
(16 pp., free on Cambridge Core).

Abstract (VERBATIM):
> "In this paper we study the number and the relative position of the limit cycles of a plane
> quadratic system with a weak focus. In particular, we prove **the limit cycles of such a system can
> never have (2,2)-distribution**, and that there is at most one limit cycle not surrounding this weak
> focus under any one of the following conditions: (i) the system has at least 2 saddles in the finite
> plane, (ii) the system has more than 2 finite singular points and more than 1 singular point at
> infinity, (iii) the system has exactly 2 finite singular points, more than 1 singular point at
> infinity, and the weak focus is itself surrounded by at least one limit cycle."

Theorems (VERBATIM):
> "**THEOREM A.** Suppose a quadratic system has 4 finite singular points: a strong focus, a weak
> focus, and two saddles. Then the system has at most one limit cycle surrounding the strong focus."
> "**THEOREM 3.** Suppose (1) satisfies (11) and Δ₀ ≥ 0, Δ ≥ 0. Then (1) has at most one limit cycle
> surrounding the strong focus O."
> "**THEOREM 4.** Suppose (1) satisfies (11) and Δ₀ ≥ 0, Δ < 0. Then, if there is a limit cycle
> surrounding the weak focus, there is at most one limit cycle surrounding the strong focus."
> "The limit cycles of a quadratic system are said to have (i, j)-distribution if there are exactly i
> limit cycles surrounding one focus and exactly j limit cycles surrounding another focus (multiple
> limit cycles being counted according to their multiplicity)."
> "**THEOREM 5.** If a quadratic system has a weak focus and a strong focus, its limit cycles cannot
> have (2,2) distribution."
> "**THEOREM B.** Suppose the finite singular points of a quadratic system include a strong focus, a
> weak focus, and at most one saddle. Suppose also that the system is positively normalised. Then if
> the strong focus is unstable, it has no limit cycle surrounding it."

### 6.2 Zhang Pingguang 2001 / 2002 — two-foci distribution
**Citations.**
- Zhang Pingguang, *On the distribution and number of limit cycles for quadratic systems with two
  foci*, (Chinese) Acta Math. Sinica **44** (2001), 37–44.
- Zhang Pingguang, *On the distribution and number of limit cycles for quadratic systems with two
  foci*, Qual. Theory Dyn. Syst. **3** (2002), 1–28. DOI 10.1007/BF02969414.
**PDF: PAYWALLED.** OpenAlex on the DOI: `{'is_oa': False, 'oa_status': 'closed'}`. No open Chinese
version located (CNKI / SciChina return no free full text; I did not attempt to bypass any paywall
or captcha). The English version is the QTDS one and is the one LS2004 recommends
("For this see the English version [53] of [52]").

**Content, VERBATIM as quoted by Llibre–Schlomiuk 2004 §7.1(vii):**
> "(vii) If there are limit cycles surrounding two foci of a quadratic system, then around one of the
> foci there is at most one limit cycle."
And LS2004 §1, VERBATIM: "a recently obtained result affirming that if a quadratic system has limit
cycles around two foci, then around one of them we cannot have more than just one limit cycle."

LS2004 also *uses* it to rule out a semistable cycle (VERBATIM, §8): "By Zhang Pingguang's result
[52, 53], crossing Θ from a node to a focus there could appear at most one limit cycle surrounding
the focus, which must be semistable near the bifurcation curve ... But if for a system X in QW3 such
a semistable limit cycle exists, then close to it there would be two limit cycles surrounding the
focus, one stable and the other unstable, ... and this is in contradiction with Zhang Pingguang's
result."

### 6.3 Zhang Pingguang 1999 — weak focus + strong focus
**Citation.** Zhang Pingguang, *Quadratic systems with a weak focus and a strong focus*, Appl. Math.
Chinese Univ. (Gaoxiao Yingyong Shuxue Xuebao) **14** (1999), 7–14. **PDF: not located, no open text.**
*(paraphrase, secondary — from search summaries only:)* proves that if a quadratic system with a weak
focus and a strong focus has limit cycles around both foci simultaneously, there is a unique limit
cycle around one of them. This is the precursor of the 2001/2002 theorem.

### 6.4 Zhang Pingguang — unbounded separatrix cycles
**Citation.** Zhang Pingguang, *Unbounded separatrix cycles in quadratic systems*, (Chinese), Proc.
Conf. Diff. Eq., Nanking, October 1993, 11 pp. — listed as item [* P-107] in J. W. Reyn's quadratic
systems preprint bibliography (`http://ta.twi.tudelft.nl/dv/staff/reyn/biblio/preprints.html`).
**PDF: not located.** Flagged because "unbounded separatrix cycle" = graphic through infinity, i.e.
exactly the object of interest; worth chasing if a Chinese-language library is available.

---

## 7. Kooij & Zegeling 1998 — the sharpest open statement about a *second-order* weak focus

**Citation.** R. E. Kooij and A. Zegeling, *Limit cycles in quadratic systems with a weak focus and a
strong focus*, Kyungpook Math. J. **38** (1998), no. 2, 323–340.
**PDF: PARTIAL — OBTAINED (open, first ~8 pp.).**
`/Users/scottg/Claude_all/papers/kooij_zegeling_1998_weak_focus_strong_focus_Kyungpook38.pdf`
(from `https://kmj.knu.ac.kr/journal/download_pdf.php?spage=323&volume=38&number=2`; the server
returned only the opening pages, but they contain all the statements below).

**Coppel's Conjecture W (VERBATIM):**
> "**Conjecture W.** If a quadratic system contains a weak singularity, i.e. a point (x0, y0) where
> P2(x0,y0) = Q2(x0,y0) = 0 and div(P2(x0,y0), Q2(x0,y0)) = 0, then there exists at most one limit
> cycle not surrounding this singularity. If such a limit cycle exists, it is hyperbolic.
> Conjecture W is still open, although partial results have been obtained."

**Main theorem (VERBATIM):**
> "**Theorem 1.3.** Consider a quadratic system with exactly two finite real singularities, a weak
> focus and a strong focus, and at least two singularities at infinity. **If the weak focus is of
> second or third order then the strong focus is surrounded by at most one limit cycle, and if it
> exists it is hyperbolic.** If the weak focus is of first order while the two foci have opposite
> stability, then the same conclusion holds."

**Two other quoted results (VERBATIM):**
> "**Theorem 1.1 (Zhang and Cai [17])** A quadratic system with exactly two finite real singularities,
> a weak focus surrounded by at least one limit cycle and a strong focus, and at least two
> singularities at infinity, has at most one limit cycle surrounding the strong focus. If it exists it
> is hyperbolic."
> "**Theorem 1.2 (Zegeling and Kooij [15])** A quadratic system with three or four finite real
> singularities, including a weak singularity and a strong focus, has at most one limit cycle
> surrounding the strong focus. If it exists then it is hyperbolic and furthermore this limit cycle is
> unique in the whole phase plane."

**And two structural facts (VERBATIM, §2):**
> "In a quadratic system a limit cycle has to surround a focus, see [12]. Furthermore, **a quadratic
> system with two weak foci has no limit cycles**, see [13]."
> "It is well-known that the order of a weak focus in a quadratic system is at most three, see [3] or [12]."

**Remark 1.2 (VERBATIM) — the one uncovered case:**
> "The case where a quadratic system has at least two singularities at infinity, a first order weak
> focus not surrounded by limit cycles and a strong focus while the two foci have the same stability,
> is neither covered by Theorem 1.1, nor by Theorem 1.3."

---

## 8. The global-classification programme (ALSV) — Birkhäuser 2021 and the open 2020 preprint

### 8.1 ALSV 2020 — OBTAINED (open preprint on ddd.uab.cat)
**Citation.** J. C. Artés, J. Llibre, D. Schlomiuk, N. Vulpe, *Global topological configurations of
singularities for the whole family of quadratic differential systems*, Qual. Theory Dyn. Syst.
**19** (2020), Paper 51.
**PDF:** `/Users/scottg/Claude_all/papers/ALSV2020_global_topological_configurations_singularities_whole_quadratic_family_ddd.pdf`
(27 pp.) — from `https://ddd.uab.cat/pub/artpub/2020/221281/ArtLliSchVul2019_preprint.pdf`.

Abstract (VERBATIM, excerpt):
> "In [1] the authors proved that there are 1765 different global geometrical configurations of
> singularities of quadratic differential systems in the plane. There are other 8 configurations
> conjectured impossible ... We prove that there are exactly 208 topologically distinct global
> topological configurations of singularities for the whole quadratic class. ... **From here the next
> goal would be to obtain a bound for the number of possible different phase portraits, modulo limit
> cycles.**"

**Relevance:** this is a classification of *singularity configurations*, not of graphics. It gives no
statement about which graphics surround a focus. Its value here is negative/structural: the global
quadratic classification has, as of 2020–2026, only reached the level of singularity configurations;
the phase-portrait-level classification (which is where graphics live) is still being assembled
family by family (QW3, QW2, QW^I, QW1SN, QW1SN11, QW1SN02, saddle-node families, codimension-1 and -2
books, …).

### 8.2 ALSV 2021 book — NOT OBTAINED
**Citation.** J. C. Artés, J. Llibre, D. Schlomiuk, N. Vulpe, *Geometric Configurations of
Singularities of Planar Polynomial Differential Systems — A Global Classification in the Quadratic
Case*, Birkhäuser/Springer, Cham, 2021. ISBN 978-3-030-50569-1.
**Paywalled (Springer monograph).** Same caveat as 8.1: singularities, not graphics.

### 8.3 ALSV 2026 — codimension (open access, EJDE)
**Citation.** J. C. Artés, J. Llibre, D. Schlomiuk, N. Vulpe, *Codimension in planar polynomial
differential systems*, Electron. J. Differential Equations **28** (2026), 1–47. EJDE is open access —
not fetched in this slice, flagged as easy to obtain from `ejde.math.txstate.edu` if the codimension
stratification is wanted.

---

## 9. The DRR programme and the 121 graphics — status as of 2026

### 9.1 The source
**Citation.** F. Dumortier, R. Roussarie and C. Rousseau, *Hilbert's 16th problem for quadratic vector
fields*, J. Differential Equations **110** (1994), 86–133.
**PDF: PAYWALLED (Elsevier).** Companion paper: F. Dumortier, R. Roussarie, C. Rousseau, *Elementary
graphics of cyclicity 1 and 2*, Nonlinearity **7** (1994), 1001–1043 — also paywalled.

### 9.2 What the 121 graphics ARE — and the answer to "which can surround a focus"

**All 121 of them surround an anti-saddle.** VERBATIM, Rousseau–Shan–Zhu, arXiv:1502.00689, §1:

> "Also, limit cycles in quadratic vector fields necessarily surround a unique singular point with
> nondegenerate linear part, and linear vector fields can have no limit cycles. Hence, it is possible
> to compactify the space of equivalence classes of quadratic vector fields with a nondegenerate
> singular point of anti-saddle type: this yields a compact parameter space K. Limit cycles in the
> compact set S²×K accumulate on graphics ... **The DRR program reduces the proof that H(2) < ∞ to the
> proof that each graphic Γ ⊂ S² surrounding a nondegenerate singular point of anti-saddle type and
> occurring for a parameter value A0 ∈ K has finite cyclicity in S²×K**, i.e. can produce only a finite
> number of limit cycles in a neighborhood U of Γ for parameter values A in a neighborhood V of A0.
> **Achieving the DRR program requires proving the finite cyclicity of 121 graphics in S²×K.**"

So the question "which of the OPEN graphics can surround a focus?" has the answer: **every one of the
121, open or closed, surrounds an anti-saddle (focus or center) by construction.** The graphics are
partitioned instead by (i) whether the surrounded anti-saddle is a *center* (the "center graphics",
handled by the Bautin trick) or a generic focus, and (ii) the type of the singular points *on* the
graphic (hyperbolic / semi-hyperbolic saddle-node / nilpotent triple point / degenerate), which is
what the DRR labels encode: `H` (hyperbolic), `H^k_j` hemicycles, `I^k_j` (through a triple nilpotent
point at infinity), `DI`, `DF`, `F`, etc.

### 9.3 The count of proved graphics

**As of 2015 (VERBATIM, Rousseau–Shan–Zhu, arXiv:1502.00689):**
> "In this paper, we prove that the two graphics through a nilpotent point of saddle type, (I¹₁₂) and
> (I¹₁₃), that do not surround a center, have finite cyclicity. **Therefore the results from this paper
> will bring the number of graphics of the program for which finite cyclicity is proved to 88.**"

=> **88 proved, 33 open, as of 2015.** I could find **no published updated count** anywhere between
2015 and 2026; the literature since then reports individual graphics, not running totals.

**Roussarie–Rousseau 2015 (arXiv:1506.07104), VERBATIM:**
> "We have a partial result for every graphic, but one (namely (H³₁₄)), through a triple point at
> infinity:
> **Theorem 1.1.** Let us consider the graphics (I¹₁₄), (I¹₆b), (H³₁₃) and (DI₂b) through a triple point
> at infinity (see Figure 1). Then for any of them, the boundary periodic limit set obtained in the
> blowing up has a finite cyclicity."
> "**Theorem 1.2.** The graphic (I¹₁₄) has a finite cyclicity inside the family of quadratic vector
> fields."
> "As for the finite cyclicity of the other graphics (I¹₆b), (H³₁₃) and (DI₂b), we intend to address the
> problem in the next future. The finite cyclicity of (H³₁₃) should be straightforward with arguments
> identical to those used for (I¹₁₄). It will be done simultaneously with the corresponding generic
> graphic (H³₁₂). Some of the limit periodic sets to be studied for (I¹₆b) will involve four Dulac maps
> of second type. ... As for the graphic (DI₂b), some of the limit periodic sets to be studied involve
> four Dulac maps of second type, two of them through the semi-hyperbolic points P1 and P2 on the
> blown-up sphere."
> "We also hope to adapt them to study the boundary graphic of the hemicycle (H³₁₄): there, the
> additional difficulty is the two semi-hyperbolic points along the equator."

=> As of 2015 the *named* still-open hard cases at a triple point at infinity are
**(I¹₆b), (H³₁₃) [expected easy], (DI₂b), and the hemicycle (H³₁₄) [hardest]**; all of them surround a
**center** (they are in the reversible stratum) — VERBATIM: "all graphics through a nilpotent point and
surrounding a center occur in the stratum of reversible systems."

**2026 claim on (H³₁₄) — PREPRINT, UNREFEREED, TREAT WITH CAUTION.**
**Citation.** Haibo Lu, *Local Uniform Finite Cyclicity of the H³₁₄ Semihyperbolic Hemicycle in
Quadratic Systems*, arXiv:2607.13785v3 (26 Aug 2026), 99 pp.
**PDF: OBTAINED.** `/Users/scottg/Claude_all/papers/h3_14_hemicycle_finite_cyclicity_2607.13785.pdf`
Abstract (VERBATIM, excerpt):
> "We prove local uniform finite cyclicity of the labelled H³₁₄ hemicycle in a full neighborhood of its
> base field in the twelve-dimensional space of planar quadratic vector fields. ... The resulting bound
> is existential and is not claimed to be sharp. **In the terminology of the quadratic finite-cyclicity
> program, this completes the labelled H³₁₄ open case.**"
This is a single-author arXiv preprint with heavy computer-assisted content (Appendix F: "Computer-assisted
calculations and data availability"). It is **not** a peer-reviewed confirmation that H³₁₄ is done.

**Bottom line on the DRR status:** the last *citable* count is 88/121 (2015). Since then, individually:
(I¹₁₄) done (Roussarie–Rousseau 2015); (I¹₁₂), (I¹₁₃) done (Rousseau–Shan–Zhu 2015); (H³₁₃), (I¹₆b),
(DI₂b) partial (boundary limit periodic set only); (H³₁₄) claimed in an unrefereed 2026 preprint. So
roughly **~30 graphics remain genuinely open in 2026**, and I found no source claiming the programme is
close to complete. Nobody claims a *large* cyclicity for any of them; the open cases are open because
the *analysis* (four Dulac maps of second type, nilpotent blow-ups with 2-D displacement maps) is hard,
not because large numbers of cycles are suspected.

### 9.4 Supporting open PDFs obtained for §9
- `/Users/scottg/Claude_all/papers/finite_cyclicity_graphics_nilpotent_saddle_quadratic_1502.00689.pdf`
  — Rousseau, Shan, Zhu, *Finite cyclicity of some graphics through a nilpotent point of saddle type
  inside quadratic systems*, arXiv:1502.00689 (= QTDS 2015). 18 pp.
- `/Users/scottg/Claude_all/papers/finite_cyclicity_center_graphics_nilpotent_quadratic_1506.07104.pdf`
  — Roussarie & Rousseau, *Finite cyclicity of some center graphics through a nilpotent point inside
  quadratic systems*, arXiv:1506.07104. 42 pp.
- `/Users/scottg/Claude_all/papers/h3_14_hemicycle_finite_cyclicity_2607.13785.pdf` — Lu 2026 (see above).
- `/Users/scottg/Claude_all/papers/birs2007_hilbert16_workshop_report.pdf` — Rousseau (ed.),
  BIRS workshop report 07w5021, *Mathematical developments around Hilbert's 16th problem*, 2007.
  VERBATIM: "An important contribution is this direction is the program started in 1991 by
  Dumortier-Roussarie-Rousseau [7] and reducing the proof that H(2) < ∞ to the proof that 121 graphics
  have finite cyclicity."
- `/Users/scottg/Claude_all/papers/yakovenko_tangential_hilbert16_lectures_math0104140.pdf` —
  S. Yakovenko, *Quantitative theory of ODEs and tangential Hilbert 16th problem*, CRM Proc./
  arXiv:math/0104140 (background on cyclicity, not quadratic-specific).
- **PAYWALLED, not obtained:** Dumortier–El Morsalani–Rousseau, *Hilbert's 16th problem for quadratic
  vector fields and cyclicity of graphics*, Nonlinear Anal. 1997; *Finite cyclicity of elementary
  graphics surrounding a focus or center in quadratic systems*, QTDS 3 (2002) (DOI 10.1007/BF02969336 —
  OpenAlex: closed); Dumortier–Rousseau, *Study of the cyclicity of some degenerate graphics inside
  quadratic systems*, CPAA 8 (2009); Rousseau–Zhu, *Finite cyclicity of HH-graphics …*.

---

## 10. Graphics through infinity surrounding a focus/center — the sharp cyclicity numbers

**Citation.** D. Marín and J. Villadelprat, *The cyclicity of hyperbolic hemicycles*, arXiv:2501.16924
(28 Jan 2025), 48 pp.
**PDF: OBTAINED** (was already in `papers/`): `/Users/scottg/Claude_all/papers/marin_villadelprat_2025_2501.16924.pdf`

This is the most directly relevant *quantitative* paper for the "graphic through infinity surrounding a
focus" question.

Abstract (VERBATIM, excerpt):
> "We consider families of planar polynomial vector fields of degree n and study the cyclicity of a type
> of unbounded polycycle Γ called **hemicycle**. Compactified to the Poincaré disc, Γ consists of an
> affine straight line together with half of the line at infinity and has two singular points, which are
> hyperbolic saddles located at infinity. ... For the other results we consider the case n = 2. More
> concretely they are addressed to the quadratic integrable systems belonging to the class Q^R₃ and
> having two hemicycles, Γu and Γℓ, surrounding each one a center."

> "**Theorem B.** Let us take any (a0, b0) ∈ (−2, 0)×(0, 2). Then **the cyclicity of Γu when we perturb
> (7) inside the whole family of quadratic differential systems is exactly 2** if a0 ≠ −1 and at least 2
> if a0 = −1. Moreover the same statement is true for Γℓ."

> "**Theorem C.** If (a0, b0) belongs to K1\{a0 = −1} (respectively, K2) then **the cyclicity of
> Π = {Γu, Γℓ} when we perturb (7) inside the whole family of quadratic differential systems is exactly
> 3** (respectively, 2). Moreover it is at least 3 for (a0, b0) ∈ {−1}×(0, 2)."

> "We stress that Theorem C deals with the simultaneous bifurcation of limit cycles from Γu and Γℓ,
> which are the outer boundaries of two period annuli. Note that if this simultaneous cyclicity is 3
> then, as a consequence of Theorem B, two limit cycles bifurcate from Γu and one from Γℓ, or vice versa."

And, on the strongest known single-polycycle cyclicity in the quadratic family (VERBATIM):
> "a result due to Swirszcz (see [32, Theorem 1]) ... he perturbs the differential system (7) but taking
> (a0, b0) ∈ S := {0 < b0 < −a0} ∩ {a0 < −2}. For these parameters the singular point (0, 1/2) is also a
> center but the polycycle at the boundary of its period annulus is not an hemicycle. It is a **bicycle**
> Γb with the two vertices at infinity, and consisting of a branch of a hyperbola together with a segment
> of ℓ∞. ... Swirszcz identifies a curve C (see Figure 4) such that **the cyclicity of Γb is 3 if
> (a0, b0) ∈ S∩C and 2 if (a0, b0) ∈ S\C**."

And on the period annulus (Iliev), VERBATIM:
> "Among other results Iliev proves that the cyclicity of the period annulus P of the center at (0, 1/2)
> of the differential system (7) is 3 for (a0, b0) = (−4, 2) and 2 for (a0, b0) = (−1, 1) ... Moreover he
> conjectures that the cyclicity of P is equal to 3 if (a0, b0) is inside the shaded area in Figure 4 and
> equal to 2 if (a0, b0) is outside."

**This is the key number for the stacking arithmetic: a graphic through infinity surrounding a
center/focus in the quadratic family has cyclicity 2 (generic hemicycle) or 3 (the Swirszcz bicycle on a
codimension-1 curve). Not more, in any case that has been computed.** And when two such polycycles are
present simultaneously, the *joint* cyclicity is 3, not 4.

---

## 11. Gasull & Santana 2025 — status of H(2)

**Citation.** A. Gasull and P. Santana, *A note on Hilbert 16th problem*, arXiv:2407.13465v2 (1 Oct 2024);
published in Proc. Amer. Math. Soc. (2025).
**PDF: OBTAINED.** `/Users/scottg/Claude_all/papers/gasull_santana_note_hilbert16_2407.13465.pdf` (9 pp.)

Status statements (VERBATIM):
> "Under this notation the second part of Hilbert's 16th problem consists in obtaining an upper bound for
> H(n) and it is yet an open problem. **Even for the quadratic case, it is not known if H(2) < ∞.**
> However, advances has been made and lower bounds for H(n) have been found. For small values of n, the
> best lower bounds so far are **H(2) ≥ 4** [3, 18], H(3) ≥ 13 [10] and H(4) ≥ 28 [15]."

Main results (VERBATIM):
> "**Theorem 1.** Given n ∈ N, it holds H(n + 1) ≥ H(n) + 1."
> "**Theorem 2.** For n ∈ N, the following statements hold. (a) If H(n) < ∞, then there is X ∈ Σⁿ_h such
> that π(X) = H(n). (b) If H(n) = ∞, then for each k ∈ N there is X_k ∈ Σⁿ_h such that π(X_k) ≥ k."

> "In particular, it follows from Theorem 1 that if H(n0) = ∞ for some n0 ∈ N, then H(n) = ∞ for every
> n ≥ n0."

**Load-bearing consequence for the search programme:** Theorem 2(a) says that if H(2) is finite, the
maximum is realised by a **structurally stable** quadratic field with **only hyperbolic** limit cycles.
So a hypothetical 5-limit-cycle quadratic system may be assumed structurally stable with 5 hyperbolic
cycles — which makes numerical/interval-arithmetic search legitimate and complete in principle (no need
to hunt semistable or degenerate configurations to find the *maximum*).

---

## 12. Villanueva & Tucker 2026

**Citation.** Y. Villanueva and W. Tucker, *Darboux-type center conditions for families of planar
polynomial vector fields*, arXiv:2602.22558v2 (2 Jul 2026), 6 pp. (dated July 7, 2026).
**PDF: OBTAINED** (already in `papers/`): `/Users/scottg/Claude_all/papers/villanueva_tucker_2026_2602.22558.pdf`

This paper is about the **center–focus problem / Bautin ideal**, not about graphics or about H(2)
directly. What it says about the status of H(2) (VERBATIM):

> "The second part of Hilbert's 16th problem asks for the Hilbert number H(n) – the maximal number of
> limit cycles (isolated periodic orbits) the family of planar polynomial ordinary differential equations
> of degree n can display. Note that H(n) should only depend on the degree n, not on the particular
> polynomial vector field itself. **This question is unresolved, even for the simplest case n = 2.** Even
> finding non-trivial lower bounds for H(n) appears to be very hard."

> "A classical result of Bautin [1, 15] states that for the class of quadratic vector fields, **the maximal
> number of small amplitude limit cycles is three: M(2) = 3.**"

Main result (VERBATIM): "For F(n) – the class of (non-homogeneous) differential systems of degree n – we
have: B(F(n)) = ⟨L1, L2, ...⟩ ⊂ ⟨v^h₃, v^h₄, ..., v^h_{n+1}⟩."

**Relevance:** it contributes nothing about graphics; its contribution to the stacking question is the
restatement of **M(2) = 3** — i.e. the local (small-amplitude) contribution from any single quadratic
weak focus is capped at 3, no matter what the order is (order 3 being the maximum weakness).

---

## 13. Other files collected in this slice

| File in `papers/` | What it is |
|---|---|
| `llibre_schlomiuk_2004_CJM56_weak_focus_third_order.pdf` | LS2004 QW3, full text, §1 |
| `pubmat41.pdf` | Artés–Llibre 1997 QW3 precursor, W1–W20, §2 |
| `zhang_pingguang_cai_suilin_1991_quadratic_systems_weak_focus_BAMS44.pdf` | (2,2)-distribution impossible, §6.1 |
| `kooij_zegeling_1998_weak_focus_strong_focus_Kyungpook38.pdf` | Conjecture W, order-2/3 uniqueness, §7 (partial scan) |
| `ALSV2020_global_topological_configurations_singularities_whole_quadratic_family_ddd.pdf` | 208 topological configurations, §8.1 |
| `finite_cyclicity_graphics_nilpotent_saddle_quadratic_1502.00689.pdf` | 88/121 count, DRR framing, §9 |
| `finite_cyclicity_center_graphics_nilpotent_quadratic_1506.07104.pdf` | (I¹₁₄) done; (I¹₆b),(H³₁₃),(DI₂b),(H³₁₄) open, §9 |
| `h3_14_hemicycle_finite_cyclicity_2607.13785.pdf` | 2026 unrefereed claim on H³₁₄, §9 |
| `birs2007_hilbert16_workshop_report.pdf` | QW2 "largest number of limit cycles" conjecture, §3 |
| `marin_villadelprat_2025_2501.16924.pdf` | hemicycle cyclicity = 2, bicycle = 3, joint = 3, §10 |
| `gasull_santana_note_hilbert16_2407.13465.pdf` | H(2) finiteness open; max realised structurally stably, §11 |
| `villanueva_tucker_2026_2602.22558.pdf` | Bautin ideal; M(2)=3, §12 |
| `yakovenko_tangential_hilbert16_lectures_math0104140.pdf` | cyclicity background |

**Confirmed unavailable open-access (all checked via OpenAlex/Semantic Scholar/ddd.uab.cat/mat.uab.cat/
worldscientific/Cambridge; no paywall or captcha bypass attempted):** ALS 2006 (IJBC 16, QW2);
ALS 2010 (IJBC 20, QW^I); Artés–Trullàs 2026 (IJBC 36); Artés–Cairó QW1SN11 and QW1SN02; Zhang Pingguang
QTDS 3 (2002) and Acta Math. Sinica 44 (2001); Zhang Pingguang, Appl. Math. Chinese Univ. 14 (1999);
DRR 1994 (JDE 110) and DRR 1994 (Nonlinearity 7); Dumortier–El Morsalani–Rousseau QTDS 3 (2002)
"Finite cyclicity of elementary graphics surrounding a focus or center"; ALSV Birkhäuser 2021.

---

## 14. WHICH STACKING RECIPES (weak focus of order k + surrounding graphic of cyclicity m) ARE ALIVE, DEAD, OR UNKNOWN

Notation: the recipe is "perturb a quadratic system that has a weak focus of order k at one point and a
graphic of cyclicity m surrounding *some* anti-saddle, and collect k + m limit cycles". Target: k + m ≥ 5.

### DEAD

**D1. k = 3, graphic surrounding the weak focus of order 3 itself. — DEAD.**
No such configuration exists: LS2004 Theorem 16(III) says the unique graphic in QW3 surrounds the
**strong** focus, and Li Chengzhi (Chinese Ann. Math. B 7 (1986), 174–190; quoted VERBATIM in LS2004
§7.1(viii) and in Artés–Llibre 1997) says there are no limit cycles around a weak focus of order 3 at all.
*Source: LS2004 §8 Thm 16; LS2004 §7.1(viii); Artés–Llibre 1997 §4.*

**D2. k = 3, graphic (through infinity) surrounding the strong focus, m ≥ 2. — DEAD (m = 1).**
LS2004 Proposition 17: the QW3 graphics (W13 on G1 — saddle + saddle-node; W15 at G2 — saddle +
saddle-node; W18 on G3 — two saddles) are **elementary** and have **cyclicity exactly 1**, with the
bifurcating cycle hyperbolic. Part (a) is analytic (via DRR 1994 Theorem 1); parts (b) for W15, W18 are
numerically verified only.
Total: 3 + 1 = 4. *Source: LS2004 Prop. 17, Thm 18, Cor. 19.*

**D3. Any (2,2) or richer two-focus distribution. — DEAD.**
Zhang Pingguang (Acta Math. Sinica 44 (2001) / QTDS 3 (2002)): with limit cycles around two foci, one of
them carries at most one. So the only distributions are (1, n). To reach 5 you need (1,4) — i.e. 4 cycles
around a single focus — or (0,5).
*Source: LS2004 §7.1(vii) VERBATIM; Zhang–Cai 1991 Thm 5 for the weak-focus special case.*

**D4. Two weak foci. — DEAD.**
"a quadratic system with two weak foci has no limit cycles". *Source: Kooij–Zegeling 1998 §2 VERBATIM.*

**D5. k = 2 or 3 + strong focus with ≥ 2 cycles (two finite singularities). — DEAD.**
Kooij–Zegeling Theorem 1.3: with exactly two finite real singularities and ≥ 2 at infinity, if the weak
focus is of order 2 or 3, the strong focus carries **at most one** limit cycle (hyperbolic). Since
Bautin caps the weak focus at k, the total is ≤ k + 1 ≤ 4.
*Source: Kooij–Zegeling 1998 Thm 1.3 VERBATIM.*

**D6. k = 2 or 3 + weak focus with ≥ 3 finite singular points. — DEAD.**
Zegeling–Kooij (quoted as Thm 1.2 in Kooij–Zegeling 1998): with three or four finite real singularities
including a weak singularity and a strong focus, the strong focus carries at most one limit cycle, and it
is unique in the whole phase plane. Also Zhang–Cai 1991 Theorem A for the 4-singularity case.
*Source: Kooij–Zegeling 1998 Thm 1.2 VERBATIM; Zhang–Cai 1991 Thm A VERBATIM.*

**D7. k ≥ 4. — DEAD.**
The order of a weak focus in a quadratic system is at most three; M(2) = 3.
*Source: Kooij–Zegeling 1998 §1 VERBATIM ("It is well-known that the order of a weak focus in a quadratic
system is at most three"); Bautin, restated VERBATIM in Villanueva–Tucker 2026.*

**D8. Two hemicycles at once, giving 2 + 2 = 4 from graphics alone. — DEAD.**
Marín–Villadelprat Theorem C: the **simultaneous** cyclicity of the pair {Γu, Γℓ} is exactly 3, not 4
(each individually has cyclicity 2). *Source: Marín–Villadelprat 2025 Thm B + Thm C VERBATIM.*

### ALIVE (i.e. not excluded by anything I found — but nobody has realised them)

**A1. k = 1 (weak focus of order 1) + a graphic of cyclicity ≥ 3 surrounding the *other* focus.**
Not excluded by Zhang Pingguang (that theorem constrains only configurations where cycles actually
surround both foci), and the k=2/k=3 uniqueness theorems of Kooij–Zegeling explicitly do **not** cover
first-order weak foci with the two foci of the same stability (their Remark 1.2, VERBATIM). The largest
polycycle cyclicity ever computed in the quadratic family is **3** (Swirszcz's bicycle on the curve C;
Iliev's period-annulus cyclicity 3 at Q₄⁺ = (−4,2)). 1 + 3 = 4. To get 5 you would need a polycycle of
cyclicity 4 — **nobody has exhibited one, and nobody has excluded one.**
*Sources: Kooij–Zegeling 1998 Remark 1.2; Marín–Villadelprat 2025 §1 (Swirszcz, Iliev).*

**A2. k = 1 + finite saddle-node / infinite saddle-node families (QW1SN, QW1SN11, QW1SN02).**
Artés–Trullàs 2026 says VERBATIM "It can be shown that some of these parts have at least two limit
cycles", and the QW1SN11/QW1SN02 abstracts report 28 and 25 subsets with limit cycles and non-algebraic
bifurcation sets of separatrix connections *and* double limit cycles. Nobody reports > 2. This is the
frontier the Artés school is currently working; the counts are not published as bounds, so the recipe is
formally **unknown/alive** here. *Sources: Artés's own abstracts, §5.*

**A3. The QW2 "corners".**
The BIRS 2007 report records VERBATIM that Artés–Llibre–Schlomiuk *conjecture* that certain corners of the
QW2 bifurcation diagram "will produce the larger number of limit cycles in the family of quadratic vector
fields". This is the single most promising documented lead in slice B, and I could not read the paper.
Requires ALS 2006, IJBC 16(11), 3127–3194. *Source: BIRS 2007 report §3.6 VERBATIM.*

**A4. Degenerate corners of QW3 whose graphics have unproven cyclicity.**
LS2004 §10 VERBATIM: "the program outlined in [15] ... has stopped short (so far) of proving the finite
cyclicity of graphics present in these more degenerate systems, such as for example the system
corresponding to the point Σ = [0:1:0] in Figure 1." That point is on the boundary of QW3 and its graphic
has *no* proved cyclicity bound at all. Alive by ignorance. *Source: LS2004 §10 VERBATIM.*

**A5. The ~30 still-open DRR graphics.**
All 121 surround an anti-saddle. Roughly 33 were open in 2015 (88 proved); (I¹₁₂), (I¹₁₃), (I¹₁₄) closed in
2015; (H³₁₃), (I¹₆b), (DI₂b) only partial; (H³₁₄) claimed 2026 in an unrefereed preprint. Every open one
could in principle have cyclicity ≥ 4. But note: "finite cyclicity unproved" ≠ "cyclicity large"; the
obstruction in every documented case is the analysis (nilpotent blow-ups, four second-type Dulac maps,
2-D displacement maps), not evidence of many cycles.
*Sources: Rousseau–Shan–Zhu 2015 (count 88); Roussarie–Rousseau 2015 (named open cases); Lu 2026 (H³₁₄).*

### UNKNOWN (need a source I could not obtain)

**U1. Does any portrait in the closure of QW2 have a graphic surrounding the weak focus of order two?**
**UNRESOLVED IN THIS SWEEP.** Requires ALS 2006 (paywalled, no repository copy, confirmed by OpenAlex,
Semantic Scholar, ddd.uab.cat and the full 2008–2012 UAB preprint archive). The author's own supplementary
page hosts only figure archives, explicitly excluding the paper "due to the copyright".
A secondary, **unverified** claim circulating in search summaries is that a QW2 system has at most one
limit cycle around the weak focus and at most two in the whole plane; if that is in the paper, U1 is
answered "no" and QW2 caps at 2 + 1 = 3 (in QW2 proper) or 4 (in its closure, where QW3 lives).

**U2. Cyclicity of the QW^I graphics.** ALS 2010 paywalled. But note the class is self-limiting: an
invariant straight line **forces the weak focus to order 1** (VERBATIM from the authors' abstract), so the
best conceivable stack there is 1 + m.

**U3. Zhang Pingguang's "unbounded separatrix cycles in quadratic systems" (1993, Chinese).**
Directly on-topic (graphics through infinity) and completely unobtained. Listed only in Reyn's preprint
bibliography.

### The arithmetic, in one line

Every closed-form route documented in the literature tops out at **4**:
`3 (Bautin, weak focus of order 3) + 1 (elementary graphic of cyclicity 1 around the strong focus)`,
or `1 (weak focus of order 1) + 3 (Swirszcz bicycle / Iliev period annulus)`, or
`3 (simultaneous hemicycle pair)` + at most one more elsewhere.
**Reaching 5 requires either a polycycle of quadratic cyclicity ≥ 4 (none known, none excluded), or a
violation of Zhang Pingguang's (1, n) theorem, or 4 cycles around a single non-weak focus.**
