# Independent Galias–Tucker sign-convention check

Source inspected: all extracted text of the eight-page final 2022 article and page images 3–4, DOI 10.1016/j.amc.2021.126691. This is a narrow source-consistency check, not a reimplementation of the proof.

Section 2, p.2 defines P by the smallest **positive** return time. Section 3.1 and Figure 1 define f(y)=y−P(y). Section 3.3, p.4 applies interval Newton to that same f and reports the derivative/stability of P in Lemma 3. The inverse map appears separately in Section 4, Lemma 5, p.5; it is not used as the meaning of P in Section 3.

Lemma 2, p.3 labels each endpoint difference as P(y)−y, while all four pairs have the signs appropriate to y−P(y) when compared with Lemma 3 and Figure 1:

| Bracket (paper numbering) | Lemma 3 type | Signs expected for P(y)−y, left/right | Printed Lemma 2 signs |
|---|---|---|---|
| y1, near 0.04268960388208 | stable | +, − | −, + |
| y2, near 6.66666e−8 | unstable | −, + | +, − |
| y3, near 2.2478e−21 | stable | +, − | −, + |
| y4, near 7.0711e−75 | unstable | −, + | +, − |

This supports a probable sign-label transcription error in Lemma 2. It supplies no reason by itself to discard the brackets or existence statement: multiplying the displacement function by −1 preserves its zeros, opposite endpoint signs, and interval-Newton operator. Stability should be quoted from Lemma 3 or independently evaluated forward multipliers, not inferred from the mislabeled endpoint display. No alternative inverse-map convention in the text resolves the discrepancy. An author-confirmed erratum was not searched for in this narrow check.
