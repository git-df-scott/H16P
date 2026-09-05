# Strike 6 replay and proof audit

Base commit: `3b94f34`. Scope: route 4, the outside-lobe Q4 four-zero
problem. All new research files are in this directory and
[Q4_SIXTH_BOUNDARY_REDUCTION.md](../../Q4_SIXTH_BOUNDARY_REDUCTION.md).

From the repository root:

```bash
python q4/sixth/check_exact.py
python q4/sixth/explore_determinant.py
python q4/sixth/boundary_diagnostic.py
```

The first command checks exact identities and creates two rigorous
one-point coefficient certificates using rational intervals and analytic
series tails. The global deductions additionally require the proofs in
the two notes. The other two commands reproduce fixed numerical controls;
they are not a proof of the remaining global signs. Their parameter
lists are fixed in the scripts. No command runs an unbounded search.

The recorded environment uses Python 3.12, SymPy 1.14, mpmath 1.3,
NumPy 2.3.5, and SciPy 1.17. Exact versions are in `environment.json`.
Rerunning changes recorded runtimes but should preserve certified signs.

## Proof dependencies and checks against overclaiming

- The coefficient ratio increases with both baseline anchors and along
  the dangerous mixture. These are distinct monotonicity statements;
  both are used explicitly.
- The corner with both anchors at zero is taken after dividing `H` by
  `t^2`; the finite confluent matrix is nonsingular.
- The supersolution residual is nonnegative at `q0=185/108` and
  `q0=167/90`. It is **not** asserted nonnegative for `q0=19/10`.
- The comparison excludes a positive first maximum even at the stated
  threshold equality, since the solution comparison is strict.
- The second-anchor determinant normalization has a positive denominator.
  Its boundary at `s=r` is a derivative row; simply repeating the row
  would create a spurious zero determinant.
- The no-two-zeros argument uses the *strict tail margin* of Theorem N
  in confluent limits. Mere continuity of a negative quantity would not
  be sufficient.
- A vanishing cofactor vector is handled separately: the determinant is
  identically zero, which never passes the required strict positivity.
- The boundary reduction includes zero boundary values. It does not
  assert that the determinant is monotone in the second anchor.
- A positive determinant still needs positive baseline momentum and the
  later Green and original-integral height gates.
- No strict-fibre or boundary sample is promoted to a universal sign.
  Neither of the two residual sign problems is solved.

## Independent numerical control

The boundary ODE integrates the determinant directly in logarithmic time.
It records a separate final-product reconstruction to expose cancellation.
The `a=1,r=.95,s=1` control instead uses 45-digit mpmath quadrature of the
explicit limiting `P` and `Phi` kernels. It agrees with the ODE to about
`7.2e-18` absolute error. This is an independent method check at one point,
not a global certificate or a finite-lift counterexample.

The derivations were checked locally against the inherited definitions,
symbolically replayed where applicable, and numerically cross-checked as
described. No new sub-agent or outside referee audit was performed.
