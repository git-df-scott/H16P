# Strike 7 replay

See [the proof and its scope](../../Q4_SEVENTH_LIMITING_FACE.md).
The finite-lift route remains open.

From the repository root:

```sh
python q4/seventh/a1_determinants.py
python q4/seventh/finite_basis.py
python q4/seventh/parameter_connection.py
```

These use SymPy 1.14.0 and exact rational arithmetic. The first checks four
ODEs, eight center data, two determinant factorizations, the Jensen
concavity identity, and all 29 positive polynomial coefficients. The
second checks the finite-lift closed subspace. The third derives a rational
parameter connection and checks its differential compatibility and center
data; it does not establish a sign comparison.

`explore_upper.py` is a bounded numerical diagnostic of a failed upper-bound
argument. Its 18 frozen points include positive upper bounds paired with
negative actual determinants. It imports the Strike-6 evaluator and uses
the versions recorded in `../sixth/environment.json`. These samples neither
prove nor disprove a global determinant inequality.

The JSON files are direct replay outputs. No numerical samples enter the
limiting-face or closed-subspace proofs. No independent audit is claimed.
