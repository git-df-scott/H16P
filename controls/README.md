# Numerical control

`reproduce_four_cycle.py` reproduces four Poincare-map fixed points in the
five-parameter quadratic family used by Kuznetsov, Kuznetsova, and Leonov.

Run:

```bash
python --version  # control generated with Python 3.12.13
python -m pip install -r controls/requirements.txt
python controls/reproduce_four_cycle.py
```

The script writes `four_cycle_control.json` and `four_cycle_control.png`.
Every output is explicitly **non-rigorous**. Floating-point integration and a
small return residual do not prove a periodic orbit. The control exists to
test candidate-generation code before candidates enter the interval gate in
[`../RIGOROUS_CERTIFICATION.md`](../RIGOROUS_CERTIFICATION.md).
