# Lane 1 VALIDATION
engine: lane1/retmap1.c sha256: d03c1c850573872ebbf298e1a7392fba89fca1c7913c9c008271664af34fec2c
defaults: {"rtol": 1e-12, "Tmax": 400.0, "Rmax": 10000.0, "nstep": 4000000}

## Cherkas-Artes-Llibre 2003, rows 1-8

| row | phi | s_max | #brackets | cycle x (this engine) | published x | max |dx| | ref-engine agrees |
|---|---|---|---|---|---|---|---|
| 1 | +x | >=1000 (unbounded on the grid) | 3 | 1.2809, 2.0070, 4.0193 | 1.26, 1.98, 3.95 | 6.93e-02 | True |
| 2 | +x | 40.2120 | 3 | 1.1935, 2.0596, 3.0896 | 1.4, 1.9, 3.1 | 2.06e-01 | True |
| 3 | -x | 0.7075 | 3 | 0.3235, 0.6155, 0.8482 | 0.32, 0.66, 0.8 | 4.82e-02 | True |
| 4 | -x | 0.5337 | 3 | 0.5569, 0.7466, 0.8523 | 0.56, 0.75, 0.87 | 1.77e-02 | True |
| 5 | -x | 0.4324 | 3 | 0.6278, 0.8018, 0.8729 | 0.63, 0.8, 0.88 | 7.10e-03 | True |
| 6 | +x | >=1000 (unbounded on the grid) | 3 | 1.0504, 1.1817, 1.5415 | 1.05, 1.16, 1.5 | 4.15e-02 | True |
| 7 | +x | >=1000 (unbounded on the grid) | 3 | 1.5108, 2.2316, 4.4763 | 1.28, 2.15, 4.43 | 2.31e-01 | True |
| 8 | +x | >=1000 (unbounded on the grid) | 3 | 1.3573, 2.3071, 4.1455 | 1.29, 2.22, 4.63 | 4.84e-01 | True |

## How much of the deviation is the paper's own coefficient rounding?

The published a, a20, a11, a01, a10 are printed to a fixed number of decimals.
For each row every printed coefficient is moved by half of its last printed
digit, one at a time, and the induced motion of the three cycle abscissae is
recorded.  This is the smallest honest error bar on a comparison with the table.

| row | this engine | published | deviation | rounding envelope | within envelope+0.02 |
|---|---|---|---|---|---|
| 1 | 1.2809, 2.0070, 4.0193 | 1.26, 1.98, 3.95 | 0.0209, 0.0270, 0.0693 | 0.1320, 0.3260, 1.2259 | True |
| 2 | 1.1935, 2.0596, 3.0896 | 1.4, 1.9, 3.1 | 0.2065, 0.1596, 0.0104 | 0.1072, 0.5773, 0.8454 | False |
| 3 | 0.3235, 0.6155, 0.8482 | 0.32, 0.66, 0.8 | 0.0035, 0.0445, 0.0482 | 0.0658, 0.1354, 0.0708 | True |
| 4 | 0.5569, 0.7466, 0.8523 | 0.56, 0.75, 0.87 | 0.0031, 0.0034, 0.0177 | 0.0027, 0.0076, 0.0036 | True |
| 5 | 0.6278, 0.8018, 0.8729 | 0.63, 0.8, 0.88 | 0.0022, 0.0018, 0.0071 | 0.0022, 0.0083, 0.0076 | True |
| 6 | 1.0504, 1.1817, 1.5415 | 1.05, 1.16, 1.5 | 0.0004, 0.0217, 0.0415 | 0.0034, 0.0100, 0.0108 | False |
| 7 | 1.5108, 2.2316, 4.4763 | 1.28, 2.15, 4.43 | 0.2308, 0.0816, 0.0463 | 0.2238, 0.5153, 0.6580 | True |
| 8 | 1.3573, 2.3071, 4.1455 | 1.29, 2.22, 4.63 | 0.0673, 0.0871, 0.4845 | 0.0737, 0.4834, 0.6305 | True |

## KKL control (Kuznetsov et al.)
- origin nest, ray +x: s_max=1e+05, brackets=3, r=['0.6832', '2.1837', '15.9628'], ref agrees=True
- origin nest, ray -x: s_max=1e+05, brackets=3, r=['0.3377', '0.5199', '0.7074'], ref agrees=True
- published origin cycles r = [0.6832, 2.1837, 15.9628] (ray +x) -- exact match to 4 d.p.
- second focus B = (-6.259641, 7.449768)
- remote nest, ray -x from B: brackets=1, s=['3706'], crossing x = ['-3712'], ref agrees=True
- repo's recorded remote section coordinate -3711.56 (same cycle, different section: theirs is the x-axis from the origin, mine the -x ray from B)

## Andronov-Hopf curves beta*(s): interior extrema

| seed | phi | s range | resolved | interior extrema | prominences | height range |
|---|---|---|---|---|---|---|
| cherkas1 | +x | [1e-3, 1000] | 139/240 | 2 | 3.44e-01, 1.24e-01 | 0.0003026 |
| cherkas2 | +x | [1e-3, 40.21] | 240/240 | 2 | 1.72e-02, 1.72e-02 | 0.0003562 |
| cherkas3 | -x | [1e-3, 0.7075] | 239/240 | 2 | 3.81e-01, 4.17e-01 | 0.0003461 |
| cherkas4 | -x | [1e-3, 0.5337] | 240/240 | 2 | 6.33e-02, 9.60e-02 | 0.001227 |
| cherkas5 | -x | [1e-3, 0.4324] | 240/240 | 2 | 2.72e-01, 2.72e-01 | 0.0008709 |
| cherkas6 | +x | [1e-3, 1000] | 119/240 | 2 | 5.60e-02, 1.73e-01 | 0.0001012 |
| cherkas7 | +x | [1e-3, 1000] | 232/240 | 2 | 2.02e-03, 2.02e-03 | 0.002083 |
| cherkas8 | +x | [1e-3, 1000] | 240/240 | 2 | 3.39e-03, 3.39e-03 | 3.647e-05 |

## Row 4 vs the published Andronov-Hopf polynomial
- published degree-6 fit on [0.6,0.9] has 2 interior extrema at x = [0.6238, 0.8056]
- this engine's beta*(x) on [0.6,0.9] has 2 interior extrema at x = [0.6207, 0.7981]

## Two-engine agreement at the bracket endpoints (PROTOCOL rule 2)
- 24 bracket endpoints re-integrated with SciPy DOP853 (dense event location);
  worst relative difference in D between the two engines: 2.58e-05
- all Cherkas brackets reproduced: True

wall time 76.5s
