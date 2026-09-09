# Negative-perturbation fold search and exact center identification

No counterexample or outer double-zero candidate was found. Six negative-perturbation searches completed with large normalized residuals. They do not exclude a fold elsewhere.

The same three interior return equalities were retained while epsilon varied over [-0.8,-0.01], with moduli r=0.5,1,2 and two starting radii per modulus. The initial epsilon was -0.2. The proposed outer gap divided by the outer anchor was bounded between 0.1 and 9. Complete half-return gates and the event-corrected radial derivative from the previous experiment were used.

To avoid accepting the known reversible-center locus, the search rejected control triples within 1e-4 of (tau,w,r+u)=(0,0,0). Its objective also divided out the distance to this locus relative to the reference field at epsilon=-0.1, in addition to the previous amplitude and reference-displacement normalization. This handles one known center locus; it is not a complete quadratic-center classification.

All six searches moved toward epsilon=-0.01 and the largest allowed proposed radius. The best maximum absolute normalized fold residuals were approximately 0.8285, 0.8847, and 0.9485 for the three moduli. None approached the 1e-6 candidate threshold. The smallest evaluated distance to the known center locus was 0.3073, so the failed construction was not merely convergence onto that locus. There were 5,472 parameter-sensitivity half-passages and 636 radial-variation half-passages. No five-cycle claim follows from optimizer termination.

## Exact algebraic result

On tau=w=0 and u=-r, the field is

P=-y-(2+r²)x²+y², Q=x(1+B y), where B=-1-3r²-epsilon.

For B nonzero, the affine coordinates X=-Bx/2 and Y=(1+By)/2, with unchanged time, give

X'=(b-2)/4+(1-b)Y+aX²+bY², Y'=-2XY,

a=-2(r²+2)/(epsilon+3r²+1), b=2/(epsilon+3r²+1).

Both transformed equations were verified as exact symbolic identities. Thus this limiting family is contained in the reversible family already investigated; it is not a new family of isolated cycles. The original equilibria (0,0) and (0,1) have zero trace and determinants 1 and epsilon+3r². Reversibility and ellipticity give local centers at both when epsilon>-3r². This statement gives no global cycle bound for nearby fields.

The present finite-Q4 fold attempts supply no promising double cycle. A distinct mechanism worth examining is compatibility between finite-annulus return zeros and the simultaneous boundary bifurcations of the reversible family. Their separate possible counts cannot simply be added: the same perturbation coefficients must realize them together.

Scripts and JSON records preserve all fits, normalized residuals, center distances, and the exact coordinate identity. The pinned source repository was not modified.
