# Q4 numerical controls

All numbers here are floating-point regression data, not interval proofs.
Regenerate them with:

    python q4/q4_controls.py --output q4/data/controls.json

## Positive control: three simple zeros

Zhao proves existence of three zeros through the center-endpoint hierarchy
\(0<|\nu_1|\ll|\nu_2|\ll|\nu_3|\ll|\nu_4|\), but publishes no finite
decimal/rational vector. We therefore constructed a finite control in the
same four-function space rather than mislabeling it as published data.

At \(\kappa=4\), setting \(\mu_4=1\) and solving
\(I(2)=I(3)=I(3.7)=0\) at 70 decimal digits gives

\[
\mu=\begin{pmatrix}
-10.879545578484009594356916559918594165243840226638\\
 17.285882261952706527844435811180895953807462874156\\
-18.538781607158024636097957631356781574918503088031\\
1
\end{pmatrix}.
\]

Independent grid/bracket isolation returned

\[
 s=2.000000000002993,\quad2.999999999979691,\quad
 3.700000000013739.
\]

The 70-digit construction residuals at the three forced values are all below
\(1.8\times10^{-69}\).

The maximum reproduced Q4 zero count is therefore **three**, matching the
published lower bound. This control is outside the \(\kappa<85/23\)
five-candidate strip; it validates integral evaluation, not the filter.

## Exact negative control

Take \(\mu=(1,0,0,0)\). Then \(I(h)=hI_{00}(h)<0\) throughout the open
annulus because \(h<0\) and \(I_{00}>0\). It has exactly zero interior zeros.

## Degenerate control

At \(\kappa=4,s_0=2.5\), a numerical solve of
\(I(s_0)=I'(s_0)=0\), with \(\mu_3=0,\mu_4=1\), gives

\[
 \mu=(-13.875035391538406,-3.186357934775721,0,1).
\]

The stored residual is zero to double precision and the centered curvature
estimate is \(I''(s_0)=0.02488010908408\ne0\). It exercises tangency and
multiplicity diagnostics; it is not a rigorous double-root certificate.

## Independent evaluator check

At \(\kappa=4,s=2.5\), the four basis values from 60-digit adaptive area
quadrature are

\[
(-0.1394948031904048,\ 0.25107246911387726,\
 0.2547155035116732,\ -1.1354885770378151).
\]

Fixed Gauss area quadrature differs by at most \(1.32\times10^{-14}\). A
separate DOP853 Hamiltonian-orbit integration, converting enclosed areas to
Green line integrals, differs by at most \(5.54\times10^{-13}\).

The machine-readable record is q4/data/controls.json.
