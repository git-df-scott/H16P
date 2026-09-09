# Correction: exact higher-order focus routes

**The proposed exact third-order-focus-plus-existing-outer-cycle route is ruled out by a published theorem.** I should have checked this before extending the numerical branch search.

Chengzhi Li's 1986 paper proves that no limit cycle surrounds a third-order weak focus of a real quadratic system. I retrieved the original 17-page scan and visually checked its opening statement; I have not independently audited the full proof. [Original journal record](https://camath.fudan.edu.cn/cambcn/ch/reader/view_abstract.aspx?file_no=7B205&flag=1).

Llibre and Schlomiuk's 2004 classification explicitly invokes that result on printed page 328. Its Theorem 18 carefully distinguishes proved neighborhood bounds from numerical evidence for other phase portraits; it must not be cited as an unconditional four-cycle bound for every nearby field. [Classification paper](https://doi.org/10.4153/CJM-2004-015-2).

Zhang and Zhao's 2001 publisher abstract also states that, with a second-order weak focus and a strong focus, the system has at most two cycles, distributed at most one around each. I verified the abstract, not the subscription-only proof. This obstructs the earlier exact second-order-focus-plus-two-existing-outer-cycles precursor when its strong-focus hypothesis applies. [Publisher abstract](https://link.springer.com/article/10.1007/s11766-001-0018-y).

The numerical calculations remain reproducible, but their proposed search directions are superseded. Near-zero coefficients in rational approximations do not establish exact higher-order foci; the theorems also must not be extended automatically to arbitrary finite perturbations.

**Research decision:** stop searching these exact-focus surfaces for the prohibited existing outer cycles. Keep the original Q4 closed-annulus problem separate. Any further weak-focus construction must identify a mechanism not excluded by the applicable theorem before computation begins. No counterexample has been found.
