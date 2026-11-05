# The survival probability of a heavy-tailed critical Derrida-Retaux system

Yitzchak Shmalo

For the Derrida-Retaux recursion X_{n+1} = (X_n^(1) + ... + X_n^(m) - 1)^+ in the critical regime (m-1) E[X_0 m^{X_0}] = E[m^{X_0}] < infinity, Chen, Derrida, Hu, Lifshits and Shi (*A max-type recursive model: some properties and open questions*, arXiv:1705.04787, Conjecture 1) conjectured that

    P(X_n != 0) ~ 4 / ((m-1)^2 n^2).

The paper shows that this fails under the critical moment condition alone. For m = 2 and the nonconstant critical law

    P(X_0 = 0) = cB,  P(X_0 = 1) = 0,  P(X_0 = k) = c 2^{-k} k^{-3}  (k >= 2),

with B = zeta(2) - zeta(3) and c the normalizing constant, one has E[X_0^3 2^{X_0}] = infinity and

    liminf n^2 P(X_n != 0) <= 3/2.

The value 3/2 is the constant alpha(alpha-2)/2 at alpha = 3 that Chen, Hu and Shi (PTRF 182 (2022), (1.5)) predict for this tail. The proof combines an exact recursion for E[2^{X_n}] with the upper bound of Chen and Shi (*The stable Derrida-Retaux system at criticality*, Theorem 1.1) on the product of E[2^{X_i}], i < n. Conjecture 1 under the additional condition E[X_0^3 m^{X_0}] < infinity is not addressed.

## Contents

- `paper/main.tex`, `paper/main.pdf`: the paper.
- `code/checks.py`: the constants of the law and the criticality identity to 50 digits; in exact rational arithmetic, the identity behind the recursion and the one-step barrier inequality for two parameter pairs (16,800 cases); the truncated third moment identity. Output in `code/checks_out.txt` and `code/checks.json`.
- `lean/Barrier.lean`: the polynomial identity behind the recursion and the barrier inequality for (a, beta) = (3, 3/2), checked with Lean 4 and Mathlib. The analysis and the cited theorem are not formalized.

## Running

    cd code
    python checks.py > checks_out.txt

The program needs `mpmath` and runs in a few seconds. The Lean file compiles inside any Lean 4 project with Mathlib (`lake env lean Barrier.lean`).

## License

Paper: CC BY 4.0. Code: MIT.
