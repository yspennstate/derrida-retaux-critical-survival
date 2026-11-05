"""Checks for "The survival probability of a heavy-tailed critical Derrida-Retaux system".

1. The constants of the law P(X0=0)=cB, P(X0=1)=0, P(X0=k)=c 2^-k k^-3 (k>=2) to 50 digits,
   and the criticality identity E[2^X0] = E[X0 2^X0].
2. In exact rational arithmetic: the identity
       (1+e)^2/2 + (1-q)^2/2 - 1 = e + e^2/2 - q + q^2/2
   behind the recursion for e_n = E[2^X_n] - 1, and the barrier step
       d_{n+1} <= d_n - gamma/n^2,  d_n = e_n - beta/n,  gamma = a - beta - beta^2/2,
   whenever 0 <= e_n <= beta/n and q_n - q_n^2/2 >= a/n^2, for (a, beta) = (3, 3/2) and (7/4, 11/10).
3. The truncated moment E[X0^3 2^X0 1{X0<=k}] = c(k-1) for 2 <= k <= 600 (exact, in units of c).

Writes checks.json next to this file. Needs mpmath for part 1.
"""
import json
import pathlib
from fractions import Fraction as F

ROOT = pathlib.Path(__file__).resolve().parent
out = {}

# 1. constants
import mpmath as mp

mp.mp.dps = 50
A = mp.polylog(3, mp.mpf("0.5")) - mp.mpf("0.5")
B = mp.zeta(2) - mp.zeta(3)
c = 1 / (A + B)
m0 = c * B + c * (mp.zeta(3) - 1)
m1 = c * (mp.zeta(2) - 1)
# direct partial sums as an independent check of the closed forms
A_sum = mp.nsum(lambda k: mp.mpf(2) ** (-k) / k ** 3, [2, mp.inf])
B_sum = mp.nsum(lambda k: (k - 1) / mp.mpf(k) ** 3, [2, mp.inf])
assert abs(A - A_sum) < mp.mpf(10) ** -40 and abs(B - B_sum) < mp.mpf(10) ** -40
assert abs(c * (A + B) - 1) < mp.mpf(10) ** -45
assert abs(m0 - m1) < mp.mpf(10) ** -45
out["constants"] = {
    "A": mp.nstr(A, 30), "B": mp.nstr(B, 30), "c": mp.nstr(c, 30),
    "P(X0=0)": mp.nstr(c * B, 30), "P(X0=2)": mp.nstr(c / 32, 30),
    "E[2^X0]": mp.nstr(m0, 30), "E[X0 2^X0]": mp.nstr(m1, 30),
}
print("constants", out["constants"], flush=True)

# 2. exact identity and barrier step
cases = 0
for a, beta in [(F(3), F(3, 2)), (F(7, 4), F(11, 10))]:
    gamma = a - beta - beta * beta / 2
    assert gamma > 0
    for n in range(3, 203):
        r, rn = F(1, n), F(1, n + 1)
        b, bn = beta * r, beta * rn
        for q in (a * r * r + a * a * r ** 4, 2 * a * r * r):
            assert 0 <= q <= 1 and q - q * q / 2 >= a * r * r
            for j in range(21):
                e = b * F(j, 20)
                e_next = e + e * e / 2 - q + q * q / 2
                assert (1 + e) ** 2 / 2 + (1 - q) ** 2 / 2 - 1 == e_next
                assert e_next - bn <= (e - b) - gamma * r * r
                cases += 1
out["barrier_cases"] = cases
print("exact identity and barrier step:", cases, "cases passed", flush=True)

# 3. truncated third moment, in units of c: sum_{j=2}^{k} j^3 2^j * 2^-j j^-3 = k-1
for k in range(2, 601):
    assert sum(F(j ** 3 * 2 ** j, 2 ** j * j ** 3) for j in range(2, k + 1)) == k - 1
out["truncated_moment_k_max"] = 600
print("truncated moment identity: k = 2..600 passed", flush=True)

(ROOT / "checks.json").write_text(json.dumps(out, indent=2), encoding="utf8")
print("all checks passed", flush=True)
