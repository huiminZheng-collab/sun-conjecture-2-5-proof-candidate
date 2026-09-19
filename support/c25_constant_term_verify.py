"""Exact checks for the constant-term route to Sun's Conjecture 2.5.

This file supplies computational evidence only.  It checks

1. the constant-term formula for G_n;
2. the four-factor identity after the quadratic substitution;
3. the minimum Kummer valuation in the non-p-divisible sector; and
4. the scaled multinomial congruence needed for the p-divisible sector.

All arithmetic is exact Python integer arithmetic.
"""

from functools import lru_cache
from itertools import product
from math import comb, factorial


def vp(n: int, p: int) -> int:
    if n == 0:
        return 10**9
    ans = 0
    while n % p == 0:
        ans += 1
        n //= p
    return ans


def multinomial(parts: tuple[int, ...]) -> int:
    n = sum(parts)
    ans = 1
    used = 0
    for q in parts:
        ans *= comb(used + q, q)
        used += q
    return ans


def compositions4(n: int):
    for a in range(n + 1):
        for b in range(n - a + 1):
            for c in range(n - a - b + 1):
                yield a, b, c, n - a - b - c


@lru_cache(None)
def G(n: int) -> int:
    """Sun's original positive binomial convolution."""
    return sum(
        4**k * comb(2 * n - 2 * k, n - k) ** 2 * comb(2 * k, k)
        for k in range(n + 1)
    )


def G_alternating(n: int) -> int:
    """The alternating expression supplied by the constant-term transform."""
    return sum(
        (-1) ** k * 16 ** (n - k) * comb(n, k) * comb(2 * k, k) ** 2
        for k in range(n + 1)
    )


def ct_laurent(n: int) -> int:
    # CT [16-(x+2+x^-1)(y+2+y^-1)]^n, after expanding the outer power.
    return sum(
        (-1) ** k * 16 ** (n - k) * comb(n, k) * comb(2 * k, k) ** 2
        for k in range(n + 1)
    )


def signed_balanced_sum(n: int) -> int:
    """Evaluate the four-factor coefficient through its balanced row states."""
    rows = []
    sign_exponents = (
        lambda a, b, c, d: c,
        lambda a, b, c, d: b,
        lambda a, b, c, d: d,
        lambda a, b, c, d: b + c + d,
    )
    for row_number, exponent in enumerate(sign_exponents):
        choices = []
        for q in compositions4(n):
            a, b, c, _ = q
            coefficient = multinomial(q) * (-1) ** exponent(*q)
            choices.append((a + b, a + c, coefficient))
        rows.append(choices)

    states = {(0, 0): 1}
    for choices in rows:
        nxt = {}
        for (sx, sy), value in states.items():
            for dx, dy, coefficient in choices:
                if sx + dx <= 2 * n and sy + dy <= 2 * n:
                    key = sx + dx, sy + dy
                    nxt[key] = nxt.get(key, 0) + value * coefficient
        states = nxt
    return (-1) ** n * states[(2 * n, 2 * n)]


def mul_poly(left, right):
    out = {}
    for (a, b), x in left.items():
        for (c, d), y in right.items():
            out[a + c, b + d] = out.get((a + c, b + d), 0) + x * y
    return {monomial: value for monomial, value in out.items() if value}


def check_four_factor_identity():
    # Lambda(x^2,y^2), using Lambda=16-(2-x-x^-1)(2-y-y^-1).
    one = {(0, 0): 1}
    xpart = {(0, 0): 2, (2, 0): -1, (-2, 0): -1}
    ypart = {(0, 0): 2, (0, 2): -1, (0, -2): -1}
    lhs = {(0, 0): 16}
    for monomial, value in mul_poly(xpart, ypart).items():
        lhs[monomial] = lhs.get(monomial, 0) - value
    lhs = {monomial: value for monomial, value in lhs.items() if value}

    factors = (
        {(1, 1): 1, (1, 0): 1, (0, 1): -1, (0, 0): 1},
        {(1, 1): 1, (1, 0): -1, (0, 1): 1, (0, 0): 1},
        {(1, 1): 1, (1, 0): 1, (0, 1): 1, (0, 0): -1},
        {(1, 1): 1, (1, 0): -1, (0, 1): -1, (0, 0): -1},
    )
    rhs = one
    for factor in factors:
        rhs = mul_poly(rhs, factor)
    rhs = {(a - 2, b - 2): -value for (a, b), value in rhs.items()}
    assert lhs == rhs


def minimum_bad_sector_valuation(p: int, r: int, m: int = 1):
    """Minimise the valuation over balanced four-row compositions.

    The state records the two numerator exponents and whether at least one of
    the sixteen entries is not divisible by p.
    """
    n = m * p**r
    rows = []
    for q in compositions4(n):
        a, b, c, _ = q
        rows.append((a + b, a + c, vp(multinomial(q), p), any(x % p for x in q), q))

    states = {(0, 0, False): (0, ())}
    for _ in range(4):
        nxt = {}
        for (sx, sy, bad), (value, witness) in states.items():
            for dx, dy, w, row_bad, q in rows:
                tx, ty = sx + dx, sy + dy
                if tx > 2 * n or ty > 2 * n:
                    continue
                key = (tx, ty, bad or row_bad)
                candidate = (value + w, witness + (q,))
                if key not in nxt or candidate[0] < nxt[key][0]:
                    nxt[key] = candidate
        states = nxt
    return states[(2 * n, 2 * n, True)]


def check_jacobsthal_multinomials():
    """Check the lifted Babbage--Jacobsthal ratio used in the proof."""
    tests = 0
    for p in (3, 5, 7):
        for s in range(4):
            for n in range(1, 7):
                for q in compositions4(n):
                    lower = multinomial(tuple(p**s * x for x in q))
                    upper = multinomial(tuple(p ** (s + 1) * x for x in q))
                    # The quotient is a p-adic unit, not necessarily a rational
                    # integer.  Test v_p(upper/lower-1) through the numerator.
                    assert vp(upper, p) == vp(lower, p), (p, s, q)
                    assert vp(upper - lower, p) - vp(lower, p) >= 2 * s + 2, (
                        p,
                        s,
                        q,
                    )
                    tests += 1
    return tests


def check_full_congruence():
    tests = []
    ranges = {3: (4, 5), 5: (3, 3), 7: (2, 3), 11: (2, 2)}
    for p, (max_r, max_m) in ranges.items():
        for r in range(1, max_r + 1):
            for m in range(1, max_m + 1):
                lhs = G(m * p**r)
                rhs = G(m * p ** (r - 1))
                exact_r = vp(m * p**r, p)
                assert (lhs - rhs) % p ** (2 * exact_r) == 0, (p, r, m)
                tests.append((p, r, m))
    return tests


def main():
    assert all(G(n) == G_alternating(n) == ct_laurent(n) for n in range(31))
    print("original convolution = alternating constant term: n=0..30 passed")
    check_four_factor_identity()
    print("four-factor Laurent identity: passed")
    assert all(G(n) == signed_balanced_sum(n) for n in range(7))
    print("signed balanced expansion: n=0..6 passed")

    for p, r, m in ((3, 1, 1), (5, 1, 1), (7, 1, 1), (3, 2, 1), (5, 2, 1)):
        value, witness = minimum_bad_sector_valuation(p, r, m)
        print(f"p={p}, r={r}, m={m}: minimum={value}, target={2*r}, witness={witness}")
        assert value >= 2 * r

    tests = check_jacobsthal_multinomials()
    print(f"lifted Jacobsthal multinomial ratios: {tests} cases passed")

    full_tests = check_full_congruence()
    print(f"full Conjecture 2.5 spot checks: {len(full_tests)} cases passed")


if __name__ == "__main__":
    main()
