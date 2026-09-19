"""Independent exact checks for Sections 4--6 of manuscript.tex.

This deliberately does not import the project's discovery/verifier scripts.
It checks the multinomial Jacobsthal quotient as an integer-difference
valuation, exhaustive paired arrays for small row sums, and the full
coefficient theorem (including zero and out-of-support targets).
"""

from functools import lru_cache
from itertools import product
from math import factorial


def vp(z: int, p: int) -> int:
    if z == 0:
        return 10**9
    z = abs(z)
    e = 0
    while z % p == 0:
        z //= p
        e += 1
    return e


def comps(total: int, length: int = 4):
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in comps(total - first, length - 1):
            yield (first,) + tail


@lru_cache(maxsize=None)
def multinomial(parts) -> int:
    ans = factorial(sum(parts))
    for z in parts:
        ans //= factorial(z)
    return ans


def scaled(parts, scale: int):
    return tuple(scale * z for z in parts)


def check_multinomial_jacobsthal():
    tested = 0
    sharp = []
    for p, max_s in ((3, 3), (5, 2), (7, 2), (11, 1)):
        for s in range(max_s + 1):
            for total in range(1, 9):
                for v in comps(total):
                    lo = multinomial(scaled(v, p**s))
                    hi = multinomial(scaled(v, p ** (s + 1)))
                    assert vp(lo, p) == vp(hi, p), (p, s, v, lo, hi)
                    excess = vp(hi - lo, p) - vp(lo, p)
                    assert excess >= 2 * s + 2, (p, s, v, excess)
                    if excess == 2 * s + 2:
                        sharp.append((p, s, v))
                    tested += 1
    return tested, sharp[:8]


def sign_exponent(rows) -> int:
    return (
        rows[0][2]
        + rows[1][1]
        + rows[2][3]
        + rows[3][1]
        + rows[3][2]
        + rows[3][3]
    )


def check_paired_arrays():
    tested = 0
    for p in (3, 5, 7):
        # L is the lower-level row sum N/p.  Exhaust L <= 3.
        for lower_n in range(1, 4):
            rows = list(comps(lower_n))
            for u in product(rows, repeat=4):
                lower_a = sum(row[0] + row[1] for row in u)
                lower_b = sum(row[0] + row[2] for row in u)
                rminus1 = min(vp(lower_n, p), vp(lower_a, p), vp(lower_b, p))
                R = 1 + rminus1
                flat_nonzero = [z for row in u for z in row if z]
                s = min(vp(z, p) for z in flat_nonzero)
                assert s <= R - 1, (p, lower_n, lower_a, lower_b, s, R)

                lo = 1
                hi = 1
                for row in u:
                    lo *= multinomial(row)
                    hi *= multinomial(scaled(row, p))
                assert vp(hi - lo, p) >= 2 * R, (
                    p,
                    lower_n,
                    lower_a,
                    lower_b,
                    s,
                    R,
                    vp(hi - lo, p),
                )
                assert sign_exponent(u) % 2 == sign_exponent(
                    tuple(scaled(row, p) for row in u)
                ) % 2
                tested += 1
    return tested


def check_primitive_lemma():
    """Min-plus exhaustive check, independent of coefficient cancellation."""
    checked_targets = 0
    minima = []
    for p, N in ((3, 3), (3, 6), (3, 9), (5, 5), (7, 7)):
        row_data = []
        for row in comps(N):
            row_data.append(
                (
                    row[0] + row[1],
                    row[0] + row[2],
                    vp(multinomial(row), p),
                    any(z % p for z in row),
                )
            )
        # State: (x exponent, y exponent, globally primitive) -> min valuation.
        states = {(0, 0, False): 0}
        for _ in range(4):
            new = {}
            for (A, B, primitive), oldv in states.items():
                for da, db, rowv, row_primitive in row_data:
                    key = (A + da, B + db, primitive or row_primitive)
                    val = oldv + rowv
                    if val < new.get(key, 10**9):
                        new[key] = val
            states = new
        local_min = 10**9
        for (A, B, primitive), val in states.items():
            if not primitive or A % p or B % p:
                continue
            R = min(vp(N, p), vp(A, p), vp(B, p))
            assert val >= 2 * R, (p, N, A, B, R, val)
            local_min = min(local_min, val - 2 * R)
            checked_targets += 1
        minima.append((p, N, local_min))
    return checked_targets, minima


def multiply(f, g):
    out = {}
    for (a, b), ca in f.items():
        for (c, d), cb in g.items():
            key = (a + c, b + d)
            out[key] = out.get(key, 0) + ca * cb
    return {k: v for k, v in out.items() if v}


def power(f, n):
    ans = {(0, 0): 1}
    base = f
    while n:
        if n & 1:
            ans = multiply(ans, base)
        n //= 2
        if n:
            base = multiply(base, base)
    return ans


def build_P():
    factors = [
        {(1, 1): 1, (1, 0): 1, (0, 1): -1, (0, 0): 1},
        {(1, 1): 1, (1, 0): -1, (0, 1): 1, (0, 0): 1},
        {(1, 1): 1, (1, 0): 1, (0, 1): 1, (0, 0): -1},
        {(1, 1): 1, (1, 0): -1, (0, 1): -1, (0, 0): -1},
    ]
    ans = {(0, 0): 1}
    for factor in factors:
        ans = multiply(ans, factor)
    return ans


def check_full_coefficients():
    P = build_P()
    tested = 0
    # Include N=0, zero coordinates, boundary coordinates, and one point
    # beyond each side of the support square.
    for p, max_lower in ((3, 4), (5, 2), (7, 1)):
        powers = {n: power(P, n) for n in range(max_lower + 1)}
        powers.update({p * n: power(P, p * n) for n in range(max_lower + 1)})
        mod = p * p
        for n in range(max_lower + 1):
            lo = powers[n]
            hi = powers[p * n]
            max_ab = 4 * n
            targets = range(0, max_ab + 2)
            for a in targets:
                for b in targets:
                    lhs = hi.get((p * a, p * b), 0)
                    rhs = lo.get((a, b), 0)
                    assert (lhs - rhs) % mod == 0, (p, n, a, b, lhs, rhs)
                    tested += 1
    return tested


def main():
    jac, sharp = check_multinomial_jacobsthal()
    primitive, primitive_slack = check_primitive_lemma()
    paired = check_paired_arrays()
    coeff = check_full_coefficients()
    print(f"multinomial Jacobsthal cases: {jac}")
    print(f"sample cases sharp at 2s+2: {sharp}")
    print(f"primitive-sector target minima checked: {primitive}")
    print(f"minimum slack over 2R by (p,N): {primitive_slack}")
    print(f"exhaustive paired-array cases: {paired}")
    print(f"full coefficient cases incl. zero/outside support: {coeff}")
    print("all independent checks passed")


if __name__ == "__main__":
    main()
