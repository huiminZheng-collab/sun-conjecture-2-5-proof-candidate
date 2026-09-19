"""Independent exact checks for the balanced-residue lemma.

This script was written for the adversarial referee pass.  It does not use
the project's coefficient-congruence verifier.  All computations use exact
integers.
"""

from itertools import product
from math import factorial


def vp(n, p):
    if n == 0:
        return 10**9
    out = 0
    while n % p == 0:
        out += 1
        n //= p
    return out


def compositions4(n):
    for a in range(n + 1):
        for b in range(n - a + 1):
            for c in range(n - a - b + 1):
                yield a, b, c, n - a - b - c


def multinomial(row):
    n = sum(row)
    out = factorial(n)
    for e in row:
        out //= factorial(e)
    return out


def residue_state_minimum(modulus):
    """Minimum positive T for four residue rows satisfying both balances."""
    row_states = set()
    for row in product(range(modulus), repeat=4):
        mass = sum(row)
        if mass % modulus == 0:
            a, b, c, _ = row
            row_states.add(((a + b) % modulus, (a + c) % modulus,
                            mass // modulus, mass > 0))

    states = {(0, 0, False): 0}
    for _ in range(4):
        nxt = {}
        for (x, y, nonzero), old_t in states.items():
            for dx, dy, dt, row_nonzero in row_states:
                key = ((x + dx) % modulus, (y + dy) % modulus,
                       nonzero or row_nonzero)
                value = old_t + dt
                if value < nxt.get(key, 10**9):
                    nxt[key] = value
        states = nxt
    return states[(0, 0, True)]


def exact_level_identity(row, p):
    """Check Legendre valuation equals the sum of residue-level carries."""
    n = sum(row)
    actual = vp(multinomial(row), p)
    carries = []
    modulus = p
    while modulus <= n:
        residue_mass = sum(e % modulus for e in row) - (n % modulus)
        assert residue_mass % modulus == 0
        carries.append(residue_mass // modulus)
        modulus *= p
    assert actual == sum(carries), (row, p, actual, carries)
    return actual, carries


def all_target_minima(n, p, r):
    """Minimise valuations for every p^r-balanced primitive target."""
    scale = p**r
    rows = []
    for row in compositions4(n):
        a, b, c, _ = row
        rows.append((a + b, a + c, vp(multinomial(row), p),
                     any(e % p for e in row)))

    states = {(0, 0, False): 0}
    for _ in range(4):
        nxt = {}
        for (x, y, primitive), old_v in states.items():
            for dx, dy, dv, row_primitive in rows:
                key = x + dx, y + dy, primitive or row_primitive
                value = old_v + dv
                if value < nxt.get(key, 10**9):
                    nxt[key] = value
        states = nxt

    checked = 0
    minimum = 10**9
    witnesses = []
    for (a, b, primitive), value in states.items():
        if primitive and a % scale == 0 and b % scale == 0:
            assert value >= 2 * r, (n, p, r, a, b, value)
            checked += 1
            if value < minimum:
                minimum = value
                witnesses = [(a, b)]
            elif value == minimum:
                witnesses.append((a, b))
    return checked, minimum, witnesses


def sharpness_check(p, r):
    n = p**r
    rows = (
        (0, 0, 0, n),
        (0, 0, 1, n - 1),
        (n - 1, 1, 0, 0),
        (n, 0, 0, 0),
    )
    assert all(sum(row) == n for row in rows)
    a = sum(row[0] + row[1] for row in rows)
    b = sum(row[0] + row[2] for row in rows)
    value = sum(vp(multinomial(row), p) for row in rows)
    assert (a, b, value) == (2 * n, 2 * n, 2 * r)


def multiply_polynomials(left, right):
    out = {}
    for (a, b), x in left.items():
        for (c, d), y in right.items():
            key = a + c, b + d
            out[key] = out.get(key, 0) + x * y
    return {key: value for key, value in out.items() if value}


def exact_p_polynomial():
    factors = (
        {(1, 1): 1, (1, 0): 1, (0, 1): -1, (0, 0): 1},
        {(1, 1): 1, (1, 0): -1, (0, 1): 1, (0, 0): 1},
        {(1, 1): 1, (1, 0): 1, (0, 1): 1, (0, 0): -1},
        {(1, 1): 1, (1, 0): -1, (0, 1): -1, (0, 0): -1},
    )
    out = {(0, 0): 1}
    for factor in factors:
        out = multiply_polynomials(out, factor)
    return out


def support_wording_counterexample():
    """Exact-support holes need not persist under odd-prime dilation."""
    base = exact_p_polynomial()
    power = {(0, 0): 1}
    fourth = twelfth = None
    for exponent in range(1, 13):
        power = multiply_polynomials(power, base)
        if exponent == 4:
            fourth = power.get((2, 2), 0)
        if exponent == 12:
            twelfth = power.get((6, 6), 0)
    assert fourth == 0
    assert twelfth == -291456
    assert twelfth % 9 == 0
    return fourth, twelfth


def main():
    # The residue proof uses no primality: test every modulus 2,...,18.
    residue_results = {}
    for modulus in range(2, 19):
        residue_results[modulus] = residue_state_minimum(modulus)
        assert residue_results[modulus] == 2
    print("residue minima:", residue_results)

    # Exhaust every primitive target with both balances divisible by p^r.
    regimes = ((3, 3, 1), (5, 5, 1), (7, 7, 1),
               (6, 3, 1), (10, 5, 1), (9, 3, 2))
    for n, p, r in regimes:
        checked, minimum, targets = all_target_minima(n, p, r)
        print((n, p, r), "targets", checked, "minimum", minimum,
              "first sharp targets", targets[:8])

    # Exact Legendre-level identity, including levels above R.
    level_checks = 0
    for p in (2, 3, 5, 7):
        for n in range(0, 19):
            for row in compositions4(n):
                exact_level_identity(row, p)
                level_checks += 1
    print("Legendre level identities:", level_checks)

    for p in (2, 3, 5, 7):
        for r in range(0, 5):
            sharpness_check(p, r)
    print("sharpness witnesses: 20 passed")

    print("exact-support wording counterexample:",
          support_wording_counterexample())


if __name__ == "__main__":
    main()
