"""Exact checks for the multivariate strengthening behind Conjecture 2.5.

For P=(xy+x-y+1)(xy-x+y+1)(xy+x+y-1)(xy-x-y-1), check

    [x^(a p^r) y^(b p^r)] P^(n p^r)
      == [x^(a p^(r-1)) y^(b p^(r-1))] P^(n p^(r-1))
         (mod p^(2r))

through every coefficient in several small regimes.
"""


def multiply(left, right, modulus=None):
    out = {}
    for (a, b), x in left.items():
        for (c, d), y in right.items():
            key = a + c, b + d
            out[key] = out.get(key, 0) + x * y
            if modulus is not None:
                out[key] %= modulus
    return {key: value for key, value in out.items() if value}


def polynomial_power(base, exponent, modulus):
    answer = {(0, 0): 1}
    power = {key: value % modulus for key, value in base.items()}
    while exponent:
        if exponent & 1:
            answer = multiply(answer, power, modulus)
        exponent //= 2
        if exponent:
            power = multiply(power, power, modulus)
    return answer


def four_factor_polynomial():
    factors = (
        {(1, 1): 1, (1, 0): 1, (0, 1): -1, (0, 0): 1},
        {(1, 1): 1, (1, 0): -1, (0, 1): 1, (0, 0): 1},
        {(1, 1): 1, (1, 0): 1, (0, 1): 1, (0, 0): -1},
        {(1, 1): 1, (1, 0): -1, (0, 1): -1, (0, 0): -1},
    )
    answer = {(0, 0): 1}
    for factor in factors:
        answer = multiply(answer, factor)
    return answer


def check_regime(p, r, n):
    modulus = p ** (2 * r)
    scale = p**r
    lower_scale = p ** (r - 1)
    base = four_factor_polynomial()
    upper = polynomial_power(base, n * scale, modulus)
    lower = polynomial_power(base, n * lower_scale, modulus)

    checks = 0
    for a in range(4 * n + 1):
        for b in range(4 * n + 1):
            lhs = upper.get((a * scale, b * scale), 0)
            rhs = lower.get((a * lower_scale, b * lower_scale), 0)
            assert (lhs - rhs) % modulus == 0, (p, r, n, a, b, lhs, rhs)
            checks += 1
    return checks


def main():
    regimes = (
        (3, 1, 1),
        (3, 1, 2),
        (3, 2, 1),
        (3, 2, 2),
        (5, 1, 1),
        (5, 1, 2),
        (5, 2, 1),
        (7, 1, 1),
    )
    total = 0
    for regime in regimes:
        checks = check_regime(*regime)
        total += checks
        print(f"{regime}: {checks} coefficients passed")
    print(f"total: {total} multivariate coefficient congruences passed")


if __name__ == "__main__":
    main()
