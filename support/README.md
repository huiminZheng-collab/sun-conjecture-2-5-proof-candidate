# Exact verification programs

These four programs provide exact-integer checks for the transformations,
valuation bounds, paired multinomial terms, and finite coefficient regimes
reported in the manuscript.  They are independent error-detection aids and
are not part of the logical proof.

## Requirements

- Python 3.10 or later
- no third-party Python packages

## Run

From this directory, execute:

```text
python c25_constant_term_verify.py
python c25_multivariate_verify.py
python referee_carry_lemma_check.py
python referee_jacobsthal_pairing_check.py
```

Every program terminates by assertion failure if a checked identity or
congruence fails.

## Coverage

- Three independent formulas for the sequence through `n=30`.
- The four-factor Laurent identity and signed balanced-array expansion.
- Five sharp primitive-sector valuation regimes.
- 2,508 lifted multinomial Jacobsthal ratios.
- 39 direct instances of Sun's congruence.
- 368 simultaneously scaled multivariate coefficients.
- 29,260 Legendre digit-level identities and 20 sharpness witnesses.
- 5,928 multinomial Jacobsthal cases.
- 510,768 paired divisible-array cases.
- 840 complete coefficient cases, including zero and boundary indices.

All arithmetic is exact and no random sampling is used.

