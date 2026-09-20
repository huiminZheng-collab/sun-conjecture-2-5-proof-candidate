# Proof of Z.-H. Sun's Conjecture 2.5

This content-anonymous release contains a manuscript proving Conjecture 2.5
in Z.-H. Sun, *New congruences involving Apéry-like numbers*,
arXiv:2004.07172v2.  In fact, the manuscript proves a multivariate
coefficient supercongruence that implies the conjecture.

## Main result

For the sequence

```text
G_n = sum_{k=0}^n 4^k binom(2n-2k,n-k)^2 binom(2k,k),
```

the paper proves, for every odd prime `p` and all positive integers `m,r`,

```text
G_(m p^r) = G_(m p^(r-1))  (mod p^(2r)).
```

## Contents

- `paper/main.pdf`: content-anonymous manuscript.
- `paper/main.tex`: sanitized manuscript source.
- `support/*.py`: exact-integer diagnostic checks used during verification.
- `support/README.md`: requirements, commands, and coverage of the checks.
- `RELEASE-MANIFEST.sha256`: SHA-256 hashes of the frozen public payload.
- `VERIFY.ps1`: fail-closed manifest verifier.
- `ANONYMITY-AUDIT.md`: release-level anonymity statement.

The programs require Python 3.10 or later and use only the standard library.

## Status and provenance

Version 0.3 is a writing-focused revision dated 20 September 2026.  It
reorganizes the introduction around the source of the second power of p,
states the relation to existing Dwork and multivariate frameworks more
precisely, and adds signposts at the two decisive valuation steps.  The
theorem and proof mechanism are unchanged.  This repository is public
dissemination, not a journal submission.

The public payload was assembled in a new staging directory, without the
submission correspondence, cover letter, portal metadata, or inherited Git
history.  The final public tree was scanned for identity strings and local
paths, and the PDF was separately checked for visible identity and metadata.

## AI-use disclosure

The proof presented in this article was found by OpenAI Codex.

## Non-claims

Public availability does not imply acceptance, peer review, correctness, or
priority.  The manuscript omits the author's identity, but the hosting account
and public Git history may be identity-linked.  This is therefore a
**content-anonymous** release, not an author-unlinkable double-blind release.
