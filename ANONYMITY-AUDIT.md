# Anonymity audit

## Classification

This is a **content-anonymous** release.  The manuscript and payload omit the
author's name, affiliation, email address, ORCID, acknowledgments, submission
identifiers, correspondence, and local filesystem paths.  The GitHub account
and commit metadata may nevertheless be linked to the author, so this release
does not claim author-unlinkable double blindness.

## Audit surface

The release audit covers the visible PDF, extracted PDF text, PDF metadata,
TeX source and comments, Python programs, Markdown files, manifest, filenames,
Git tree, and archive/local-path patterns.  The PDF was rendered page by page
for visual inspection.  The source package contains no inherited `.git`
directory, cover letter, referee correspondence, submission-portal metadata,
or local build auxiliaries.

The GitHub web interface confirms that the repository is public and displays
the intended tree and public commit.  A post-push fresh clone was unavailable
because outbound terminal connections to `github.com:443` timed out; remote
byte-for-byte verification is therefore partial for this version.
