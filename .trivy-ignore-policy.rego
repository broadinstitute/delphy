#
# .trivy-ignore-policy.rego
#
# Rego policy for class-based CVE filtering in Trivy security scans.
# This is a skeleton/no-op policy that can be filled in later if scan
# noise requires systematic filtering.
#
# CONTEXT:
#   Delphy is a C++20 Bayesian phylogenetics tool that:
#   - Runs as a command-line tool (no network listeners, no web UI)
#   - Processes genomic data files (FASTA, DPHY format)
#   - Is distributed both as native binaries and Docker containers
#   - Uses dependencies from Conan (boost, eigen, ms-gsl), git submodules
#     (abseil, flatbuffers, ctpl, cxxopts), and vendored headers (cppcoro)
#
# USAGE:
#   The security-scan.yml workflow automatically passes this policy to
#   Trivy via --ignore-policy. You do not need to manually specify it.
#
# TO ADD FILTERS:
#   Uncomment and modify the example rules below, or add new ones.
#   Common patterns for bioinformatics/computational tools:
#   - Filter CVEs requiring physical access (AV:P)
#   - Filter CVEs requiring adjacent network access (AV:A)
#   - Filter CVEs requiring local + user interaction (AV:L + UI:R)
#   - Filter CVEs with availability-only impact (C:N + I:N)
#
# CVSS VERSION SUPPORT:
#   Rules should handle both CVSS v3.1 and v4.0 vector strings since
#   Trivy is transitioning to v4.0 for newer advisories.
#
# LAST REVIEWED: 2026-05-22
# REVIEW CADENCE: As needed when false positives accumulate

package trivy

default ignore = false

# Example rule (currently disabled - uncomment and modify as needed):
#
# ignore {
#   # Physical access required (AV:P) - cloud-hosted binaries are never physically accessible
#   input.CVSS != null
#   contains(input.CVSS, "AV:P")
# }
