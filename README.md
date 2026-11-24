# bsd-compute

Free GitHub Actions compute farm for the BSD (Birch–Swinnerton-Dyer) campaign.

The campaign's local box is small and saturated by its own pipelines
(Phase 1 certification, congruent sieve). This public repository exists so
CPU-bound, embarrassingly parallel PARI/GP batch work can be offloaded to
GitHub-hosted runners — free for public repositories.

## Usage

Actions → **bsd-compute** → Run workflow → choose a mode:

- `proof` — proof of life: E_34 (y² = x³ − 1156x) conductor, rank,
  a_p table (p ≤ 50), Tamagawa artefact check; also measures and prints
  the runner's actual CPU/RAM/disk from inside the job.
- `batch` — the real payload shape: for the 20 certified curves in
  `curves.json`, compute a_p tables (p ≤ 500) and analytic Sha via
  `ellanalyticrank`, and report throughput in curves/hour.

Results land as a job artifact and in the job summary. No secrets are
used or needed. Nothing runs on push — manual dispatch only.

## Layout

- `.github/workflows/bsd-compute.yml` — the workflow
- `scripts/proof_e34.gp` — proof-of-life GP script
- `scripts/batch.py` — batch payload (cypari2, GP-string evaluation)
- `curves.json` — batch input: label, certified rank, minimal
  a-invariants, conductor (from the campaign's certified set)
