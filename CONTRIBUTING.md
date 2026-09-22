# Contributing to CASCADE

Thanks for interest in CASCADE. This project is computational methods software
(Apache-2.0). Please read [docs/DUAL_USE.md](docs/DUAL_USE.md) before proposing
features that expand biological capability beyond the documented intended use.

## Development

1. Fork and clone the repository.
2. Create a Python 3.11+ environment and install `requirements.txt` plus `pytest`.
3. Run CPU tests: `pytest tests/ -q` (no GPU required).
4. Do not commit secrets (`.env`), model weights, or large raw `outputs/` trees.
   Use `scripts/build_zenodo_bundle.py` for archival deposits.

## Pull requests

- Prefer small, reviewable PRs with a clear problem statement.
- Add or update tests for behavior changes in fitness, stitch, mining filters,
  or API surfaces.
- Keep marketing claims out of the README; application language belongs in the
  unvalidated application / dual-use docs.
- Update `THIRD_PARTY_NOTICES.md` if you add a runtime dependency that ships
  or downloads weights.

## Issues

Use GitHub Issues for bugs and feature requests. Include OS, GPU (if any),
`EVAL_CMD`, and a minimal reproduction.

## Code of conduct

Be respectful. Harassment or attempts to solicit harmful biological misuse
will result in blocked collaboration.
