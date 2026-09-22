# Building the JOSS draft PDF

## GitHub Action

A workflow under `.github/workflows/joss-draft.yml` compiles `paper/paper.md`
with the Open Journals `inara` image on each change to `paper/`.

## Local Docker

```bash
docker run --rm \
  --volume "$PWD/paper:/data" \
  --user "$(id -u):$(id -g)" \
  --env JOURNAL=joss \
  openjournals/inara
```

On Windows (PowerShell), mount the absolute `paper` path and omit `--user` if
needed.

## Before submission checklist

- [ ] Add ORCID to `paper/paper.md` author metadata
- [ ] Confirm affiliation wording
- [ ] Upload Zenodo bundle; paste DOI into `paper.bib` (`cascade_audit2026`),
      `CITATION.cff`, and optionally the paper text
- [ ] Fix any `orcid: # ADD...` placeholder (YAML must be valid — remove the
      comment line or set a real ORCID)
- [ ] Open a JOSS pre-submission inquiry / submit at https://joss.theoj.org/
- [ ] Ensure the GitHub repo is public with Apache-2.0 `LICENSE`
- [ ] Point reviewers at `README.md`, `scripts/smoke_gpu.sh`, and `pytest tests/`
