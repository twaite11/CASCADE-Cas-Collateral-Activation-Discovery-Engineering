# Building the JOSS draft PDF

## Figures

Regenerate plots from local mining + RL outputs:

```bash
python scripts/build_paper_figures.py
```

Outputs land in `paper/figures/` (PNG plots + optional SVGs copied from
`Desktop/patent_figures` when present). Edit `paper/paper.md` captions as needed.

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

## Before submission checklist

- [ ] Add ORCID in `paper/paper.md`
- [ ] Confirm affiliation
- [ ] Paste Zenodo DOI into `CITATION.cff` / `paper.bib` when available
- [ ] Skim figures: only the three *L. booriae* Cas13a are claimed as verified
- [ ] Submit at https://joss.theoj.org/
