# Constraint Lease Audit Protocol

**A freeze-first, mechanically-mapped audit protocol for identifiability of published computational models.**

[![DOI: pending](https://img.shields.io/badge/DOI-pending-blue)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Companion repository for the manuscript submitted to *Bioinformatics* (under review).

## What is this?

Published computational biology models are heavily reused, but the provenance of
their parameters is rarely examined systematically: which values are experimentally
anchored, which are borrowed across species, which are fitted constructs, and which
are outright assumptions?

This repository provides:

1. **A freeze-first audit protocol** — criteria are frozen before any distribution
   statistic is computed; per-item classification is then a mechanical table lookup
   (zero manual adjudication across 1,217 judgements).
2. **Per-item judgement tables** for six published models spanning whole-cell ODE,
   FBA, and agent-based types:
   - Karr et al. 2012 whole-cell model (154 parameters)
   - O'Hara–Rudy 2011 cardiac AP model (262 parameters)
   - Wodke et al. 2013 *M. pneumoniae* metabolism (306 reactions)
   - Gaspari et al. 2020 serum-free adaptation (489 reactions)
   - E-Cell 1999 hypothetical cell (3 model-level entries)
   - BAAIWorm 2024 *C. elegans* (3 layer-level entries)
3. **A gap-morphology library** — eight recurring gap forms (F1–F8) inventoried
   across all models.
4. **Schema v1.5** — a four-state uniqueness taxonomy:
   `identified` / `one-of-many` / `non-identified` / `paradigmatic`
   (the fourth state covers deep-learning weights: intentional non-identifiability
   by training paradigm, functionally equivalent and reproducible — distinct from
   *unintended* non-identifiability i.e. source debt).

## Repository layout

```
schema/       constraint_census_schema v1.4 & v1.5 (four-state)
judgements/   per-item judgement tables + verdicts (six models)
demo/         consumer-side demo + figure pipeline
figures/      manuscript figures (PNG)
```

## Quick start

```bash
python3 demo/lease_consumption_demo.py   # consumer-side validation
python3 schema/constraint_census_schema_v15.py  # schema self-test
```

## Key findings (see manuscript)

- Identifiability profiles differ by model type: experimentally anchored tables
  (Karr 58.4% identified) vs constraint-based networks (Wodke 51.0%) vs
  hypothetical cells (E-Cell 33.3%).
- Gap forms follow a type-consistent methodology fingerprint.
- **Model repositories do not currently carry the maintenance-history structure
  that identifiability auditing requires** — a reportable property of the
  ecosystem and a governance opportunity.

## Extending the audit (v1.5 four-state)

The four-state taxonomy has been applied beyond ODE/FBA to logic models
(Saez-Rodriguez 2007), agent-based frameworks (PhysiCell 2018), discriminative
deep learning (DeepBind 2015), structure-prediction transformers (AlphaFold2),
and generative diffusion models (RFdiffusion) — see `ten`/`eleven`-model
comparison archives in the manuscript supplementary data.

## Citation

If you use this protocol or the judgement tables, please cite the manuscript
(citation pending publication) and this repository.

## License

MIT — see [LICENSE](LICENSE).
