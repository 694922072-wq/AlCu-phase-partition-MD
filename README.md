# Al4Cu9/Al2Cu phase-partition sensitivity in Cu/Al interfacial multilayers

Data and code release for the associated Computational Materials Science manuscript, “Al4Cu9/Al2Cu phase-partition sensitivity in Cu/Al interfacial multilayers: a molecular dynamics study of load-bearing and candidate localization.”

## Purpose

This release provides an auditable route to the retained simulation inputs, nine controlled mechanical trajectories, processed mechanical statistics, figure/supplementary source data, analysis code, model structures, and archived state files. It preserves the manuscript's evidence boundaries and does not claim to reproduce unavailable atomistic settings.

## Models

- **M3_SYM** — Cu/Al2Cu/Al/Al2Cu/Cu; thick-Al2Cu symmetric control.
- **M4_RATIO** — Cu/Al4Cu9/Al2Cu/Al/Al2Cu/Al4Cu9/Cu; phase-partition bridge.
- **M4_SYM** — Cu/Al4Cu9/Al2Cu/Al/Al2Cu/Al4Cu9/Cu; symmetric composite-IMC model.
- **M4_LIT** — Cu/Al4Cu9/Al2Cu/Al; literature-like topology reference outside the controlled n=3 chain.

The public manuscript label is `M4_RATIO`; see `MODEL_NAME_MAPPING.md` for its historical internal alias.

## Simulation workflow

```text
construction -> minimization -> 300 K equilibration -> tension -> analysis
```

## Reproducibility scope

The controlled mechanical comparison uses nine independent velocity-seed trajectories: three each for M3_SYM, M4_RATIO, and M4_SYM. The seeds and raw SHA-256 values are listed in `processed_data/AUTHORITATIVE_NINE_TRAJECTORIES.csv`. M4_LIT is outside this n=3 chain.

## Folder description

- `code/` — model construction, mechanical extraction, candidate-overlap, cluster, PTM, DXA, and plotting scripts.
- `inputs/` — unmodified LAMMPS inputs organized under public model names.
- `model_structures/` — lightweight model metadata, crystallographic inputs, and Zenodo links/hashes.
- `processed_data/` — final mechanical metrics and authoritative trajectory index.
- `figure_source_data/` — lightweight per-figure source tables, metadata, and plotting scripts.
- `supplementary_source_data/` — source tables supporting the Supplementary Information.
- `zenodo_archive/` — complete expanded payload and the ZIP64 archive; excluded from GitHub by `.gitignore`.

## Data availability

The GitHub layer and Zenodo archive are technically prepared but have not been uploaded by this task. After deposition, replace `[GITHUB_URL]` and `[ZENODO_DOI_OR_RECORD_URL]` in `DATA_AVAILABILITY_STATEMENT.md`. No URL or DOI has been invented.

## Potential

The third-party potential `CuAgAuNiPdPtAlPbFeMoTaWMgCoTiZr_Zhou04.eam.alloy` is not redistributed. Obtain it from the NIST Interatomic Potentials Repository and verify SHA-256 `cf667915dcf1327d3be5379e59005c82e2ee22c6379de2d105e2c88467217a1c`.

## Licenses

- Original code and workflow scripts: MIT (`LICENSE_CODE.txt`).
- Original data and data documentation: CC BY 4.0 (`LICENSE_DATA.txt`).
- Third-party potential: excluded and governed by its source terms.

## Citation

Use `CITATION.cff`. Add the final article DOI and Zenodo DOI after assignment.

## Limitations and audit

Read `KNOWN_LIMITATIONS.md` and `PUBLIC_RELEASE_AUDIT_REPORT.md` before release. Technical audit status: **PASS**.
