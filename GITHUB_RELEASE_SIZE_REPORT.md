# GitHub release size report

- Recommended GitHub upload set: every file in this directory except `zenodo_archive/`.
- The final measured GitHub-layer file count and byte size are recorded in `RELEASE_SIZE_METRICS.json` to avoid self-referential report-size drift.
- `.gitignore` excludes `zenodo_archive/`, restart files, dump files, and generic binary files.

## Large files in GitHub layer

- None at or above 10 MiB.

## Excluded file groups

- full model data: Zenodo only; GitHub has metadata/hash: 4 files, 141713550 bytes (135.15 MiB)
- large complete figure source: Zenodo only: 37 files, 71331007 bytes (68.03 MiB)
- large final/restart state: Zenodo only: 84 files, 3839143154 bytes (3.58 GiB)
- raw trajectory: Zenodo only: 10 files, 5394997 bytes (5.15 MiB)

Full per-file exclusions are listed in `GITHUB_EXCLUDED_FILES.csv`.
