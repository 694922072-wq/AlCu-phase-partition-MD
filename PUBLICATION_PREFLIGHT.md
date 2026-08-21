# Publication Preflight

- Overall: **PASS**
- Immutable source: `D:\COMMAT_PUBLIC_RELEASE_v1`
- Scope: validation only; no scientific payload was modified or recomputed.

## Checks

| Check | Status | Evidence |
|---|---|---|
| Release source exists | **PASS** | D:\COMMAT_PUBLIC_RELEASE_v1 |
| Manifest coverage and SHA256 | **PASS** | 385 rows verified; manifest self excluded |
| Nine authoritative trajectories | **PASS** | 9 model/replica/seed/raw-file SHA256 rows verified |
| Zenodo archive identity | **PASS** | 2233375981 bytes; SHA256 9a28974a4011641e8eae1797034a4ce7f34b0c4b08d2250075ec92eefaec43cb |
| GitHub lightweight-layer exclusions | **PASS** | 124 files; 753571 bytes; no Zenodo archive, state dump, potential binary, or >100 MiB file |

## Failures

- None.

## Scientific immutability

PASS. This preflight performed read-only hashing, metadata inspection, and release-layer classification only.
