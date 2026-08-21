# Reproducibility guide

## 1. Verify the release

Verify every file listed in `MANIFEST_SHA256.csv`. The manifest excludes itself. For the Zenodo payload, also verify `zenodo_archive/ZENODO_FILE_INVENTORY.csv`.

## 2. Obtain the potential

Download `CuAgAuNiPdPtAlPbFeMoTaWMgCoTiZr_Zhou04.eam.alloy` from the NIST Interatomic Potentials Repository and verify SHA-256 `cf667915dcf1327d3be5379e59005c82e2ee22c6379de2d105e2c88467217a1c`. Atom type 1 is Al and atom type 2 is Cu.

## 3. Select a model

Full construction-ready data files are in `zenodo_archive/payload/model_structures/`. Public model names are M3_SYM, M4_RATIO, M4_SYM, and M4_LIT. See `MODEL_NAME_MAPPING.md` for the internal historical M4_RATIO alias retained in unmodified source text.

## 4. Run workflow

```text
construction
  -> minimization
  -> 300 K equilibration (NPT then NVT)
  -> tension
  -> analysis
```

The LAMMPS files in `inputs/` are unmodified scientific inputs. Their original relative `read_data` and output paths are preserved. Prepare a separate run directory with the expected `structure/`, `outputs/`, `logs/`, and `dumps/` locations rather than editing the archived release copy.

## 5. Independent seeds

The controlled mechanical dataset contains nine independent Maxwell-Boltzmann velocity-seed trajectories: three each for M3_SYM, M4_RATIO, and M4_SYM. M4_LIT remains a separate single-trajectory topology reference.

## 6. Analysis

Use `code/SCRIPT_INDEX.md` to identify mechanical, candidate-overlap, cluster, PTM, DXA, and plotting scripts. Some scripts retain original absolute workstation paths; configure paths in a separate working copy. Do not overwrite released files.

## 7. Interpretation boundary

This package supports inspection and rerunning within the archived scope, but it does not claim to reproduce unavailable PTM settings or an unarchived initial interface gap. Atomistic descriptors must not be upgraded from representative candidate/indicator evidence to n=3 mechanistic proof.
