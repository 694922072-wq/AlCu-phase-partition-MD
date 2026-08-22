# Historical pre-publication audit notice

This file records the pre-publication release audit completed before
the GitHub and Zenodo records were made public.

The public release has since been completed successfully.

Current public identifiers:

GitHub:
https://github.com/694922072-wq/AlCu-phase-partition-MD

Zenodo:
https://doi.org/10.5281/zenodo.22043296

For the current publication state, refer to the repository README,
CITATION.cff, DATA_AVAILABILITY_STATEMENT.md, and public Zenodo record.

---

# Public release audit report

- Audit timestamp: 2026-08-20T15:20:46.221115+00:00
- Scientific immutability records: 345
- Zenodo payload inventory rows: 258
- Final measured release, GitHub-layer, expanded-payload, and ZIP sizes: see `RELEASE_SIZE_METRICS.json`
- Overall technical gate: PASS

| Check | Status | Evidence |
|---|---|---|
| Scientific files unchanged | PASS | All copied scientific payloads match their source SHA-256 values; only public-facing metadata and paths were generated. |
| Nine trajectories preserved | PASS | Authoritative controlled raw trajectories present: 9. |
| Seed consistency | PASS | Seeds match the accepted nine-trajectory list. |
| License status | PASS | MIT for code/workflows and CC BY 4.0 for original data; scopes are separated. |
| Historical public-path cleanup | PASS | Prohibited public path components found: 0. |
| M4_RATIO public naming | PASS | Public directory and metadata label are M4_RATIO; historical internal label is documented. |
| GitHub readiness | PASS | GitHub layer contains no restart/dump/binary state payload; Zenodo directory is ignored. |
| Zenodo inventory | PASS | Inventory rows: 258; bad rows: 0. |
| Zenodo models | PASS | Four full construction-ready model data files are present. |
| Zenodo binary states | PASS | Retained binary/final state files: 84; final data/restart pair present for 9/9 trajectories. |
| Zenodo archive CRC | PASS | ZIP64 archive created and CRC-tested. |
| Third-party potential excluded | PASS | Potential name/hash/source are documented; bytes are absent. |

## Remaining author actions

1. Create the GitHub repository and upload every release-root file except `zenodo_archive/`; confirm `.gitignore` is honored.
2. Create a Zenodo draft, upload `zenodo_archive/COMMAT_D_26_03451_ZENODO_v1.zip`, and reserve a DOI.
3. Insert the GitHub URL and Zenodo DOI/reviewer link into README/Data Availability and the CMS revision materials.
4. Confirm that all authors/institution have authority to apply MIT and CC BY 4.0, then approve the release tag/version.
5. Optionally connect the Zenodo record to the GitHub release after the repository is public.

## Release decision

**PASS: the package can enter the GitHub/Zenodo publication stage, subject to the author actions above.**

No LAMMPS calculation, seed change, stress-strain edit, raw-trajectory edit, model regeneration, or manuscript scientific revision was performed.
