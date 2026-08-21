# Data scope

## GitHub layer

The GitHub-ready layer consists of all files outside `zenodo_archive/`. It contains original code, unmodified LAMMPS inputs, lightweight model metadata, processed mechanical tables, figure-level processed source data and plotting scripts, supplementary source tables, documentation, and audit files. It excludes raw trajectories, full construction-ready model data, restart files, and large figure assets.

## Zenodo layer

`zenodo_archive/` contains an expanded deposition payload and a ZIP archive. The payload contains exactly nine authoritative controlled raw mechanical trajectories, the M4_LIT topology-reference trajectory, four full construction-ready model data files, retained final/equilibrated/binary restart states, complete figure source data, processed and supplementary source data, code, inputs, and public metadata.

## Third-party potential

The file `CuAgAuNiPdPtAlPbFeMoTaWMgCoTiZr_Zhou04.eam.alloy` is not redistributed. Obtain it from the NIST Interatomic Potentials Repository and verify SHA-256 `cf667915dcf1327d3be5379e59005c82e2ee22c6379de2d105e2c88467217a1c`.

## Evidence boundary

The n=3 replication applies only to global mechanical metrics. Candidate localization, candidate-cluster, PTM, and DXA results remain representative-trajectory evidence under the present potential and loading protocol.
