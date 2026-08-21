# Known limitations

1. The numerical PTM RMSD cutoff and the associated base configuration are unavailable in the retained project record. No default value was inferred.
2. The initial interface gap was not archived. No value was reconstructed.
3. Retained full binary restart/final-state dump payloads are stored only in `zenodo_archive/` and are excluded from the GitHub upload set. A continuous atom-by-atom dump series is not claimed where it was not present in the audited source.
4. `M4_LIT` is a single-trajectory literature-like topology reference, not a controlled n=3 model.

These limitations must remain visible in any public release and in the final Data Availability wording.
