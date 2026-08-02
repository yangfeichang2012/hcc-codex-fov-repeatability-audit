# Within-biopsy CODEX FOV repeatability audit

This archive reconstructs the three frozen SPARC-MIF ratios from the exact public source objects in Zenodo record `19123188`, then reproduces the patient-level ICC(1,1), 2,000-replicate patient-cluster bootstrap intervals, and leave-one-patient-out ICC diagnostics reported for PRJ-2026-009.

## Scope

- Population: the source-selected HCC spatial subset in the public record.
- Primary set: pretreatment CODEX observations with two available FOVs per patient.
- Biological unit: patient; cells and FOVs do not increase the inferential sample size.
- Features: SPARC-MIF-01, SPARC-MIF-02, and SPARC-MIF-03 only.
- Main estimator: balanced one-way random-intercept ANOVA ICC(1,1), with negative values truncated to zero as frozen before analysis.
- Uncertainty: 2,000 patient-cluster bootstrap replicates, seed `20260729`, percentile 95% intervals.

The archive does not analyse the seven non-executable transcriptomic programs, the technical benchmark numerically, or the CAR-T cohort numerically. It does not train a prediction model and does not establish biopsy-to-biopsy, whole-tumour, mechanistic, or clinical reliability.

## Reproduction

Python 3.12 is recommended.

```bash
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements.txt
.venv/Scripts/python code/run_source_to_result.py --download --output-dir reproduced_outputs
.venv/Scripts/python code/validate_release.py --output-dir reproduced_outputs
```

On Linux/macOS, use `.venv/bin/python` instead. The `--download` route retrieves only three immutable public source objects and verifies their byte counts and SHA-256 hashes before processing. To use a pre-existing source cache, replace `--download` with `--source-dir PATH`.

The analysis extracts only the 48 cell-type tables and one source microenvironment-annotation object into a temporary directory. No expression matrix is extracted or read. No cell-level, patient-level, FOV-level, or paired-ratio file is written to the release outputs.

## Expected outputs

- `primary_icc.tsv`: the three primary point estimates and intervals.
- `bootstrap_icc_draws.tsv`: all 6,000 aggregate ICC bootstrap draws used for the intervals.
- `leave_one_patient_out.tsv`: anonymized patient-slot influence diagnostics.
- `source_preflight.json`: verified source hashes and structural checks.
- `analysis_manifest.json`: execution boundary and result summary.
- `SHA256SUMS.txt`: hashes of deterministic outputs other than itself.

The frozen expected outputs are under `expected/`. `validate_release.py` compares a fresh run byte-for-byte with that directory and checks the registered numerical claims.

## Data provenance and licensing

The source objects remain hosted by their creators in [Zenodo record 19123188](https://zenodo.org/records/19123188) under CC BY 4.0. They are not redistributed here. See `SOURCE_DATA.md`, `NOTICE`, and `LICENSE_DATA_AND_OUTPUTS.md`.

Code is licensed under BSD-3-Clause. Derived aggregate tables, documentation, and figures are licensed under CC BY 4.0 unless a file states otherwise.

## Citation

Use the version DOI in `CITATION.cff` for this archive and cite the source dataset separately.
