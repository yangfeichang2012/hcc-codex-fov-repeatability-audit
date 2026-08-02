# Source data

This repository does not redistribute the public source archives.

The deterministic reconstruction uses three files from Zenodo record `19123188`:

| File | Bytes | SHA-256 | Role |
|---|---:|---|---|
| `metadata.csv` | 3,556 | `9e87b7c374db5b06de3bc22e2a1c295232cf61e964c2b0548cf1c9fc3121b882` | Patient, state, response, and treatment mapping |
| `CODEX-MIF.zip` | 863,912,545 | `afd294100a878f51ecac6f7e79be3a191f8a44cd7440f049dbcae1bc134c8373` | Forty-eight source cell-type tables; expression tables are not extracted |
| `notebook_and_resources_v2.zip` | 29,863,743 | `f3fa9ab603709522c5e6eee2c112be267720e66e6e1c4075fa2d1db3ed351750` | Source-fixed microenvironment annotations |

Record DOI: `10.5281/zenodo.19123188`  
Record licence: CC BY 4.0  
Source citation: Wu et al., *Spatial multi-omics and deep learning reveal fingerprints of immunotherapy response and resistance in hepatocellular carcinoma*.

The source record contains de-identified public identifiers. The release pipeline keeps those identifiers in memory only and writes anonymized leave-one-patient-out slots (`P01`, `P02`, ...). It writes no region-level ratios or cell-level records.
