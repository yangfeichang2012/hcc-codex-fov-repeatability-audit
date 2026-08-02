# Protocol history and analysis boundary

## Frozen design

The internal protocol and SAP were frozen on 29 July 2026 before feature-value reading. The result-blind feature-definition audit found three executable MIF ratios and seven non-executable transcriptomic programs. The prespecified 1-4-feature branch therefore narrowed the analysis to per-feature reliability estimates and removed the proposed `stable_feature_fraction` summary.

The frozen primary analysis used pretreatment CODEX observations, patient as the biological unit, two available FOVs per eligible patient, ICC(1,1), 2,000 patient-cluster bootstrap replicates, seed `20260729`, and leave-one-patient-out influence diagnostics.

## Review-driven presentation correction

External review identified a part-whole dependency in the earlier selected-FOV-versus-two-FOV-mean rank correlation. That statistic and all related qualification language were retired. The current manuscript and this archive use ICC(1,1) as the primary estimator. The scientific scope is within-biopsy paired-FOV repeatability, not biopsy-to-biopsy or whole-tumour reliability.

## Source-to-result route

Version 1.0.0 implements the author-approved deterministic public-source reconstruction route. It starts from the three hash-pinned objects in Zenodo record `19123188`, extracts no expression matrix, reconstructs the three ratios in memory, regenerates the full 6,000-draw bootstrap series, and writes only aggregate or anonymized influence outputs.

The 20-to-13 source-selection history, block/section provenance, FOV distances, and patient-level pathological response reconstruction remain unavailable from the released objects. All estimates remain conditional on the source-selected spatial subset, paired-FOV availability, and feature estimability.
