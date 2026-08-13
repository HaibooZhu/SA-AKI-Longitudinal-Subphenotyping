# Reproducible literature-search record

- Search date: 2026-08-09
- Database: PubMed (NCBI E-utilities)
- Unique records retrieved: 303
- Screening status: key overlapping records screened; full retrieved set retained for audit
- Positioning decision: remove any unqualified claim that this is the first SA-AKI subphenotyping study

## Exact queries

### sa_aki_trajectory_subphenotype

```text
("sepsis-associated acute kidney injury"[Title/Abstract] OR "septic acute kidney injury"[Title/Abstract] OR (sepsis[Title/Abstract] AND "acute kidney injury"[Title/Abstract])) AND (trajectory[Title/Abstract] OR trajectories[Title/Abstract] OR subphenotype*[Title/Abstract] OR phenotype*[Title/Abstract]) AND ("2023/01/01"[Date - Publication] : "3000"[Date - Publication])
```

Hits: 142; retrieved: 142

### critical_illness_aki_trajectory

```text
("acute kidney injury"[Title/Abstract]) AND (trajectory[Title/Abstract] OR trajectories[Title/Abstract] OR subphenotype*[Title/Abstract]) AND (critical care[Title/Abstract] OR critically ill[Title/Abstract] OR intensive care[Title/Abstract]) AND ("2023/01/01"[Date - Publication] : "3000"[Date - Publication])
```

Hits: 128; retrieved: 128

### adqi_28_sa_aki_subphenotype

```text
("sepsis-associated acute kidney injury"[Title/Abstract] OR (sepsis[Title/Abstract] AND "acute kidney injury"[Title/Abstract])) AND (subphenotype*[Title/Abstract] OR phenotype*[Title/Abstract] OR trajectory[Title/Abstract] OR trajectories[Title/Abstract]) AND (ADQI[Title/Abstract] OR "Acute Disease Quality Initiative"[Title/Abstract] OR "ADQI 28"[Title/Abstract])
```

Hits: 3; retrieved: 3

### furosemide_stress_or_response_aki

```text
("acute kidney injury"[Title/Abstract]) AND ("furosemide stress test"[Title/Abstract] OR "furosemide responsiveness"[Title/Abstract] OR "diuretic responsiveness"[Title/Abstract])
```

Hits: 70; retrieved: 70

## Eligibility criteria

Include adult ICU/critical-care AKI or SA-AKI trajectory/subphenotype studies and directly relevant consensus or furosemide-response methods papers. Exclude pediatric-only, non-human, biomarker-only without trajectory relevance, conference-only, protocol-only, and preprint records.

The CSV preserves all retrieved records. Citation decisions must be documented separately after relevance screening; retrieval alone does not establish eligibility.

## Positioning implication

The updated search identified adult SA-AKI subphenotyping and sepsis-AKI creatinine-trajectory studies published before this revision, including PMID 38445412, PMID 38730421, and PMID 41557579. The manuscript must therefore avoid an unqualified first-ever novelty claim. The defensible contribution is narrower: cross-database reproducibility of three multivariable kidney-function trajectories over the first seven days after SA-AKI onset, with explicit differences in available inputs across cohorts.

Key inclusion/exclusion decisions are recorded in `key_literature_screening.csv`.
