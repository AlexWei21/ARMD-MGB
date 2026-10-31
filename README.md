# ARMD-MGB

The Antibiotic Resistance Microbiology Dataset – Mass General Brigham (ARMD-MGB) is a de-identified resource for antimicrobial resistance (AMR) research derived from electronic health records. It covers more than 10 years of blood, urine and respiratory cultures from adult patients across the Mass General Brigham healthcare system. Each culture is linked to its antibiotic susceptibility results and to the clinical context in which it was collected. That context includes prior antibiotic exposure, earlier cultures and resistant organisms, comorbidities, procedures, laboratory results, vital signs, care setting and neighborhood deprivation.

Resistance models and risk factors estimated in one health system can depend on its local prescribing, testing and documentation practices. ARMD-MGB provides an independently collected population with a culture-centered structure shared with the Stanford ARMD dataset. This makes it possible to develop and evaluate methods for antimicrobial stewardship, causal inference and clinical decision support across health systems.

The dataset is available on [PhysioNet](https://physionet.org/content/armd-mgb/). This repository contains the code that builds each released table from MGB source data.

---

## Repository Layout

```
scripts/
├── microbiology_cultures_*.ipynb      # Cohort, demographics, care setting, ADI scores
├── time-to-event-augmented-queries/   # Prior clinical history features
├── labs.{py,R}, vitals.{py,R}         # Laboratory and vital-sign summaries
├── implied_susceptibility.py          # Rule-based implied susceptibility
├── prior_abx_class_subtype_lookup.ipynb
└── post_pulling/                      # Date shifting, deduplication and conflict filtering
```

---

## Released Tables

Script paths are relative to `scripts/`.

### Microbiology Cultures

| Table | Script | Description |
|-------|--------|-------------|
| `microbiology_cohort.csv` | `microbiology_cultures_cohort_query.ipynb` | Defines the study cohort. Selects blood, urine and respiratory cultures from patients older than 18 at admission and links each culture to the encounter during which it was collected. A culture is kept if it is the patient's first, or if it was collected more than 14 days after the patient's previous culture. Susceptibility results are collapsed into Susceptible, Intermediate and Resistant. The notebook also assigns the de-identified patient, encounter and order IDs and each patient's random date offset. |
| `implied_susceptibility.csv` | `implied_susceptibility.py` | Infers results for antibiotics that were not tested on a culture. It uses intrinsic-resistance rules and results for related antibiotics on the same culture and organism (for example, doripenem from meropenem in *Acinetobacter*). Covers *Acinetobacter*, *Citrobacter*, *Enterobacter*, *Escherichia coli*, *Klebsiella*, *Proteus*, *Pseudomonas aeruginosa* and *Serratia*. Inferred values are kept in a separate `implied_susceptibility` column. |

### Patient and Encounter Context

| Table | Script | Description |
|-------|--------|-------------|
| `demographics.csv` | `microbiology_cultures_demographics.ipynb` | Age at culture collection in groups from 18–24 years to above 90, and sex (`gender`: Female, Male, Others). |
| `ward_type.csv` | `microbiology_cultures_ward_info.ipynb` | Care setting at the time of the culture: inpatient, outpatient, emergency department and ICU indicators. |
| `adi_scores.csv` | `microbiology_cultures_adi_scores.ipynb` | Area Deprivation Index score and state rank for the patient's 9-digit ZIP code (Neighborhood Atlas). |

### Prior Clinical History

Scripts are in `time-to-event-augmented-queries/`.

| Table | Script | Description |
|-------|--------|-------------|
| `prior_abx.csv` | `medication_exposure.ipynb` | Antibiotics ordered before the culture, with drug category and time from order to culture. |
| `prior_org.csv` | `previous_infecting_organisms.ipynb` | Organisms grown in the patient's earlier cultures, grouped into categories, with days before the current culture. |
| `microbial_resistance.csv` | `microbial_resistance.ipynb` | Resistant results from the patient's earlier cultures, by organism and antibiotic, with time before the current culture. |
| `comorbidity.csv` | `comorbidities.ipynb` | ICD-10 diagnosis codes mapped to Elixhauser comorbidity and AHRQ CCSR categories. |
| `prior_procedures.csv` | `prior_procedures.ipynb` | Prior mechanical ventilation, surgery, central venous catheter, dialysis and urinary catheter, with days before the culture. |
| `nursing_home_visits.csv` | `nursing_home_visits.ipynb` | Prior nursing home or skilled nursing facility contact, identified from visit codes, admission source and discharge disposition. |

### Laboratory Results and Vital Signs

| Table | Script | Description |
|-------|--------|-------------|
| `labs.csv` | `labs.py` (R version: `labs.R`) | 25th and 75th percentiles, median, first and last value over the 14, 30 and 180 days before the culture. Labs: WBC, neutrophils, lymphocytes, hemoglobin, platelets, sodium, bicarbonate, BUN, creatinine, lactate and procalcitonin. |
| `vitals.csv` | `vitals.py` (R version: `vitals.R`) | The same summaries over the 48 hours before the culture for heart rate, respiratory rate, temperature, and systolic and diastolic blood pressure. |

### Others

| Script | Description |
|--------|-------------|
| `prior_abx_class_subtype_lookup.ipynb` | Builds a lookup table mapping each antibiotic to its class and subtype. |
| `time-to-event-augmented-queries/antibiotic_class_exposure.ipynb` | Prior exposure summarized by antibiotic class. |
| `time-to-event-augmented-queries/antibiotic_subtype_exposure.ipynb` | Prior exposure summarized by antibiotic subtype. |

---

## De-identification and Quality Control

| Script | Description |
|--------|-------------|
| `post_pulling/data_deduplication.py` | Removes duplicate rows from each table. |
| `post_pulling/filter_different_patients.py` | Removes culture orders whose accession number was issued to two different patients. |
| `post_pulling/filter_deid_conflicts.py` | Applies the same check to tables that already carry de-identified IDs. |
| `post_pulling/filter_multi_encounter_cultures.py` | Removes cultures linked to more than one encounter because their collection time falls inside overlapping encounter windows. |
| `post_pulling/deidentify_time.py` | Shifts every date by the patient's random offset. Shifted columns carry a `_shifted` suffix. `deidentify_time.ipynb` is the notebook version. |

---

## Requirements

- Python 3 with the packages in `requirements.txt`.
- R with `data.table` for the R versions of the laboratory and vital-sign scripts.

---

## Notes

- **Relationship to ARMD-Stanford.** ARMD-MGB follows the table structure and feature naming of the Stanford Antibiotic Resistance Microbiology Dataset (ARMD) so that analyses can be compared across the two sites. The scripts were adapted from the ARMD-Stanford SQL pipeline ([HealthRex CDSS, Aim 1a](https://github.com/HealthRex/CDSS/tree/master/scripts/antibiotic-susceptibility/Aim%201a)) and mostly keep the names of the corresponding queries.
- **Source data.** The scripts read MGB source extracts that are not distributed with this repository. Input and output locations are set in a configuration block at the top of each script or notebook.

---

## Citation

If you use ARMD-MGB, please cite:

> Wei, Z., & Kanjilal, S. (2025). Antibiotic Resistance Microbiology Dataset Mass General Brigham (ARMD-MGB). PhysioNet. RRID:SCR_007345. https://doi.org/10.13026/2r5k-b955

---

## References

- ARMD-MGB on PhysioNet: <https://physionet.org/content/armd-mgb/>
- Nateghi Haredasht, F. et al. Antibiotic Resistance Microbiology Dataset (ARMD): A Resource for Antimicrobial Resistance from EHRs. *Scientific Data* 12, 1299 (2025). <https://doi.org/10.1038/s41597-025-05649-7>
- [HealthRex CDSS – Aim 1a](https://github.com/HealthRex/CDSS/tree/master/scripts/antibiotic-susceptibility/Aim%201a)

## License

Released under the [MIT License](LICENSE).
