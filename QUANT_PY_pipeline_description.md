# `QUANT_PY` pipeline — Description

<img src="GNH%20TD%20Logo.png" alt="G&H Team Data logo" width=25%>

## Authors

* Stuart Rison
* Mike Samuels

Based on the "Quantitive Traits processing pipeline Jan 2025 redo" `pipeline_jan_2025.md` R-script written by Ben Jacobs with contributions from Saeed Bidi, Sarah Finer, Sam Hodgson, Stravoula Kanoni, Rohini Mathur, Caroline Morton, Daniel Stow, David van Heel, Julia Zollner.

## Summary

The Genes & Health (G&H) `QUANT_PY` pipeline extracts and processes quantitative data from G&H phenotype data.  It creates files and covariate files suitable for `regenie` \[G/Ex\]WAS analysis as well as generic files for each quatitative trait at a _per individual_ level (one row per individual summarising the individual's values for the trait) and a _per result_ level (one line per individual-result).  The pipeline processes HES data to identify admitted patient care (APC) episodes.  From this, three versions of the created files are generated: 1) All data, 2) Out of hospital data (without APC + a buffer), 3) In hospital data (within APC + a  buffer).

## Input data
### Trait files
The pipeline requires 3 trait input files: 1. `trait_features.csv`, 2. `trait_aliases_long.csv`, 3. `unit_conversions.csv`.

#### _`trait_features.csv`_
This file lists all the quantitative traits currently extracted from the phenotype data.  The .csv file has 4 columns: trait,target_units,min,max.
Extract:
```
trait,target_units,min,max
2h postprandial glucose,millimol/L,0.6,45.0
AFP,kU/L,1e-100,10000.0
ALP,units/L,8.0,1500.0
ALT,units/L,5.0,1500.0
APTT,seconds,1e-100,100.0
AST,units/L,3.0,1000.0
Albumin,g/L,10.0,80.0
Alcohol units per week,units/week,0.0,350.0
```
So, for example, there is a "2h postprandial glucose" trait which is reported in "millimol/L", excluding any value less than 0.6 millimol/L or over 45.0 millimol/L.

#### _`trait_aliases_long.csv`_
SNOMED codes are missing for some of the G&H data pulls, this means that quantitative trait extraction is based on free-text trait descriptions (aka `original_term` within the script).  This file assigns all valid trait descriptions (aliases) to a trait.
Extract:
```
trait,alias
2h_postprandial_glucose,2h postprandial glucose
2h_postprandial_glucose,"Glucose tolerance test, 2 hour post prandial (procedure)"
...
Alcohol_units_per_week,Alcohol units per week
Alcohol_units_per_week,Alcohol units/week (qualifier value)
...
Blood_ketones,Blood ketone level (observable entity)
Blood_ketones,POCT Blood Ketones
...
Platelets,Platelet count
Platelets,Platelet count (observable entity)
...
creatinine,Creatinine Serum
creatinine,Creatinine level (observable entity)
```
#### _`trait_aliases_long.csv`_
The same trait may be measured in different units depending of the setting (e.g. primary vs secondary care) or the data source (trust 1 vs trust 2).  This file allows unit concersions to see if a trait in a valid but undesired unit can be converted to a target_unit (as defined in `trait_features.csv`).
Extract:
```
result_value_units,target,multiplication_factor
%,%,1.0
*10^9/l,10^9/L,1.0
g/L,mg/L,1000.0
IU/L,units/L,1.0
Kg,grams,1000.0
mg/L,g/L,0.001
miu/L,milliunits/L,1.0
nmol/L,nanomol/L,1.0
Units/Day,units/week,7.0
```
### Hospital admission data

`QUANT_PY` uses NHS England Digital Hospital Episode Statistics (HES) APC data to identify periods of hospitalisation --as certain hospital day treatments are logged as APC events (for example immunotherapy infusions), **only APC >2 calendar days are considered as hospitalisation**.  At present, `QUANT_PY` does not exclude data obtained during a A&E episode (HES AE + ECDS) unless this leads to a hospital admission (in which case it is "subsumed" by an APC episode) however, the script is written such that it could be accommodate these if needed/desired.

`QUANT_PY` **extends the hospitalisation episode by a 2-week buffer on either side of the APC event** on the basis that individiual admitted to hospital are typically unwell in the day leading to the hospitalisation and may be discharged recovering but prior to a return to their baseline status.

### Phenotype data
The pipeline imports G&H phenotype data from `/library-red/phenotypes_rawdata/`, that is, from the following sources:
1. **DSA__BartHealth_NHS_Trust**: Secondary care data from the Barts Health NHS Trust \[North East London: ~40,000 individuals with data\]
2. **DSA__BradfordTeachingHospitals_NHSFoundation_Trust**: Secondary care data from the Bradford Teaching Hospitals NHS Trust \[Bradford and environs: ~1,700 individuals with data\]
3. **DSA__Discovery_7CCGs**: Primary care data from the North East London ICS \[North East London: ~45,000 individuals with data\]
4. **DSA_NHSDigital**: Data from the National Diabetes Audit (NDA) \[England-wide: ~13,000 individuals with data]

All files phenotypes processed are listed in [Appendix A]

#### Notes about the processing of qunatitative data
1. `QUANT_PY` **does not use incremental data generation**. Everytime a release is produced, all current and historically collected data a read in, concatenated and **then** deduplicated.
2. The aim of all phenotype processing steps is to produce a combined dataframe with the following columns:

| pseudo_nhs_number | test_date | original_term | result | result_value_units | provenance | source | hash |
|-------------------|-----------|---------------|--------|--------------------|------------|--------|------|
| 64-char pseudo NHS number | YYYY-MM-DD | free-text term e.g. "Mean Cell Volume" | result (flot64), e.g. 83.8 | unit of result, e.g. "Femtolitre" | file(s) data extracted from e.g. "2023_05_Barts_path" | source ("primary care" or "secondary care" | A hash value used to deduplicte data (unsigned int64) |

3. All `QUANT_PY` output files are derived from the above. 

## Pipeline steps
It is advisable to run the pipeline on a VM with lots of memory, typically an `n2d-highmem` 32 processor VM with 256Gb memory.
### STEP 0: Transfer phenotype data to `ivm`
Phenotype data is large in both size and number of files and stored in different direcotries at different directory depth.  Buffering issues affect processing of data directly from the `/library-red/` Google Cloud bucket.  It is therefore simpler to copy all phenotype file to the `ivm` running `QUANT_PY`.  This transfer can be effected within the pipeline by setting a pipeline flag.

### STEP 1: Import phenotype files with appropriate pre-processing
`R` is very good at handling "raggedness" but in doing so, it makes assumptions.  This can lead to the "wrong" data ending in a column.  Python can also import .csv/.tsv/.tab files and make assumptions about the seprators/raggedness/column data type but in `QUANT_PY` this is intentionally and explicitily avoided.  This means that some files need to be pre-processed.  This take the form of one of the following pre-processing operations:
1. **Commas in double-quotes stripping**: Exclude any row with double quoted text with one or more commas in in.  This excludes <2% of rows but mean that the parsing behaviour is consistent and predictable.
2. **Exclude "unterminated" double-quotes**: Some rows include a double-quote not paired with a second double-quote before the next separator.  In such cases, the importing functions often "glob" all text in subsequent rows until another double-quote is found.  Therefore these rows are excluded.
3. **Excluded rows with non-standard number of fields**: some rows may have additional/fewer separators either intentionally or erroneously creating additional/deleting fields.  `QUANT_PY` rejects any lines with a non-standard number of separators.
4. **Strip double-quote**: This can be applied to non comma-delimited data files.  In some such files, double-quotes can appear singly ("), doubly ("") or even triply (""")
Processed files are listed in [Appendix A]

Per provenance `.arrow` files are prepared at each step.  These can be found in the following directories (with an examplar `.arrow` file listed for each directory:
* **`.../data/primary_care/arrow/`**: `2024_12_Discover_path.arrow`
* **`.../data/secondary_care/arrow`**: `2023_05_Bradford_measurements.arrow`
* **`.../data/secondary_care/nda`**: `2025_04_formatted_nda.arrow`

### STEP 2: Progressively merge files

Files are merged in the following order (with de-duplication after every merge operation).  Again, intermediate file directories and files are listed:
1. Primary care data + NDA data (primary care data): **`.../data/combined_datasets/`**: `2025_04_Combined_primary_care.arrow`
2. Barts data + Bradford data (secondary care data): **`.../data/combined_datasets/`**: `2025_04_Combined_secondary_care.arrow`

Finally, primary and secondary care data are merged.  The final output of the multiple merges is considered a key output of the pipeline is available from:
**`../outputs/`**: `YYYY_MM_Combined_all_sources.arrow`



### 

# Appendix A  -- List of processed phenotype files
```
# Primary care
.../DSA_Discovery_7CCGs/2022_04_Discovery/GNH_thw/nech-phase2-outfiles_merge/GNH_thw/nech_observations_output_dataset_20220423.csv
.../DSA_Discovery_7CCGs/2022_04_Discovery/GNH_bhr-phase2-outfiles_merge/GNH_bhr_observations_output_dataset_20220412.csv
.../DSA_Discovery_7CCGs/2022_12_Discovery/GNH_thw/nech-phase2-outfiles_merge/gh2_observations_output_dataset_20221207.csv
.../DSA_Discovery_7CCGs/2022_12_Discovery/GNH_bhr-phase2-outfiles_merge/gh2_observations_dataset_20221207.csv
.../DSA_Discovery_7CCGs/2023_04_Discovery/gh3_observations.csv
.../DSA_Discovery_7CCGs/2023_11_Discovery/gh3_observations.csv
.../DSA_Discovery_7CCGs/2024_04_Discovery/gh3_observations.csv
.../DSA_Discovery_7CCGs/2024_12_Discovery/gh3_observations.csv

# NDA
.../DSA_NHSDigitalNHSEngland/2024_10/NDA/NIC338864_NDA_BMI.txt
.../DSA_NHSDigitalNHSEngland/2024_10/NDA/NIC338864_NDA_CHOL.txt
.../DSA_NHSDigitalNHSEngland/2024_10/NDA/NIC338864_NDA_HBA1C.txt
.../DSA_NHSDigitalNHSEngland/2024_10/NDA/NIC338864_NDA_BP.txt

# Secondary care -- Bradford
.../DSA_BradfordTeachingHospitals_NHSFoundation_Trust/2023_05_BTHFT/1578_gh_lab_results_2023-06-09_noCR.ascii.redacted.tab
.../DSA_BradfordTeachingHospitals_NHSFoundation_Trust/2024_12_BTHFT/1578_gh_lab_results_2024-12-05.ascii.redacted.tab
.../DSA_BradfordTeachingHospitals_NHSFoundation_Trust/2022_06_BTHFT/1578_gh_cerner_measurements_2022-06-10_redacted.tsv
.../DSA_BradfordTeachingHospitals_NHSFoundation_Trust/2024_12_BTHFT/1578_gh_cerner_measurements_2024-12-05.ascii.redacted.tab

# Secondary care -- Barts
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/Candida_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/SHBG_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/Fasting_Glucose.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/WhiteBloodCellCount_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/Monocytes_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/Basophil_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/HDL_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/DHEA_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/Random_Glucose.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/TotalCholesterol_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/Insulin_Antibodies.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/BileAcidSerum_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/ApolipoproteinB100_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/TriglycerideSerum_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/Testosterone_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/RBC_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/Lymphocytes_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/Eosinophil_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/MCV_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/LDL_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/ApolipoproteinA1_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/Neutrophils_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/Haematocrit_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/MCHC_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/Platelet_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/Oestradiol_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/GAD_Antibodies.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/HbA1c_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2021_04_PathologyLab/Prolactin_April2021.csv
.../DSA_BartsHealth_NHS_Trust/2022_03_ResearchDatasetv1.5/GH_Pathology_20220319T143_redacted_noHistopathologyReport.csv
.../DSA_BartsHealth_NHS_Trust/2023_05_ResearchDatasetv1.5/GH_Pathology_20230517T061.ascii_redacted.nohisto.tab
.../DSA_BartsHealth_NHS_Trust/2023_12_ResearchDatasetv1.6/GH_Pathology_20231218.ascii.nohisto.redacted2.tab
.../DSA_BartsHealth_NHS_Trust/2024_09_ResearchDataset/RDE_Pathology_ascii.nohisto.redacted2.csv
.../DSA_BartsHealth_NHS_Trust/2023_05_ResearchDatasetv1.5/GandH_Measurements_20230512T304.ascii.redacted.tab
.../DSA_BartsHealth_NHS_Trust/2023_12_ResearchDatasetv1.6/GandH_Measurements_20240423.ascii.redacted2.tab
.../DSA_BartsHealth_NHS_Trust/2024_09_ResearchDataset/RDE_Measurements.ascii.redacted2.tab
```
