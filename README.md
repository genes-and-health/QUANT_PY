# `QUANT_PY` pipeline — Description

<img src="GNH%20TD%20Logo.png" alt="G&H Team Data logo" width=25%>

## Authors

* Stuart Rison
* Mike Samuels

Based on the "Quantitive Traits processing pipeline Jan 2025 redo" `pipeline_jan_2025.md` R-script written by Ben Jacobs with contributions from Saeed Bidi, Sarah Finer, Sam Hodgson, Stravoula Kanoni, Rohini Mathur, Caroline Morton, Daniel Stow, David van Heel, Julia Zollner.

## Summary

The Genes & Health (G&H) `QUANT_PY` pipeline extracts and processes quantitative data from G&H phenotype data.  It creates files and covariate files suitable for `regenie` \[G/Ex\]WAS analysis as well as generic files for each quantitative trait at a _per individual_ level (one row per individual summarising the individual's values for the trait) and a _per result_ level (one line per individual-result).  The pipeline processes HES data to identify admitted patient care (APC) episodes.  From this, three versions of the created files are generated: 1) All data, 2) Out of hospital data (without APC + a buffer), 3) In hospital data (within APC + a  buffer).

## Input data
### Trait files
The pipeline requires 3 trait input files: 1. `trait_features.csv`, 2. `trait_aliases_long.csv`, 3. `unit_conversions.csv`.

#### _`trait_features.csv`_
This file lists all the quantitative traits currently extracted from the phenotype data.  The .csv file has 4 columns: trait,target_units,min,max.

<details>
   
<summary>"trait_features.csv" file extract</summary>
  
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

</details>
So, for example, there is a "2h postprandial glucose" trait which is reported in "millimol/L", excluding any value less than 0.6 millimol/L or over 45.0 millimol/L.

#### _`trait_aliases_long.csv`_
SNOMED codes are missing for some of the G&H data pulls, this means that quantitative trait extraction is based on free-text trait descriptions (aka `original_term` within the script).  This file assigns all valid trait descriptions (aliases) to a trait.

<details>
   
<summary>"trait_aliases_long.csv" file extract</summary>
  
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
creatinine,Creatinine Serum
creatinine,Creatinine level (observable entity)
```

</details>

#### _`unit_conversions.csv`_
The same trait may be measured in different units depending of the setting (e.g. primary vs secondary care) or the data source (trust 1 vs trust 2).  This file allows unit conversions if a trait is in a valid but undesired unit which can be converted to a target_unit (as defined in `trait_features.csv`).  It also acts as a synonym dictionary to standardise unit terminology, for example, `nmol/L` is converted into the preferred term `nanomol/L`. Such conversions can be identified by a `multiplication_factor` of 1.0. 

<details>
   
<summary>"unit_conversions.csv" file extract</summary>
  
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
</details>

### Hospital admission data

`QUANT_PY` uses NHS England Digital Hospital Episode Statistics (HES) Admitted Patient Care (APC) data to identify periods of hospitalisation --as certain hospital day treatments are logged as APC events (e.g. immunotherapy infusions), **only APC episodes >2 calendar days are considered as hospitalisation**.  At present, `QUANT_PY` does not exclude data obtained during a A&E episode (HES AE + ECDS) unless this leads to a hospital admission (in which case it is "subsumed" by an APC episode).  However, the script is written such that it could accommodate these if needed/desired.

`QUANT_PY` **extends the hospitalisation episode by a 2-week buffer on either side of the APC event** on the basis that individuals admitted to hospital are typically unwell in the days leading to hospitalisation and may be discharged recovering, but prior to a return to their baseline status.

### Phenotype data
The pipeline imports G&H phenotype data in `/library-red/phenotypes_rawdata/`.  These data are from the following sources:
1. **DSA__BartHealth_NHS_Trust**: Secondary care data from the Barts Health NHS Trust \[North East London: ~40,000 individuals with data\]
2. **DSA__BradfordTeachingHospitals_NHSFoundation_Trust**: Secondary care data from the Bradford Teaching Hospitals NHS Trust \[Bradford and environs: ~1,700 individuals with data\]
3. **DSA__Discovery_7CCGs**: Primary care data from the North East London ICS \[North East London: ~45,000 individuals with data\]
4. **DSA_NHSDigital**: Data from the National Diabetes Audit (NDA) \[England-wide: ~13,000 individuals with data]

Phenotype files processed are listed in [Appendix A](#appendix-a-list-of-processed-phenotype-files).

#### Notes about the processing of quantitative data
1. `QUANT_PY` **does not use incremental data generation**. Every time a release is produced, all current and historically collected data are read in, concatenated and **then** deduplicated.
2. The aim of all phenotype processing steps is to produce a combined dataframe with the following columns:

| pseudo_nhs_number | test_date | original_term | result | result_value_units | provenance | source | hash |
|-------------------|-----------|---------------|--------|--------------------|------------|--------|------|
| 64-char pseudo NHS number | YYYY-MM-DD | free-text term e.g. "Mean Cell Volume" | result (flot64), e.g. 83.8 | unit of result, e.g. "Femtolitre" | file(s) data extracted from e.g. "2023_05_Barts_path" | source ("primary care" or "secondary care" | A hash value used to deduplicate data (unsigned int64) |

3. All `QUANT_PY` output files are derived from the above. 

## Pipeline steps
It is advisable to run the pipeline on a VM with lots of memory, typically an `n2d-highmem` 32 processor VM with 256Gb memory.

> [!TIP]
> All intermediary files are available in [`.arrow` format](https://arrow.apache.org/overview/)
>
> _(This link does not automatically open in a new window. Use CTRL+click (on Windows and Linux) or CMD+click (on MacOS) to open the link in a new window)_
> 

### STEP 0: Transfer phenotype data to `ivm`
Phenotype data is large in both size and number of files, and stored in different directories at different directory depth.  Buffering issues affect processing of data directly from the `/library-red/` Google Cloud bucket.  It is therefore simpler to copy all phenotype file to the `ivm` running `QUANT_PY`.  This transfer can be effected within the pipeline by setting a pipeline flag.

### STEP 1: Import phenotype files with appropriate pre-processing
`R` is very good at handling "raggedness" but in doing so, it makes assumptions.  This can lead to the "wrong" data ending in a column.  Python can also import .csv/.tsv/.tab files and make assumptions about the seprators/raggedness/column data type but in `QUANT_PY` this is intentionally and explicitly avoided.  This means that some files need to be pre-processed.  This take the form of one or more of the following pre-processing operations:

<details>

   <summary>Pre-processing operations</summary>
   
   1. **Commas in double-quotes stripping**: Exclude any row with double quoted text with one or more commas in in.  This excludes <2% of rows but means that the parsing behaviour is consistent and predictable.
   2. **Exclude "unterminated" double-quotes**: Some rows include a double-quote not paired with a second double-quote before the next separator.  In such cases, the importing functions often "glob" all text in subsequent rows until another double-quote is found.  Therefore, these rows are excluded.
   3. **Excluded rows with non-standard number of fields**: some rows may have additional/fewer separators either intentionally or erroneously creating additional/deleting fields.  `QUANT_PY` rejects any lines with a non-standard number of separators.
   4. **Strip double-quote**: This can be applied to non comma-delimited data files.  In some such files, double-quotes can appear singly ("), doubly ("") or even triply (""")
      
</details>

Processed files are listed in [Appendix A](#appendix-a-list-of-processed-phenotype-files).

#### Intermediary `.arrow` files
These can be useful for debugging purposes or for researchers interested in phenotypic data from a specific source.

<details>

   <summary>Per provenance `.arrow` files</summary>
   
   These can be found in the following directories (with an examplar `.arrow` file listed for each directory:
   * **`.../data/primary_care/arrow/`**: `2024_12_Discover_path.arrow`
   * **`.../data/secondary_care/arrow`**: `2023_05_Bradford_measurements.arrow`
   * **`.../data/secondary_care/nda`**: `2025_04_formatted_nda.arrow`

</details>

### STEP 2: Progressively merge files to create COMBO file

Files are merged in the following order (with deduplication after every merge operation).  Again, intermediate file directories are available as listed:
1. Primary care data + NDA data (primary care data): **`.../data/combined_datasets/`**: `YYYY_MM_Combined_primary_care.arrow`
2. Barts data + Bradford data (secondary care data): **`.../data/combined_datasets/`**: `YYYY_MM_Combined_secondary_care.arrow`

Finally, primary and secondary care data are merged.

> [!TIP]
> The final output of the multiple merges is considered a key output of the pipeline and is therefore stored in the **`.../outputs/`** directory:
> 
>  **`../outputs/`**: `YYYY_MM_Combined_all_sources.arrow`
>
> The dataframe is referred to as the **COMBO**.

On 2025-04-01, the **COMBO** `2025_04_Combined_all_sources.arrow` was `78,582,474` rows long.

### STEP 3: Add hospitalisation status column

#### Import HES data

Admitted Patient Care (APC) episodes are extracted from HES data pulls of 2021-09, 2023-07, 2024-10, and 2025-03. HES APC data are imported, cleaned up and deduplicated.

> [!TIP]
> The output of the HES APC pulls merges and deduplication are in the **`.../data/combined_datasets/`** directory:
> 
>  **`.../data/combined_datasets/`**: `YYYY_MM_Combined_HES.arrow`

#### Flagging COMBO result dates

The HES dataset is used to define three APC `region_types`:
* **APC**: a date span for the APC
* **buffer_before**: a date span of 14d prior to addmission date
* **buffer_after**: a date span of 14d after discharge date

**COMBO** test results are flagged to none (`null`) if they fall out of the above listed three region types, or to one or more of the region types, by joining the APC data to **COMBO**.  For example, a date may exist within an APC period (flagged as `["APC"]`), or within an APC and a buffer_before (for example if the date falls both within an APC and within the buffer_before of a subsequent APC; flagged as `["APC", "buffer_before"]`).

By extension, test result dates can be classifed in one of 11 (some non-mutually exclusive) categories.

<details>
   <summary>Test result data categories</summary>
   
   1. **`IN_APC_ONLY`**: test results collected in an actual hospitalisation episode not overlapping with the buffer of another APC
   2. **`IN_APC_ANY`**: test results collected in an actual hospitalisation episode (which may overlap another APC's 14d buffer)
   3. **`IN_BUFFER_BEFORE_ONLY`**: test results collected within the 14d prior to a hospitalisation episode (not overlapping with another APC or another APC's buffer)
   4. **`IN_BUFFER_BEFORE_ANY`**: test results collected within the 14d prior to a hospitalisation episode (may overlap with an APC)
   5. **`IN_BUFFER_AFTER_ONLY`**: test results collected within the 14d following a hospitalisation episode (not overlapping with another APC or another APC's buffer)
   6. **`IN_BUFFER_AFTER_ANY`**: test results collected within the 14d following a hospitalisation episode (may overlap with an APC)
   7. **`IN_BUFFERS_ONLY`**: test results collected _either_ within the 14d prior to, or following, a hospitalisation episode _and_ not during the actual hospitalisation period
   8. **`IN_BUFFERS_ANY`**: test results collected _either_ within the 14d prior to, or following, a hospitalisation episode (may overlap with one or more APCs)
   9. **`IN_TOTAL_EXCLUSION_ZONE`**: test results collected within a period from 14 days prior to a hospitalisation episode to 14 days following a hospitalisation episode _including_ the hospitalisation period
   10. **`OUT_OF_APC`**: test results collected outside of any APC hospitalisation episode
   11.  **`OUT_OF_TOTAL_EXCLUSION_ZONE`**: test results collected outside of any buffered hospitalisation episode (hospitalisation episode + 14 days either side)   

</details>

In practice, only three possible date statuses are considered: "all" (all quantitative results), "out_hospital" (i.e. `OUT_OF_TOTAL_EXCLUSION_ZONE`) and "in_hospital" (i.e. `IN_TOTAL_EXCLUSION_ZONE`) 

### STEP 4: Perform unit conversions and flag out-of-range **COMBO** results

COMBO is joined to a denormalised traits dataframe (`traits_features` x `trait_aliases`) which identifies COMBO row with traits to extract, their target units and their valid range.  **Unit conversions are performed where possible and applicable** and final results are flagged as:

* **`below_min`**: result lower than the minimum value set for this trait.  These will subsequently be excluded.
* **`ok`**: result within valid range for this trait.
* **`above_max`**: result higher than the maximum value set for this trait.  These will subsequently be excluded.

> [!NOTE]
> HbA1c values in `%` (percentages) are converted to values in `millimol/mol` using the following equation: $mmol\/mol \[IFCC\] =  (10.93*percentage \[NGSP/UKPD\]) - 23.50$  
> See [National Glycohemoglobin Standardization Program](https://ngsp.org/ifccngsp.asp) 
> _(This link does not automatically open in a new window. Use CTRL+click (on Windows and Linux) or CMD+click (on MacOS) to open the link in a new window)_ 
> 
> Because percentages apply to traits other than HbA1c, this conversion cannot be performed using the `unit_conversions.csv` and needs to be "hard-coded" in the pipeline.  

### Step 5: COMBO restricted to valid pseudoNHS numbers and valid demographics

When volunteers take part in stage 1 of Genes & Health, their questionnaire and consent form is labelled with the ID number on the Oragene saliva tube (style: `15001502031604`). These Oragene IDs are then used to label genetic samples (e.g. GSA chip or exome seq). They also label the Questionnaire (aka `S1QT`). Some people have taken part twice (or more than twice) over the years in Genes & Health, and will have a different Oragene ID each time.  The **OrageneID** is the link to genetic data, the **pseudoNHSnumber** is the link to phenotypic data.

Step 5 uses a `YYYY_MM_DD_MegaLinkage_forTRE.csv`&trade; source file to allow these linkages.

<details>
   <summary>MegaLinkage&trade; file columns</summary>
   
   * **OrageneID**: 14 digit unique OrageneID
   * **Number of OrageneIDs with this NHS number (i.e. taken part twice or more)**: 1, 2, 3 or 4.  Typically 1 (single participation), 2 (~10% individuals), 3 (~0.7% individuals), 4 (~0.05% individuals)
   * **S1QST gender**: 1=male, 2=female
   * **HasValidNHS**: "yes", "no"
   * **pseudonhs_2024-07-10**: pseudoNHSnumber
   * **51176GSA-T0PMedr3 Jan2024release**: GSA identifier (OrageneID_GSAID_RunID; OrageneID as above + '\_' + GSAID = 12digit identified + '\_' + RunID = row/col in format RxxCxx where x = digit)
   * **44628exomes_release_2023-JUL-07**: Broad Institute Exome sequencing ID (= GNH-+OrageneID)
   * **55273exomes_release_2024-OCT-08**: Broad Institute Exome sequencing ID (= GNH-+OrageneID)
   
</details>

<details>
   <summary>Extract from <code>2025_02_10_MegaLinkage_forTRE.csv</code>&trade; [redacted]</summary>

   ```
   OrageneID,Number of OrageneIDs with this NHS number (i.e. taken part twice or more),S1QST gender,HasValidNHS,pseudonhs_2024-07-10,51176GSA-T0PMedr3 Jan2024release,44628exomes_release_2023-JUL-07,55273exomes_release_2024-OCT-08
   1xxxxxxxxxx2,1,1,yes,2xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx0,1xxxxxxxxxx2_2xxxxxxxxxx2_Rxxxx2,GH-1xxxxxxxxxx2,GMH-1xxxxxxxxxx2
   1xxxxxxxxxx0,1,2,yes,8xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx0,1xxxxxxxxxx0_2xxxxxxxxxx2_Rxxxx2,GH-1xxxxxxxxxx2,GMH-1xxxxxxxxxx0
   1xxxxxxxxxx2,1,1,yes,9xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx8,1xxxxxxxxxx2_2xxxxxxxxxx2_Rxxxx2,GH-1xxxxxxxxxx2,GMH-1xxxxxxxxxx2
   1xxxxxxxxxx0,1,2,yes,8xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx5,1xxxxxxxxxx0_2xxxxxxxxxx2_Rxxxx2,GH-1xxxxxxxxxx2,GMH-1xxxxxxxxxx0
   1xxxxxxxxxx0,1,1,yes,0xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx8,1xxxxxxxxxx0_2xxxxxxxxxx2_Rxxxx2,GH-1xxxxxxxxxx2,GMH-1xxxxxxxxxx0
   1xxxxxxxxxx7,1,1,yes,8xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx0,1xxxxxxxxxx7_2xxxxxxxxxx2_Rxxxx1,GH-1xxxxxxxxxx7,GMH-1xxxxxxxxxx7
   1xxxxxxxxxx8,1,2,yes,4xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx4,1xxxxxxxxxx8_2xxxxxxxxxx2_Rxxxx2,GH-1xxxxxxxxxx8,GMH-1xxxxxxxxxx8
   ```

</details>

#### Filter to valid pseudoNHS number
There are approximately 1,000 rows excluded by this pseudoNHS validation.
Possible reasons:

1. subject asked to be removed/withdrawn
2. subject died
3. subject had multiple pseudoNHSnumbers which have been merged

#### Filter to valid demographics
Use `MegaLinkage file`$trade; to link to Stage 1 questionnaire data found in `/library-red/genesandhealth/phenotype_raw_data/QMUL__Stage1Questionnaire`.  This gives patient MONTH/YEAR of birth.  All volunteers are assumed to be born on the first day of a month.  Results are excluded if age at test is less than or equal to zero (i.e. result before birth of volunteer), or the result is dated beyond the date of the run execution (i.e. results dated to the future).

Results obtained prior to 16 years of age are also excluded.

#### Filter to within-range values
Only rows with a `range_position` (as defined in [STEP 4](#step-4-perform-unit-conversions-and-flag-out-of-range-combo-results)) equal to `ok` (cf. `below_min` and `above_max`) are kept.

### STEP 6: Window data in 10-day windows
Because the same quantitative result can come from multiple sources with a similar but non-identical date (e.g. a secondary care result is registered in an individual's primary care record with the date it was received in primary care rather than the actual test result date), the quantitative data are "windowed".

Non-overlapping 10-day windows are applied and any identical test results within these windows are deduplicated even if the result dates differ.  The earliest instance of the result is kept.

### STEP 7: Generate output files
These can all be found in the **`.../outputs/`** directory

## OUTPUT FILES
The following files are generated from the QCed **COMBO** generated in **STEP 6**

1. **per trait files** \[`../outputs/individual_trait_files/`; subdirectories: `in_hospital`, `out_hospital`, `all`\]:
     - **`_{trait}_readings_at_unique_timepoints.csv`**: one validated result per row (columns: `pseudo_nhs_number, trait, unit, value, date, gender, age_at_test, minmax_outlier`) 
     - **`_{trait}_per_individual_stats.csv`**: one row per volunteer (`pseudo_nhs_number, trait, median, mean, max, min, earliest, latest, number_observations`)
2. **per trait plots** \[`../outputs/individual_trait_plots/`; subdirectories: `in_hospital`, `out_hospital`, `all`\]:
      - **`_{trait}_{setting}.svg`**: Histograms of trait log10(values) for trait separated M and F listing median, mean, min, max, number individuals, number observations
3. **regenie files** \[`../outputs/regenie/`; subdirectories: `in_hospital`, `out_hospital`, `all` and `covariate_files`\]:
      - **`_{trait}_{setting}_[regenie_51|regenie_55].tsv`**: regenie files for 51kGWAS and 55kExome analyses
      - **`./covariate_files/_{setting}_[regenie_51|regenie_55]_megawide.tsv`**: regenie covariate files allowing age at test analyses (cf. age on joining Genes and Health)
4. **reference COMBO files** \[`../outputs/reference_combo_files/`\]:
      - **`_Combined_all_sources.arrow`**: the "raw" merger of primary, secondary and NDA data.  No QC, no restriction to the 111 traits extracted in `version010_2025_04`
      - **`_Combined_traits_NHS_and_demographics_restircted_pre_10d_windowing`**: above file processed to limit to valid NHS number, valid demographics and valid values but _not_ windowed (end of **STEP 5**)
      - **`_Combined_traits_NHS_and_demographics_restircted_post_10d_windowing`**: above file processed to limit to valid NHS number, valid demographics and valid values _and_ windowed (end of **STEP 6**)


# Appendix A: List of processed phenotype files
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
