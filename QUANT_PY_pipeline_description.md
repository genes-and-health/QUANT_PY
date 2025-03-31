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

## Phenotype data
The pipeline imports G&H phenotype data from `/library-red/phenotypes_rawdata/`, that is, from the following sources:
1. DSA__BartHealth_NHS_Trust: Secondat care data from the Barts Health NHS Trust
