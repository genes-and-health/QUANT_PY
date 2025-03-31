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
Here is an extract:
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
So, for example, there is a "2h postprandial glucose" trait which is reported in "millimol/L" excluding any value <0.6 millimol/L or over 45.0 millimol/L.

## Phenotype data
The pipeline imports G&H phenotype data from `/library-red/phenotypes_rawdata/`, that is, from the following sources:
1. DSA__BartHealth_NHS_Trust: Secondat care data from the Barts Health NHS Trust
