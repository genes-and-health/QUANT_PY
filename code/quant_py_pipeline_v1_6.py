#!/usr/bin/env python
# coding: utf-8

# ## `/quant_py_pipeline`
# 
# ## Contributors
# * Stuart Rison
# * Mike Samuels
# 
# * Ben Jacobs and colleagues from G&H Team Data 
# 
# ## Motivation & brief overview
# > "I want these done yesterday" - D. van Heel, 2023.  
# > "Could I give you some QC tasks? Might come bit by bit. Lets do them stepwise" - D. van Heel, 2024.
# 
# > "I want these the day after the new HES data come in." - D. van Heel, 2025.      
# > "How about within 1 week?" - S. Rison, 2025.
# 
# This pipeline generates clean quantitative trait data from electronic healthcare records linked to the participants in Genes & Health. In brief, we do the following:
# 
# * Perform dataset-level collection and quality control
# * Combine the datasets
# * Perform trait-specific collection and quality control
# 
# ## Note
# Branched by SR/MS on 2025-03-27 from `quant_py_pipeline_admission_exclusion_thresholder_v0_1.ipynb`.
# 
# ---
# 
# ### Notable Versions
# - v1.00 (2025-03-27):
#     - Created the script
# - v1.10 (2025-03-28):
#     - Run-through to `combo` file completed (whooo!)    
# - v1.20 (2025-03-31):
#     - Inclusion of MS script to define and include/exclude/ignore hospitalisations periods from HES data
# - v1.30 (2025-04-02):
#     - Latest HES data
#     - Tweaks to hospitalisation handling
# - v1.40 (2025-04-03):
#     - Last version pre- HES region categorisation
# - v1.50 (2025-04-03):
#     - Last pre-production release
# - v1.60 (2025-04-04):
#     - This is the script used to create QUANT_PY `version010_2025_04`
# 

# # Code

# ## Set the run details

# In[ ]:


version = "version010"
mon = "04"
yr = "2025"


# In[ ]:


primary_keys = [
    "2022_04_Discovery_path",
    "2022_12_Discovery_path",
    "2023_03_Discovery_path",
    "2023_11_Discovery_path",
    "2024_07_Discovery_path",
    "2024_12_Discovery_path",
]


# In[ ]:


nda_keys = [
    "2024_10_NHSD_NHSE_NDA_path",
    #"2025_XX_NHSD_NHSE_NDA_path", # place holder for 2025 NHSED pull; but not requested so unlikely to appear
]


# In[ ]:


barts_keys = [
    "2021_04_Barts_path",
    "2022_03_Barts_path",
    "2023_05_Barts_path",
    "2023_05_Barts_measurements",
    "2023_12_Barts_path",
    "2023_12_Barts_measurements",
    "2024_09_Barts_path",
    "2024_09_Barts_measurements",
]


# In[ ]:


bradford_keys = [
    "2022_06_Bradford_measurements",
    "2023_05_Bradford_path",
    "2024_12_Bradford_measurements",    
    "2024_12_Bradford_path",
]


# In[ ]:


GNH_PALETTE = {
    "COBALT_BLUE": "#32449b",
    "EMERALD_GREEN": "#45c086",
    "MAGENTA": "#c44887",
    "PEACH_ORANGE": "#ff8070",
    "INDIGO_PURPLE": "#312849",
}


# In[ ]:


ALL_PROVENANCE_OPTIONS = primary_keys + nda_keys + barts_keys + bradford_keys


# In[ ]:


ALL_SOURCE_OPTIONS = ["primary_care", "secondary_care"]


# In[ ]:


# The script is run in ivm on ivm storage; in this process, most source data are copied from
# /library-red/ to /home/ivm/. The intermediary files are also initially stored on /home/ivm
# during the running of the pipeline.
# However, at completion of the pipeline, all useful intermediaries and final outputs are moved to
# /library-red/phenotypes_curated/version###_YYYY_MM/

PIPELINE_NAME = "QUANT_PY"
VERSION_FOLDER_NAME = f"{version}_{yr}_{mon}"

ROOT_FOLDER_LOCATION = "/home/ivm"
RED_FOLDER_LOCATION = f"/genesandhealth/red/{PIPELINE_NAME}/{VERSION_FOLDER_NAME}"

PIPELINE_VERSION_LOCATION = f"{ROOT_FOLDER_LOCATION}/{PIPELINE_NAME}/{VERSION_FOLDER_NAME}"
PIPELINE_DATA_LOCATION = f"{PIPELINE_VERSION_LOCATION}/data"
PIPELINE_INPUTS_LOCATION = f"{PIPELINE_VERSION_LOCATION}/inputs"
PIPELINE_HELPERS_LOCATION = f"{PIPELINE_VERSION_LOCATION}/helpers"
PIPELINE_LOGS_LOCATION = f"{PIPELINE_VERSION_LOCATION}/logs"

LIBRARY_RED_DATA_LOCATION = "/genesandhealth/library-red/genesandhealth/phenotypes_rawdata" 
NHSE_SUBLICENSE_DATA_LOCATION = "/genesandhealth/nhsdigital-sublicence-red"

PIPELINE_OUTPUTS_LOCATION = f"{PIPELINE_VERSION_LOCATION}/outputs"

PIPELINE_INDIVIDUAL_TRAIT_FILES_LOCATION = f"{PIPELINE_OUTPUTS_LOCATION}/individual_trait_files"
PIPELINE_INDIVIDUAL_TRAIT_PLOTS_LOCATION = f"{PIPELINE_OUTPUTS_LOCATION}/individual_trait_plots"


# In[ ]:


# The running of the pipeline takes place in a ivm.  The input files are stored in the ivm for this
# Once the run is officialised and released, the input files will be copied to /red/QUANT_PY.
# This means that for subsequent re-running(s) of the pipeline, the input files will have to be copied
# back from /red/QUANT_PY back into PIPELINE_INPUTS_LOCATION.
TRAIT_FEATURES_LOCATION = [
    PIPELINE_INPUTS_LOCATION,
    "trait_features.csv"
]


# In[ ]:


TRAIT_ALIASES_LONG_LOCATION = [
    PIPELINE_INPUTS_LOCATION,
    "trait_aliases_long.csv"
]


# In[ ]:


UNIT_CONVERSIONS_LOCATION = [
    PIPELINE_INPUTS_LOCATION,
    "unit_conversions.csv"
]


# In[ ]:


MEGA_LINKAGE_LOCATION = [
    "/",
    "genesandhealth",
    "library-red",
    "genesandhealth",
    "2025_02_10__MegaLinkage_forTRE.csv"
]


# In[ ]:


COMBO_POST_10D_WINDOWING_LOCATION = (
    PIPELINE_OUTPUTS_LOCATION,
    "reference_combo_files",
)


# In[ ]:


# ADMISSIONS_DATA_LOCATION = [
#     "/",
#     "genesandhealth",
#     "red",
#     "QUANT_PY",
#     "other",
#     "hes_data.parquet"
# ]


# ## Import relevant packages

# In[ ]:


import polars as pl
from cloudpathlib import AnyPath
import datetime
import dateutil
import glob
from enum import Enum
from collections import defaultdict
import subprocess
from itertools import chain, combinations
import gc


# In[ ]:


# Separate imports for file copy; can be commented out when not performing copy. 
# Once copied it stays copied, until manually moved or deleted.

from tqdm import tqdm
import shutil


# In[ ]:


# Plotting packages
import altair as alt
alt.data_transformers.enable("vegafusion")
alt.renderers.enable("svg")


# In[ ]:


# Polars namespace additions
# In subsequent version this code may be integrated with the establisted TRE Tools package

@pl.api.register_lazyframe_namespace("TRE")
class TRETools:
    def __init__(self, lzdf: pl.LazyFrame) -> None:
        self._lzdf = lzdf
        
    def unique_with_logging(self, *args, label: str = "Unique", **kwargs) -> pl.LazyFrame:
        before = self._lzdf.collect().height
        filtered_lzdf = self._lzdf.unique(*args, **kwargs)
        after = filtered_lzdf.collect().height
        
        if before > 0:
            change = ((after - before) / before) * 100
            change_str = f" ({'+' if change > 0 else ''}{change:.1f}%)"
        
        unchanged = " (row count unchanged)" if after == before else ""
        
        print(f"[{label}: on {args}] Before unique: {before} rows, After unique: {after} rows{unchanged}{change_str}")
        return filtered_lzdf    
    
    def filter_with_logging(self, *args, label: str = "Filter", **kwargs) -> pl.LazyFrame:
        before = self._lzdf.collect().height
        filtered_lzdf = self._lzdf.filter(*args, **kwargs)
        after = filtered_lzdf.collect().height
        
        if before > 0:
            change = ((after - before) / before) * 100
            change_str = f" ({'+' if change > 0 else ''}{change:.1f}%)"
        
        unchanged = " (row count unchanged)" if after == before else ""
        print(f"[{label}] Before filter: {before} rows, After filter: {after} rows{unchanged}{change_str}")
        return filtered_lzdf
    
    def join_with_logging(
        self,
        other: pl.LazyFrame,
        *args,
        how: str = "inner",
        label: str = "Join",
        **kwargs
    ) -> pl.LazyFrame:
        left_before = self._lzdf.collect().height
        right_before = other.collect().height
        joined_lzdf = self._lzdf.join(other, *args, how=how, **kwargs)
        after = joined_lzdf.collect().height
        
        if left_before > 0:
            change = ((after - left_before) / left_before) * 100
            change_str = f" ({'+' if change > 0 else ''}{change:.1f}%)"
        
        unchanged = " (row count unchanged)" if after == left_before else ""
        print(f"[{label}] Join type: {how.upper()}")
        print(f"[{label}] Left: {left_before} rows, Right: {right_before} rows -> After: {after} rows{unchanged}{change_str}")
        return joined_lzdf


# ## Utility Functions
# 
# These may be deleted or formalised in subsequent versions of the Pipeline. 
# 
# Many were written as exploratory functions during pipeline design.

# ### HES functions

# In[ ]:


BUFFER_BEFORE_DAYS = 14
BUFFER_AFTER_DAYS = 14


# In[ ]:


def add_buffers(
    lf: pl.LazyFrame, 
    buffer_before_duration_in_days: pl.UInt16 = BUFFER_BEFORE_DAYS,
    buffer_after_duration_in_days: pl.UInt16 = BUFFER_AFTER_DAYS,
    padding_before_duration_in_days: pl.UInt16 = 0,
    padding_after_duration_in_days: pl.UInt16 = 0,
    id_column: pl.Utf8 = "id",
    hospital_start_date_column: pl.Utf8 = "start_date", 
    hospital_end_date_column: pl.Utf8 = "end_date"
) -> pl.LazyFrame:
    """
    Adds Buffer Before|After, Padding Before|After to the hospital stay periods.
    
    :param lf: input LazyFrame containing pseudo_nhs_number, start_date, end_date representing hospital admission/discharge dates
    :param buffer_before_duration_in_days: Days before hospital stay considered as buffer
    :param buffer_after_duration_in_days: Days after hospital stay considered as buffer
    :param padding_before_duration_in_days: Padding between hospital stay and buffer before
    :param padding_after_duration_in_days: Padding between hospital stay and buffer after
    :param hospital_start_date_column: Column name for hospital start date
    :param hospital_end_date_column: Column name for hospital end date
    :return: LazyFrame with additional region types
    """

#     lf = lf.rename({id_column: "pseudo_nhs_number"}, strict=False)

    lf = lf.rename({
            id_column: "id",
            hospital_start_date_column: "start_date",
            hospital_end_date_column: "end_date",
        })
    
    lf_extended = (
        lf.with_columns([
            (pl.col("start_date") - pl.duration(days=buffer_before_duration_in_days + padding_before_duration_in_days)).alias("buffer_before_start"),
            (pl.col("start_date") - pl.duration(days=padding_before_duration_in_days)).alias("padding_before_end"),
            (pl.col("end_date") + pl.duration(days=padding_after_duration_in_days)).alias("padding_after_start"),
            (pl.col("end_date") + pl.duration(days=buffer_after_duration_in_days + padding_after_duration_in_days)).alias("buffer_after_end"),
        ])
    )
    
    buffer_before = lf_extended.select([
        pl.col("id"),
        pl.col("buffer_before_start").alias("start_date"),
        pl.col("padding_before_end").alias("end_date"),
        pl.concat_list(pl.lit("buffer_before", dtype=region_types_enum)).alias("region_types")
#         pl.lit(["buffer_before"], dtype=region_types_enum).alias("region_types")
    ])
    
    hospital_stay = lf.select([
        pl.col("id"),
        pl.col("start_date"),
        pl.col("end_date"),
        pl.concat_list(pl.lit("APC", dtype=region_types_enum)).alias("region_types")
#         pl.lit(["APC"], dtype=region_types_enum).alias("region_types")
    ])
    
    buffer_after = lf_extended.select([
        pl.col("id"),
        pl.col("padding_after_start").alias("start_date"),
        pl.col("buffer_after_end").alias("end_date"),
        pl.concat_list(pl.lit("buffer_after", dtype=region_types_enum)).alias("region_types")
#         pl.lit(["buffer_after"], dtype=region_types_enum).alias("region_types")
    ])
    
    return (
        pl.concat([
            buffer_before, 
            hospital_stay, 
            buffer_after
        ])
        .sort(["id", "start_date"])
        .rename({
            "id":id_column
        })
        
    )


# In[ ]:


def split_overlapping_intervals_and_remerge( 
    lf: pl.LazyFrame,
    id_column: pl.Utf8 = "pseudo_nhs_number",
    start_date_column: pl.Utf8 = "start_date", 
    end_date_column: pl.Utf8 = "end_date",
) -> pl.LazyFrame:
    lf = lf.rename({
        id_column: "id",
        start_date_column: "start_date",
        end_date_column: "end_date",
    })
    
    return ( 
        lf
        .sort(["id", "start_date", "end_date"], descending=[False, False, True])
        .select([
            pl.col("id"),
            pl.col("start_date").alias("start_boundary"),
            pl.col("end_date").alias("end_boundary")
        ])
        .unpivot(index=["id"], on=["start_boundary", "end_boundary"], variable_name="type", value_name="boundary")
        .select(["id", "boundary"])
        .unique()
        .sort(["id", "boundary"])
        .with_columns(
            pl.col("boundary").shift(-1).over("id").alias("end_date")
        )
        .drop_nulls()
        .join_where(
            lf,
            pl.col("id").eq(pl.col("id_right")) &
            (pl.col("boundary") >= pl.col("start_date")) &
            (pl.col("end_date") <= pl.col("end_date_right")),
            suffix="_right"
        )
        .group_by(["id", "boundary", "end_date"], maintain_order=True) 
        .agg(pl.col("region_types").explode().unique().sort().alias("region_types"))
        .rename({"boundary": "start_date"})
        .with_columns(
            (pl.col("start_date") > pl.col("end_date").shift(1)).fill_null(True).alias("new_group") | 
            (pl.col("id") != pl.col("id").shift(1)).fill_null(True) |
            (pl.col("region_types") != pl.col("region_types").shift(1)).fill_null(True)
        )
        .with_columns(
            pl.col("new_group")
            .cum_sum()
            .alias("group")
        )
        .group_by(["id", "group", "region_types"], maintain_order=True)
        .agg(
            pl.col("start_date").min(),
            pl.col("end_date").max()
        )
        .select(
            pl.col("id").alias(id_column),
            pl.col("start_date").alias(start_date_column),
            pl.col("end_date").alias(end_date_column),
            pl.col("region_types"),
        )
    ) 


# ### Non-HES or general functions

# In[ ]:


def add_valid_test_date_from_candidate_columns(lf: pl.LazyFrame, date_cols: list[str]) -> pl.LazyFrame:
    """Adds a 'test_date' column selecting the first valid date from the provided date columns in order."""

    def safe_date_cast(column: pl.Expr) -> pl.Expr:
        """Attempts to cast a column to date; returns None if invalid."""
        return (
            pl.when(
                # Some years start with 4XXX, this corrects this so e.g. 4022 --> 2022
                # One this is corrected it attempts to cast to date assuming "%Y-%m-%d %H:%M" format
                column
                .str.replace("^4(\d{3})(.*?)$","2$1$2")
                .str.to_date(format="%Y-%m-%d %H:%M", strict=False).is_not_null()
            )
            .then(
                # Odd, strict = False should not be needed as by def only valid formats get here
                # however, without it, pipe produces a PanicExecution error
                column
                .str.to_date(format="%Y-%m-%d %H:%M", strict=False)  
            )
            .otherwise(
                None
            )
        )
    
    # Create validated date columns
    valid_cols = [safe_date_cast(pl.col(col)).alias(f"valid_{col}") for col in date_cols]
    
    return (
        lf.with_columns(valid_cols)
        .with_columns(pl.coalesce([f"valid_{col}" for col in date_cols]).alias("test_date"))  
    )
    


# In[ ]:


def display_with(df: pl.DataFrame, num_rows: int = 20, text_width: int = 80) -> None:
    """A helper function to allow fuller display of polar dataframes with limited trucation of """
    """number of rows and column width expansion"""
    df_num_rows = df.height
    if num_rows > 120:
        raise ValueError(f"Don't be daft, that's too big a num_rows ({num_rows}).  Try 120 rows or fewer.")
        return None
    if num_rows < df_num_rows:
        print(f"Too few num_rows ({num_rows}), ",end="")
        num_rows = min(df_num_rows,120)
        print(f"df is {df_num_rows} rows, setting num_rows to {num_rows}.")
    with pl.Config(tbl_rows=num_rows, fmt_str_lengths=text_width):
        display(df)
    return None


# In[ ]:


def free_variable_memory(variable_to_delete: str) -> str:
    if globals().get(variable_to_delete): 
        del globals()[variable_to_delete]
        gc.collect()
        return(f"Variable {variable_to_delete} found and deleted.")
    else:
        return(f"Variable {variable_to_delete} not found.")            


# ## Instantiate Pipeline paths

# In[ ]:


PIPELINE_HELPERS_PATH = (
    AnyPath(
        PIPELINE_HELPERS_LOCATION
    )
)


# In[ ]:


PIPELINE_INPUTS_PATH = (
    AnyPath(
        PIPELINE_INPUTS_LOCATION
    )
)


# In[ ]:


PIPELINE_OUTPUTS_PATH = (
    AnyPath(
        PIPELINE_OUTPUTS_LOCATION
    )
)


# In[ ]:


PIPELINE_OUTPUTS_REFERENCE_COMBO_FILES_PATH = (
    AnyPath(
        PIPELINE_OUTPUTS_PATH,
        "reference_combo_files"
    )
)


# In[ ]:


PIPELINE_LOGS_PATH = (
    AnyPath(
        PIPELINE_LOGS_LOCATION
    )
)


# In[ ]:


PIPELINE_INDIVIDUAL_TRAIT_FILES_PATH = (
    AnyPath(
        PIPELINE_INDIVIDUAL_TRAIT_FILES_LOCATION
    )
)


# In[ ]:


PIPELINE_INDIVIDUAL_TRAIT_PLOTS_PATH = (
    AnyPath(
        PIPELINE_INDIVIDUAL_TRAIT_PLOTS_LOCATION
    )
)


# In[ ]:


TRAIT_FEATURES_PATH = (
    AnyPath(
        *TRAIT_FEATURES_LOCATION
    )
)


# In[ ]:


TRAIT_ALIASES_LONG_PATH = (
    AnyPath(
        *TRAIT_ALIASES_LONG_LOCATION
    )
)


# In[ ]:


UNIT_CONVERSIONS_PATH = (
    AnyPath(
        *UNIT_CONVERSIONS_LOCATION
    )
)


# In[ ]:


MEGA_LINKAGE_PATH = (
    AnyPath(
        *MEGA_LINKAGE_LOCATION
    )
)


# In[ ]:


COMBO_POST_10D_WINDOWING_PATH = (
    AnyPath(
        *COMBO_POST_10D_WINDOWING_LOCATION
    )
)


# In[ ]:


# ADMISSIONS_DATA_PATH = (
#     AnyPath(
#         *ADMISSIONS_DATA_LOCATION
#     )
# )


# In[ ]:


# PROCESSED_COMBO_PATH = (
#     AnyPath(
#         *PROCESSED_COMBO_LOCATION
#     )
# )


# ## Dictionary of source files

# In[ ]:


source_files = {
    # primary care files
    "primary_care": [
        #f"{LIBRARY_RED_DATA_LOCATION}/DSA__Discovery_7CCGs/2025_XX_Discovery/gh3_observations.csv", # placeholder
        f"{LIBRARY_RED_DATA_LOCATION}/DSA__Discovery_7CCGs/2024_12_Discovery/gh3_observations.csv",
        f"{LIBRARY_RED_DATA_LOCATION}/DSA__Discovery_7CCGs/2024_07_Discovery/gh3_observations.csv",
        f"{LIBRARY_RED_DATA_LOCATION}/DSA__Discovery_7CCGs/2023_11_Discovery/gh3_observations.csv",
        f"{LIBRARY_RED_DATA_LOCATION}/DSA__Discovery_7CCGs/2023_03_Discovery/gh3_observations.csv",
        f"{LIBRARY_RED_DATA_LOCATION}/DSA__Discovery_7CCGs/2022_12_Discovery/GNH_thwfnech-phase2-outfiles_merge/cohort_gh2_observations_output_dataset_20221207.csv",
        f"{LIBRARY_RED_DATA_LOCATION}/DSA__Discovery_7CCGs/2022_12_Discovery/GNH_bhr-phase2-outfiles_merge/gh2_observations_dataset_20221207.csv",
        f"{LIBRARY_RED_DATA_LOCATION}/DSA__Discovery_7CCGs/2022_04_Discovery/GNH_thwfnech-phase2-outfiles_merge/GNH_thwfnech_observations_output_dataset_20220423.csv",
        f"{LIBRARY_RED_DATA_LOCATION}/DSA__Discovery_7CCGs/2022_04_Discovery/GNH_bhr-phase2-outfiles_merge/GNH_bhr_observations_output_dataset_20220412.csv",
        # Majority of NHSED data are secondary care but NDA is principally primary care
        f"{NHSE_SUBLICENSE_DATA_LOCATION}/DSA__NHSDigitalNHSEngland/2024_10/NDA/NIC338864_NDA_*.txt", #one of ["BMI","BP","CHOL","HBAIC","CORE"]        

    ],
    # repeat for secondary care 
    "secondary_care": [
        f"{LIBRARY_RED_DATA_LOCATION}/DSA__BartsHealth_NHS_Trust/2024_09_ResearchDataset/RDE_Pathology.ascii.nohisto.redacted2.csv",
        f"{LIBRARY_RED_DATA_LOCATION}/DSA__BartsHealth_NHS_Trust/2024_09_ResearchDataset/RDE_Measurements.ascii.redacted2.tab",
        f"{LIBRARY_RED_DATA_LOCATION}/DSA__BartsHealth_NHS_Trust/2023_12_ResearchDatasetv1.6/GandH_Measurements__20240423.ascii.redacted2.tab",
        f"{LIBRARY_RED_DATA_LOCATION}/DSA__BartsHealth_NHS_Trust/2023_12_ResearchDatasetv1.6/GH_Pathology__20231218.ascii.nohisto.redacted2.tab",
        f"{LIBRARY_RED_DATA_LOCATION}/DSA__BartsHealth_NHS_Trust/2023_05_ResearchDatasetv1.5/GandH_Measurements_202305151304.ascii.redacted.tab",
        f"{LIBRARY_RED_DATA_LOCATION}/DSA__BartsHealth_NHS_Trust/2023_05_ResearchDatasetv1.5/GH_Pathology_202305071651.ascii.redacted.nohisto.tab",
        f"{LIBRARY_RED_DATA_LOCATION}/DSA__BartsHealth_NHS_Trust/2022_03_ResearchDatasetv1.3/GandH_Pathology_202203191143_redacted_noHistopathologyReport.csv",
        f"{LIBRARY_RED_DATA_LOCATION}/DSA__BartsHealth_NHS_Trust/2021_04_PathologyLab/*.csv",
        f"{LIBRARY_RED_DATA_LOCATION}/DSA__BradfordTeachingHospitals_NHSFoundation_Trust/2022_06_BTHNFT/1578_gh_cerner_measurements_2022-06-10_redacted.tsv",
        f"{LIBRARY_RED_DATA_LOCATION}/DSA__BradfordTeachingHospitals_NHSFoundation_Trust/2023_05_BTHNFT/1578_gh_lab_results_2023-06-09_noCR.ascii.redacted.tab",
        f"{LIBRARY_RED_DATA_LOCATION}/DSA__BradfordTeachingHospitals_NHSFoundation_Trust/2024_12_BTHNFT/1578_gh_lab_results_2024-12-05.ascii.redacted.tab",
        f"{LIBRARY_RED_DATA_LOCATION}/DSA__BradfordTeachingHospitals_NHSFoundation_Trust/2024_12_BTHNFT/1578_gh_cerner_measurements_2024-12-05.ascii.redacted.tab",
    ],
}


# ## Set "ephemeral" paths
# 
# These are directories which contain file required for processing such as the raw data copied over from `/library-red/` as well as intermediary/temporary files created and used during the running of the pipeline (e.g. `combined_dataset` arrow files).
# 
# We may end up copying all to `/library-red/` or not

# In[ ]:


# This is the destination for raw data files copied from /library-red and/or /nhsdigital to 
# ivm 
PIPELINE_RAW_DATA_PATH = AnyPath(
    PIPELINE_DATA_LOCATION,
    "raw_datasets"
)


# In[ ]:


PRIMARY_ARROW_PATH = AnyPath(
    PIPELINE_DATA_LOCATION,
    "primary_care",
    "arrow"
)


# In[ ]:


NDA_ARROW_PATH = AnyPath(
    PIPELINE_DATA_LOCATION,
    "nda",
    "arrow"
)


# In[ ]:


SECONDARY_ARROW_PATH = AnyPath(
    PIPELINE_DATA_LOCATION,
    "secondary_care",
    "arrow"
)


# In[ ]:


COMBINED_DATASETS_ARROW_PATH = AnyPath(
    PIPELINE_DATA_LOCATION,
    "combined_datasets",
    "arrow"
)


# # Create all pipeline directories as needed
# These should all be needed on the first run, but subsequent runs if done do not re-create or delete the existing version.

# In[ ]:


PIPELINE_RAW_DATA_PATH.mkdir(parents=True, exist_ok=True)


# In[ ]:


PIPELINE_HELPERS_PATH.mkdir(parents=True, exist_ok=True)
PIPELINE_INPUTS_PATH.mkdir(parents=True, exist_ok=True)
PIPELINE_OUTPUTS_PATH.mkdir(parents=True, exist_ok=True)
PIPELINE_OUTPUTS_REFERENCE_COMBO_FILES_PATH.mkdir(parents=True, exist_ok=True)
PIPELINE_LOGS_PATH.mkdir(parents=True, exist_ok=True)

PIPELINE_INDIVIDUAL_TRAIT_FILES_PATH.mkdir(parents=True, exist_ok=True)
PIPELINE_INDIVIDUAL_TRAIT_PLOTS_PATH.mkdir(parents=True, exist_ok=True)


# In[ ]:


PRIMARY_ARROW_PATH.mkdir(parents=True, exist_ok=True)
NDA_ARROW_PATH.mkdir(parents=True, exist_ok=True)
SECONDARY_ARROW_PATH.mkdir(parents=True, exist_ok=True)
COMBINED_DATASETS_ARROW_PATH.mkdir(parents=True, exist_ok=True)


# ## Copy `raw_datasets` from `/library-red/` to `/home/ivm/`

# In[ ]:


# Set to True if copy from library-red is required
# For consorsium 2.0 use gcloud storage cp
PERFORM_COPY = False 


# In[ ]:


if PERFORM_COPY:
    for health_provider, files in source_files.items():
        destination_dir = AnyPath(PIPELINE_RAW_DATA_PATH)

        if not destination_dir:
            print(f'Destination not found for health provider "{health_provider}", skipping…')
            continue

        for file in tqdm(files, desc=f'Copying files for {health_provider}'):

            if '*' in file:
                matched_files = glob.glob(file)

                for file_path in matched_files:
                    try:
                        source_file = AnyPath(file_path)

                        sub_path = AnyPath(health_provider, *source_file.parts[-3:-1]) # get last two segments (exculding file name), and prepend health_provider
                        destination_dir = PIPELINE_RAW_DATA_PATH / sub_path 
                        destination_file = destination_dir / source_file.name

                        destination_dir.mkdir(parents=True, exist_ok=True)

                        # Check if file already exists in destination:
                        if destination_file.exists():
                            print(f'Skipped: {destination_file} already exists.')
                            continue

                        # Copy the file
                        shutil.copy(source_file, destination_file)
                        print(f'Copied [glob]:{source_file} -> {destination_file}') 
                    except Exception as e:
                        print(f'failed to copy {file}: {e}')

            else: 
                try:
#                     file = file.replace("/genesandhealth/library-red/genesandhealth/",
#                                  "gs://qmul-production-library-red/genesandhealth/"
#                                 )
                    source_file = AnyPath(file)

                    sub_path = AnyPath(health_provider, *source_file.parts[-3:-1]) # get last two segments (exculding file name), and prepend health_provider
                    destination_dir = PIPELINE_RAW_DATA_PATH / sub_path 
                    destination_file = destination_dir / source_file.name

                    destination_dir.mkdir(parents=True, exist_ok=True)

                    # Check if file already exists in destination:
                    if destination_file.exists():
                        print(f'Skipped: {destination_file} already exists.')
                        continue


                    # Copy the file
#                     print(f"Running: gcloud storage cp {source_file} {destination_file}")
#                     subprocess.run(f"gcloud storage cp {source_file} {destination_file}")
                    shutil.copy(source_file, destination_file)
                    print(f'Copied [single]:{source_file} -> {destination_file}')
                except Exception as e:
                      print(f'failed to copy {file}: {e}')


# ## Define column (sets) for output 

# In[ ]:


HASH_COLUMN = (
    pl.struct(
        [
            pl.col("pseudo_nhs_number"), 
            pl.col("test_date"), 
            pl.col("original_term"), 
            pl.col("result"), 
            pl.col("result_value_units")
        ]
    ).hash()
    .alias("hash")
)


# In[ ]:


TARGET_OUTPUT_COLUMNS = [
    pl.col("pseudo_nhs_number"),
    pl.col("test_date"),
    pl.col("original_term"),
    pl.col("result"),
    pl.col("result_value_units"),
    pl.col("provenance"),
    pl.col("source"),
]


# In[ ]:


TARGET_OUTPUT_COLUMNS_WITH_HASH = TARGET_OUTPUT_COLUMNS + [pl.col("hash")]


# In[ ]:


TARGET_TRAIT_LONG_COLUMNS = [
    pl.col("trait"),
    pl.col("target_units"),
    pl.col("min"),
    pl.col("max"),
    pl.col("alias"),
]


# In[ ]:


TARGET_MEGA_LINKAGE_COLUMNS = [
    pl.col("OrageneID"),
    pl.col("gender"),
    pl.col("gsa_id"),
    pl.col("exome_id"),
    pl.col("pseudo_nhs_number"),
    pl.col("FID"),
    pl.col("lane"),
]


# In[ ]:


TARGET_JOINED_LINK_FILE_AND_QUESTIONNAIRE_COLUMNS = [
    pl.col("OrageneID"), 
    pl.col("gender"), 
    pl.col("dob"), 
    pl.col("year_of_birth"), 
    pl.col("exome_id"), 
    pl.col("pseudo_nhs_number"), 
    pl.col("gsa_id"), 
    pl.col("FID"), 
    pl.col("lane"),
]


# In[ ]:


TARGET_TRAIT_RAW_ALL_COLUMNS = [
    pl.col("pseudo_nhs_number"),
    pl.col("trait"),
    pl.col("unit"),
    pl.col("value"),
    pl.col("date"),
    pl.col("gender"),
    pl.col("age_at_test"),
]


# In[ ]:


TARGET_COMBO_POST_10D_WINDOWING_COLUMNS = [
    pl.col("pseudo_nhs_number"),
    pl.col("trait"),
    pl.col("unit"),
    pl.col("value"),
    pl.col("date"),
    pl.col("gender"),
    pl.col("age_at_test"),
    pl.col("minmax_outlier"),
    pl.col("region_types")
]


# In[ ]:


TARGET_TRAIT_READINGS_AT_INDIVIDUAL_TIMEPOINTS_COLUMNS = [
    pl.col("pseudo_nhs_number"),
    pl.col("trait"),
    pl.col("unit"),
    pl.col("value"),
    pl.col("date"),
    pl.col("gender"),
    pl.col("age_at_test"),
    pl.col("minmax_outlier"),
]


# In[ ]:


TARGET_TRAIT_PER_INDIVIDUAL_STATS_COLUMNS = [
    pl.col("pseudo_nhs_number"),
    pl.col("trait"), 
    pl.col("median"),
    pl.col("mean"),
    pl.col("max"),
    pl.col("min"),
    pl.col("earliest"),
    pl.col("latest"),
    pl.col("n")
]


# ## Define filters
# 
# To be used in a polars `.filter` 

# ### Non-HES filters
# 
# Because polars filter **keeps** rows which match the criteria, if we want to exclude something we define the complement. 
# 
# For example, `EXCLUDE_NULL_UNITS` uses `.is_not_null()` which means only non_null values will be preserved; ergo, null units will be excluded, hence the naming convention.

# In[ ]:


EXCLUDE_NULL_UNITS = [
    pl.col("result_value_units").is_not_null()
]


# In[ ]:


EXCLUDE_READINGS_WITH_IMPLAUSIBLE_DATES = [
    pl.col("age_at_test") >= 0,
    pl.col("test_date") <= datetime.datetime.today(),
]


# In[ ]:


# G&H has volunteers from age 16+
EXCLUDE_READINGS_WITH_INDIVS_UNDER_SIXTEEN = [
    pl.col("age_at_test") >= 16
]


# In[ ]:


EXCLUDE_READINGS_WITH_VALUES_OUTSIDE_EXPECTED_RANGE = [
    pl.col("range_position").eq("ok")
]


# ### HES filters

# In[ ]:


IN_APC_ONLY = (
    pl.col("region_types").list.contains("APC"),
    ~pl.col("region_types").list.contains("buffer_before"),
    ~pl.col("region_types").list.contains("buffer_after"),
)

IN_APC_ANY = (
    pl.col("region_types").list.contains("APC"),
)

IN_BUFFER_BEFORE_ONLY = (
    pl.col("region_types").list.contains("buffer_before"),
    ~pl.col("region_types").list.contains("APC"),
    ~pl.col("region_types").list.contains("buffer_after"),
)

IN_BUFFER_BEFORE_ANY = (
    pl.col("region_types").list.contains("buffer_before"),
)

IN_BUFFER_AFTER_ONLY = (
    ~pl.col("region_types").list.contains("buffer_before"),
    ~pl.col("region_types").list.contains("APC"),
    pl.col("region_types").list.contains("buffer_after"),
)

IN_BUFFER_AFTER_ANY = (
    pl.col("region_types").list.contains("buffer_after"),
)

IN_BUFFERS_ONLY = (
    (
        pl.col("region_types").list.contains("buffer_before")
        | pl.col("region_types").list.contains("buffer_after")
    )
    & ~pl.col("region_types").list.contains("APC"),
)

IN_BUFFERS_ANY = (
    pl.col("region_types").list.contains("buffer_before")
    | pl.col("region_types").list.contains("buffer_after")
)

IN_TOTAL_EXCLUSION_ZONE = (
    pl.col("region_types").list.contains("APC") 
    | pl.col("region_types").list.contains("buffer_before") 
    | pl.col("region_types").list.contains("buffer_after"),
)

OUT_OF_APC = (
    ~pl.col("region_types").list.contains("APC")
    | pl.col("region_types").is_null()
)

OUT_OF_TOTAL_EXCLUSION_ZONE = (
    (
        ~pl.col("region_types").list.contains("APC") 
        & ~pl.col("region_types").list.contains("buffer_before")
        & ~pl.col("region_types").list.contains("buffer_after") 
    )
    | pl.col("region_types").is_null()
)


# In[ ]:


## RUN ALL ABOVE, Run all above, run all, run this
## A place holder, running cells above instantiates all packages/functions/path/etc. for pipeline


# In[ ]:


# Mystery red-dot!
# ​


# ## Primary Care

# ### Define primary care data paths
# 
# Because some branches to raw data can be longer than others, the globbing is defined here.

# In[ ]:


primary_care_paths = {
    "2022_04_Discovery_path": ('2022_04_Discovery', '*', '*'),
    "2022_12_Discovery_path": ('2022_12_Discovery', '*', '*'),
    "2023_03_Discovery_path": ('*', '2023_03_Discovery', '*'),
    "2023_11_Discovery_path": ('*', '2023_11_Discovery', '*'),
    "2024_07_Discovery_path": ('*', '2024_07_Discovery', '*'),
    "2024_12_Discovery_path": ('*', '2024_12_Discovery', '*'),
    # "2024_10_NHSDigitalNHSEngland_path": ('NDA','NIC338864_NDA_*','*.txt')
}


# ## Process all primary care .csv/.tsv/.tab to arrow 
# This section performs common transformation and tidy-up tasks on all primary care raw data and stores outputs as .arrow intermediates.
# 
# **Lengthy task (e.g. ~15 mins on 64 GB, 8 Cores, AMD processor)**
# 

# In[ ]:


get_ipython().run_cell_magic('time', '', 'SINK_PRIMARY_IPC=True # set to True if arrow regeneration of primary data is required.\n\nif SINK_PRIMARY_IPC:\n    for path_key, path_tuple in tqdm(primary_care_paths.items()):\n        print(f"{path_key}:")\n        (\n            pl.scan_csv(\n                AnyPath(PIPELINE_RAW_DATA_PATH, \'primary_care\', *path_tuple), \n                infer_schema=False,\n                null_values=["NULL"],\n            )\n            .filter(\n                pl.col("clinical_effective_date").is_not_null(),\n                pl.col("result_value").is_not_null(),\n                pl.col("result_value_units").is_not_null(),\n            )\n            .with_columns(\n                pl.col("original_code").cast(pl.Int64),\n                pl.col("clinical_effective_date").cast(pl.Date, strict=True).alias("test_date"),\n                pl.col("result_value").cast(pl.Float64, strict=True).alias("result"),\n                pl.col("original_term"),\n                provenance=pl.lit(path_key, pl.Enum(ALL_PROVENANCE_OPTIONS)),\n                source=pl.lit("primary_care", pl.Enum(ALL_SOURCE_OPTIONS)),\n            )\n            .with_columns(\n                HASH_COLUMN,\n            )\n            .sort("hash")\n            .unique(subset=["hash"])\n            .select(\n                *TARGET_OUTPUT_COLUMNS_WITH_HASH\n            )\n            .sink_ipc(\n                AnyPath(\n                    PRIMARY_ARROW_PATH, \n                    f"{path_key}.arrow")\n            )\n        )\n        print(f"  Arrow written.")\n')


# ### Import and concatenate all .arrow primary care data
# 
# We also include NDA data. At present due to their small size, these datasets are read directly from their Google cloud buckets. Intermediate .arrow files are written for NDA.

# In[ ]:


primary_24_arrow = (
    pl.scan_ipc(
        AnyPath(
            PRIMARY_ARROW_PATH, 
            "2024_*_Discovery_path.arrow"
        )
    )
)


# In[ ]:


primary_23_arrow = (
    pl.scan_ipc(
        AnyPath(
            PRIMARY_ARROW_PATH, 
            "2023_*_Discovery_path.arrow"
        )
    )
)


# In[ ]:


primary_22_arrow = (
    pl.scan_ipc(
        AnyPath(
            PRIMARY_ARROW_PATH, 
            "2022_*_Discovery_path.arrow"
        )
    )
)


# ### NHSE NDA

# In[ ]:


bmi = (
pl.scan_csv(
    AnyPath("/genesandhealth/nhsdigital-sublicence-red/DSA__NHSDigitalNHSEngland/2024_10/NDA/NIC338864_NDA_BMI.txt"),
    separator="|",
    )
    .filter(
        pl.col("BMI_VALUE").is_not_null()
    )
    .with_columns(
        #         STUDY_ID	BMI_DATE	AUDIT_YEAR	BMI_VALUE
        pl.col("STUDY_ID").alias("pseudo_nhs_number"),
        pl.col("BMI_DATE").str.to_date(format="%Y-%m-%d %H:%M:%S%.f").alias("test_date"),
        pl.col("BMI_VALUE").alias("result"),
        pl.lit("kg/m2").alias("result_value_units"),
        pl.lit("Body Mass Index Measured").alias("original_term"),
    )
    .select(
        pl.col("pseudo_nhs_number"),
        pl.col("test_date"),
        pl.col("original_term"),
        pl.col("result"),
        pl.col("result_value_units"),
    )
)

chol = (
pl.scan_csv(
    AnyPath("/genesandhealth/nhsdigital-sublicence-red/DSA__NHSDigitalNHSEngland/2024_10/NDA/NIC338864_NDA_CHOL.txt"),
    separator="|",
    )
    .filter(
        pl.col("CHOL_VALUE").is_not_null()
    )
    .with_columns(
        #         STUDY_ID	CHOLESTEROL_DATE	AUDIT_YEAR	CHOL_VALUE
        pl.col("STUDY_ID").alias("pseudo_nhs_number"),
        pl.col("CHOLESTEROL_DATE").str.to_date(format="%Y-%m-%d %H:%M:%S%.f").alias("test_date"),
        pl.col("CHOL_VALUE").alias("result"),
        pl.lit("mmol/L").alias("result_value_units"),
        pl.lit("Serum total cholesterol level").alias("original_term"),
    )
    .select(
        pl.col("pseudo_nhs_number"),
        pl.col("test_date"),
        pl.col("original_term"),
        pl.col("result"),
        pl.col("result_value_units"),
    )
)


hba1c = (
pl.scan_csv(
    AnyPath("/genesandhealth/nhsdigital-sublicence-red/DSA__NHSDigitalNHSEngland/2024_10/NDA/NIC338864_NDA_HBA1C.txt"),
    separator="|",
    )
    .filter(
        pl.col("HBA1C_MMOL_VALUE").is_not_null()
    )
    .with_columns(
        #         STUDY_ID	HBA1C_MMOL_VALUE	(AUDIT_YEAR)	(HBA1C_%_VALUE)	HBA1C_DATE	
        pl.col("STUDY_ID").alias("pseudo_nhs_number"),
        pl.col("HBA1C_DATE").str.to_date(format="%Y-%m-%d %H:%M:%S%.f").alias("test_date"),
        pl.col("HBA1C_MMOL_VALUE").alias("result"),
        pl.lit("mmol/mol").alias("result_value_units"),
        pl.lit("Haemoglobin A1c level").alias("original_term"),
    )
    .select(
        pl.col("pseudo_nhs_number"),
        pl.col("test_date"),
        pl.col("original_term"),
        pl.col("result"),
        pl.col("result_value_units"),
    )
)


bp = (
pl.scan_csv(
    AnyPath("/genesandhealth/nhsdigital-sublicence-red/DSA__NHSDigitalNHSEngland/2024_10/NDA/NIC338864_NDA_BP.txt"),
    separator="|",
    )
    .filter(
        pl.col("DIASTOLIC_VALUE").is_not_null(),
        pl.col("SYSTOLIC_VALUE").is_not_null()
    )
    .with_columns(
        #         STUDY_ID	(AUDIT_YEAR)	BP_Date	DIASTOLIC_VALUE	SYSTOLIC_VALUE	
        pl.col("STUDY_ID").alias("pseudo_nhs_number"),
        pl.col("BP_Date").str.to_date(format="%Y-%m-%d %H:%M:%S%.f").alias("test_date"),
        pl.col("DIASTOLIC_VALUE").alias("Diastolic arterial pressure"),
        pl.col("SYSTOLIC_VALUE").alias("Systolic arterial pressure"),
    )
    .unpivot(
         on=["Diastolic arterial pressure","Systolic arterial pressure"],
         index=["pseudo_nhs_number", "test_date"],
         variable_name="original_term",
         value_name="result"
        
    )
    .with_columns(
        pl.lit("mmHg").alias("result_value_units"),
    )
    
    .select(
        pl.col("pseudo_nhs_number"),
        pl.col("test_date"),
        pl.col("original_term"),
        pl.col("result"),
        pl.col("result_value_units"),
    )
)

(
    pl.concat([
        bmi,
        chol,
        hba1c,
        bp

    ])
    .with_columns(
        provenance=pl.lit("2024_10_NHSD_NHSE_NDA_path", pl.Enum(ALL_PROVENANCE_OPTIONS)),
        source=pl.lit("primary_care", pl.Enum(ALL_SOURCE_OPTIONS)),
    )
    .with_columns(
        HASH_COLUMN,
    )
    .unique(subset=["hash"])
    .sink_ipc(
        AnyPath(
            NDA_ARROW_PATH,
            f"{yr}_{mon}_formatted_nda.arrow"
        )
    )
)


# ## Primary Combining

# In[ ]:


nda_combined = (
    pl.scan_ipc(
    AnyPath(
        NDA_ARROW_PATH,
        f"{yr}_{mon}_formatted_nda.arrow"
        )
    )
)


# In[ ]:


get_ipython().run_cell_magic('time', '', '(\n    pl.concat([\n    primary_22_arrow\n        .unique(pl.col("hash")),\n    (\n        pl.concat([\n        primary_23_arrow\n        .unique(pl.col("hash")),\n    primary_24_arrow\n        .unique(pl.col("hash")), # up to here: 4.6 GB\n        ])\n        .unique(pl.col("hash"))\n    ),\n    nda_combined\n        .unique(pl.col("hash"))\n    ])\n    .unique(pl.col("hash")) \n    \n    .sink_ipc(\n        AnyPath(\n            COMBINED_DATASETS_ARROW_PATH,\n            f"{yr}_{mon}_Combined_primary_care.arrow"\n        )\n    )\n)\n')


# ## Secondary Care

# ### Bradford
#  "2022_06_Bradford_measurements"  
#  "2023_05_Bradford_path"  
#  "2024_12_Bradford_measurements"  
#  "2024_12_Bradford_path"

# #### Bradford Pathology
# "2023_05_Bradford_path"  
# "2024_12_Bradford_path"

# ##### `2023_05_Bradford_path`

# In[ ]:


provenance_key = "2023_05_Bradford_path"
(
    pl.scan_csv(
        AnyPath(PIPELINE_RAW_DATA_PATH, "secondary_care", '*', '*', "1578_gh_lab_results_2023-06-09_noCR.ascii.redacted.tab"),
        infer_schema=False,
        separator='\t',
    )
    
    .with_columns(
        pl.col("lab_test_performed_date").cast(pl.Date, strict=True),
    )
    .rename({
        "PseudoNHS_2023_04_24":"pseudo_nhs_number",
        "lab_test_performed_date":"test_date",
        "EVENT_DESCRIPTION":"original_term",
        "RESULT":"result",
        "RESULT_UNIT_DESC":"result_value_units",
    })
    
    .filter( # FILTER 1: these represent major filters
        pl.col("result").is_not_null(),
    )
    
    .filter( # FILTER 2: No recoverable number in `result`
        pl.col("result").ne("NA"),
        pl.col("result").ne("N/A"),
        pl.col("result").ne("NA;INS"),
        pl.col("result").ne("Error"),
        pl.col("result").ne("High"),
        pl.col("result").ne(";INS"),
        pl.col("result").ne("Negative"),
        pl.col("result").ne("TNP"),
        pl.col("result").ne("See Film Comms."),
        ~pl.col("result").str.contains("(?i)detected"),
        ~pl.col("result").str.contains("(?i)positive"),
        ~pl.col("result").str.contains("(?i)see comment"),
        ~pl.col("result").str.contains("(?i)unable to process"),
    )

    .with_columns(
        pl.col("result")
            .str.strip_prefix("less thn ")
            .str.strip_prefix("Less thn ")
            .str.strip_prefix("Lss thn ") ## NOT IN DATA
            .str.strip_prefix("Less thnn ") ## NOT IN DATA
            .str.strip_prefix("Less than ") 
            .str.strip_prefix("Less Thn ")
            .str.strip_prefix("Greater than ")
            .str.strip_prefix("Grtr thn ")
            .str.strip_prefix(" ") 
            .str.strip_prefix("<")
            .str.strip_prefix(">"),
            
        provenance=pl.lit(provenance_key, pl.Enum(ALL_PROVENANCE_OPTIONS)),
        source=pl.lit("secondary_care", pl.Enum(ALL_SOURCE_OPTIONS)),
    )
    
    .filter(
        pl.col("result").ne("")
    )
    .with_columns(
        pl.col("result").cast(pl.Float64, strict=True)
    )
    .with_columns(
        HASH_COLUMN
    )
    .unique(subset=["hash"])
    .select(
        *TARGET_OUTPUT_COLUMNS_WITH_HASH
    )
    
    .sink_ipc(
        AnyPath(
            SECONDARY_ARROW_PATH,
            f"{provenance_key}.arrow"
        )
    )

)


#  ##### `2024_12_Bradford_path`

# In[ ]:


get_ipython().run_cell_magic('time', '', 'provenance_key = "2024_12_Bradford_path"\n(\n    pl.scan_csv(\n        AnyPath(PIPELINE_RAW_DATA_PATH, \'secondary_care\', \'*\', \'*\', \'1578_gh_lab_results_2024-12-05.ascii.redacted.tab\'),\n        infer_schema=False,\n        separator=\'\\t\',\n    )\n    \n    .with_columns(\n        pl.col("lab_test_performed_date").cast(pl.Date, strict=True),\n        provenance=pl.lit(provenance_key, pl.Enum(ALL_PROVENANCE_OPTIONS)),\n        source=pl.lit("secondary_care", pl.Enum(ALL_SOURCE_OPTIONS)),\n    )\n    .rename({\n        "PseudoNHS_2024-07-10":"pseudo_nhs_number",\n        "lab_test_performed_date":"test_date",\n#         "ORDER_ID":"original_code",\n        "EVENT_DESCRIPTION":"original_term",\n        "RESULT":"result",\n        "RESULT_UNIT_DESC":"result_value_units",\n    })\n    .TRE\n    .filter_with_logging(\n        ~pl.col("result").str.contains("-No evidence of past infection."),\n        label="Exclude rows where result = \'-No evidence of past infection.\'"\n    )\n    .TRE\n    .filter_with_logging(\n        pl.col("result").is_not_null(),\n        label="Exclude rows where result is null"\n    )\n\n    .with_columns(\n        pl.col("test_date").cast(pl.Date, strict=True),\n        pl.col("result")\n            .str.strip_prefix("less thn ")\n            .str.strip_prefix("Less thn ")\n            .str.strip_prefix("Less than") \n            .str.strip_prefix("Less Thn ")\n            .str.strip_prefix("Greater than ")\n            .str.strip_prefix("Grtr thn ")\n            .str.strip_prefix(" ")\n            .str.strip_prefix("<")\n            .str.strip_prefix(">")\n            .str.strip_prefix("NA")\n            .str.strip_prefix("N/A")\n            .str.strip_prefix("Error")\n            .str.strip_prefix("High")\n            .str.strip_prefix(";INS")\n            .str.strip_prefix("Negative")\n            .str.strip_prefix("Positive")\n            .str.strip_prefix("POSITIVE")\n            .str.strip_prefix("TNP")\n            .str.strip_prefix("See Film Comms.")\n            .str.replace("(?i)detected","")\n            .str.replace("(?i)see comment","")\n            .str.replace("(?i)unable to process","")\n            .str.replace("not ","")\n            .str.replace("Not ","")\n            .str.replace("NOT ",""),\n        provenance=pl.lit(provenance_key, pl.Enum(ALL_PROVENANCE_OPTIONS)),\n        source=pl.lit("secondary_care", pl.Enum(ALL_SOURCE_OPTIONS)),\n    )\n    .with_columns(\n        HASH_COLUMN\n    )\n    .unique(subset=["hash"])\n    .select(\n        *TARGET_OUTPUT_COLUMNS_WITH_HASH\n    )\n    .filter(\n        pl.col("result").ne("")\n    )\n    .with_columns(\n        pl.col("result").cast(pl.Float64, strict=True)\n    )\n#     .collect()\n    .sink_ipc(\n        AnyPath(\n            SECONDARY_ARROW_PATH, \n            f"{provenance_key}.arrow"\n        )\n    )\n)\n')


# #### Combine Bradford Pathology Data

# In[ ]:


get_ipython().run_cell_magic('time', '', '\n(\n    pl.scan_ipc(\n        [\n            AnyPath(\n                SECONDARY_ARROW_PATH,\n                "2023_05_Bradford_path.arrow"\n            ),\n            AnyPath(\n                SECONDARY_ARROW_PATH,\n                "2024_12_Bradford_path.arrow"\n            )\n        ]\n        \n    )\n   .with_columns(\n         HASH_COLUMN\n    )\n    .unique("hash")\n \n \n    .sink_ipc(\n        AnyPath(\n            SECONDARY_ARROW_PATH,\n            f"{yr}_{mon}_Bradford_path_combined.arrow"\n        )\n    )\n)\n')


# #### Bradford Measurements
#  "2022_06_Bradford_measurements"  
#  "2024_12_Bradford_measurements"  

# #####  `2022_06_Bradford_measurements`  
# 
# Available traits: `"Body Mass Index Estimated", "Body Mass Index Measured", "Height/Length Estimated", "Height/Length Measured", "Ideal Body Weight Calculated", "Patient Stated Weight", "Weight Estimated", "Weight Measured"`

# In[ ]:


get_ipython().run_cell_magic('time', '', 'provenance_key="2022_06_Bradford_measurements"\n\n# Columns in this dataset:\n# ["PseudoNHS","age_at_measurement","date_of_measurement","EVENT_CD","EVENT_TITLE","EVENT_ANSWER"]\n# Note: no units column\n\n\n(\n    pl.scan_csv(\n        AnyPath(PIPELINE_RAW_DATA_PATH, \'secondary_care\', \'*\', \'*\', \'1578_gh_cerner_measurements_2022-06-10_redacted.tsv\'),\n        infer_schema=False,\n        separator=\'\\t\',\n    )\n    \n    .with_columns(\n        pl.col("EVENT_ANSWER").cast(pl.Float64),\n        pl.col("date_of_measurement").str.to_date(format="%d/%m/%Y"),\n        pl.col("EVENT_CD").cast(pl.Int64),\n        pl.when(pl.col("EVENT_TITLE").str.contains(r"(?i)weight"))\n        .then(pl.lit("kg"))\n        .when(pl.col("EVENT_TITLE").str.contains(r"(?i)height"))\n        .then(pl.lit("cm"))\n        .when(pl.col("EVENT_TITLE").str.contains(r"(?i)index"))\n        .then(pl.lit("kg/m2")) # BMI unit\n        .otherwise(None) #\xa0Default case\n        .alias("result_value_units"),\n        provenance=pl.lit(provenance_key, pl.Enum(ALL_PROVENANCE_OPTIONS)),\n        source=pl.lit("secondary_care", pl.Enum(ALL_SOURCE_OPTIONS)),\n        \n    )\n    .rename({\n        "PseudoNHS":"pseudo_nhs_number",\n        "date_of_measurement":"test_date",\n        "EVENT_TITLE":"original_term",\n        "EVENT_ANSWER":"result",\n    })\n    .with_columns(\n        HASH_COLUMN\n    )\n    .unique(subset=[\'hash\'])\n    .select(\n       *TARGET_OUTPUT_COLUMNS_WITH_HASH\n    )\n    .sink_ipc(\n        AnyPath(\n            SECONDARY_ARROW_PATH, \n            f"{provenance_key}.arrow")\n    )\n)\n')


#  ##### `2024_12_Bradford_measurements`  
#  
#  Potentially available: `"AVPU Conscious Level", "Any Supplemental Oxygen", "Blood Glucose, Capillary", "Body Mass Index Estimated", "Body Mass Index Measured", "Date\Time Correction", "Diastolic Blood Pressure", "EWS Category", "EWS Status", "EWS Total", "EWS Type", "Escalated score to whom?", "FiO2", "Heart Rate Monitored", "Height/Length Estimated", "Height/Length Measured", "Ideal Body Weight Calculated", "Inspired O2", "Mean Arterial Pressure, Cuff", "Oxygen Flow Rate", "Oxygen Saturation Target", "Oxygen Therapy", "Patient Stated Weight", "Respiratory Distress", "Respiratory Rate", "SBP/DBP Cuff Locations", "SpO2", "Systolic Blood Pressure", "Temperature", "Weight Estimated", "Weight Measured"`

# In[ ]:


get_ipython().run_cell_magic('time', '', 'provenance_key="2024_12_Bradford_measurements"\n\n(\n    pl.scan_csv(\n        AnyPath(PIPELINE_RAW_DATA_PATH, \'secondary_care\', \'*\', \'*\', \'1578_gh_cerner_measurements_2024-12-05.ascii.redacted.tab\'),\n        infer_schema_length=0,\n        separator=\'\\t\',\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("EVENT_ANSWER").str.contains(" - "),\n        label="result contains \' - \'"\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("EVENT_ANSWER").str.contains(r"[a-zA-Z/]"),\n        label="EVENT_ANSWER.str.contains(r\'[a-zA-Z/]"\n    )\n    .with_columns(\n        pl.col("EVENT_ANSWER").cast(pl.Float64), \n        pl.col("date_of_measurement").cast(pl.Date),#str.to_date(format="%d/%m/%Y")\n        pl.when(pl.col("EVENT_TITLE").str.contains(r"(?i)weight"))\n        .then(pl.lit("kg"))\n        .when(pl.col("EVENT_TITLE").str.contains(r"(?i)height"))\n        .then(pl.lit("cm"))\n        .when(pl.col("EVENT_TITLE").str.contains(r"(?i)index"))\n        .then(pl.lit("kg/m^2")) # BMI unit\n        .when(pl.col("EVENT_TITLE").str.contains(r"(?i)pressure"))\n        .then(pl.lit("mmHg")) # BP unit\n        .when(pl.col("EVENT_TITLE").str.contains(r"(?i)glucose"))\n        .then(pl.lit("mmol/L")) # Glucose unit\n        .otherwise(None) #\xa0Default case\n        .alias("result_value_units"),\n        provenance=pl.lit(provenance_key, pl.Enum(ALL_PROVENANCE_OPTIONS)),\n        source=pl.lit("secondary_care", pl.Enum(ALL_SOURCE_OPTIONS)),\n        \n    )\n    .rename({\n        "PseudoNHS_2024-07-10":"pseudo_nhs_number",\n        "date_of_measurement":"test_date",\n        "EVENT_TITLE":"original_term",\n        "EVENT_ANSWER":"result",\n    })\n\n    .with_columns(\n        HASH_COLUMN\n    )\n    .unique(pl.col("hash"))\n    .select(\n       *TARGET_OUTPUT_COLUMNS_WITH_HASH\n    )\n    .sink_ipc(\n        AnyPath(\n            SECONDARY_ARROW_PATH,\n            f"{provenance_key}.arrow"\n        )\n    )\n)\n\n')


# #### Combine Bradford Measurement data

# In[ ]:


get_ipython().run_cell_magic('time', '', '(\n    pl.scan_ipc(\n        AnyPath(\n            SECONDARY_ARROW_PATH,\n            "*_Bradford_measurements.arrow"\n        )\n    )\n    .with_columns(\n         HASH_COLUMN\n    )\n    .unique("hash")\n \n    .sink_ipc(\n        AnyPath(\n            SECONDARY_ARROW_PATH,\n            f"{yr}_{mon}_Bradford_measurements_combined.arrow"\n        )\n    )\n)\n')


# ## Barts

# ### Barts Pathology
# Note there are five separate datasets here:
# - April 2021
# - March 2022
# - May 2023
# - Dec 2023
# - Sept 2024

# #### `Barts 2021 04 - April 2021`
# 

# ##### Set aside some of the files
# 
# ```
#  2021_01_25_pseudoNHS_uniq.csv
#  - .csv file in directory but not raw data file
#  
#  LipoproteinA_April2021.csv
#  - No units - however can we invent these?  
#  
#  AntiMullerianHormone_April2021.csv
#  -  Date unrecoverable
#  
#  MCH_April2021.csv  
#  - no date, just age at test  
# 
#  Haemoglobin_April2021.csv
#  - no date, just age at test  
#  
#  Progesterone_April2021.csv
#  - no date, just age at test  
#  
#  'Islet Antibody.csv'
#  - non-numerical result  
#  
#  RDW_April2021.csv
#  - no date, just age at test  
# 
# ```
# Check to see if they contain recoverable data.

# The following files have the "wrong" number of fields, should be 6 (actual shown after filename) and throw the code:
# * `Haemoglobin_April2021.csv` 4
# * `LipoproteinA_April2021.csv` 5
# * `MCH_April2021.csv` 5
# * `Progesterone_April2021.csv` 5
# * `RDW_April2021.csv` 5
# * `Random Glucose..csv 7` \[note double '.'\] NB. This is because of a trailing comma which does not impact processing to this file so it is **not** considered a problem file to exclude.

# #### strip `<` from result
# - DHEA sulphate level    
# - Basophils
# 
# #### strip `>` from result
# - GAD Antibodies
# 
# #### result_value_units `is_null()`
# - MCV_April2021

# In[ ]:


problem_files = [
    '2021_01_25_pseudoNHS_uniq.csv', # not a results file **
    'Haemoglobin_April2021.csv', # 4 fields
    'LipoproteinA_April2021.csv', # 5 fields
    'MCH_April2021.csv', # 5 fields **
    'Progesterone_April2021.csv', # 5 fields
    'RDW_April2021.csv', # 5 fields
    'AntiMullerianHormone_April2021.csv', # correct number of fields but data unrecoverable (^d{2}:\d{2}\.\d)
    'Islet Antibody.csv', # non-numerical result
]


# In[ ]:


Barts_2021_04_admissible_files = [file for file in
    AnyPath(
            PIPELINE_RAW_DATA_PATH, 
            "secondary_care", 
            "DSA__BartsHealth_NHS_Trust", 
            "2021_04_PathologyLab",          
        ).glob("*.csv")
 if file.name not in problem_files
]


# In[ ]:


get_ipython().run_cell_magic('time', '', '\nprovenance_key = "2021_04_Barts_path"\n\n(\n    pl.scan_csv(\n        Barts_2021_04_admissible_files,\n        has_header=False,\n        skip_lines=1,\n        infer_schema=False,\n        new_columns=[\n            "pseudo_nhs_number",\n            "column_2",\n            "original_term",\n            "test_date",\n            "result",\n            "result_value_units",\n        ],\n        include_file_paths="file",\n        null_values=["NULL"] # Basophils\n    )\n    .filter(\n        ## Basophils and Fasting Glucose files have single rows with errors which trip casting\n        ## to float so, regretably, we need bespoke filters here\n        pl.col("result").ne("1429 at 10.40 on 28/11/14."), # Basophils: filter out 1 row\n        pl.col("result").ne("08/01/2014"), # "Fasting Glucose." filter out 1 row\n    )\n    .with_columns(\n        pl.col("file").str.strip_suffix(".csv").str.split("/").list.last()\n    )\n    .with_columns(\n        pl.col("test_date").str.to_date(format="%d/%m/%Y"),\n        pl.col("result")\n        .str.strip_prefix("<")\n        .str.strip_prefix(">")\n        .str.strip_prefix(" ")\n        .cast(pl.Float64),\n        provenance=pl.lit(provenance_key, pl.Enum(ALL_PROVENANCE_OPTIONS)),\n        source=pl.lit("secondary_care", pl.Enum(ALL_SOURCE_OPTIONS)),\n    )\n    .with_columns(\n        HASH_COLUMN\n    )\n    .unique("hash")\n    .select(\n        TARGET_OUTPUT_COLUMNS_WITH_HASH\n    )\n    .sink_ipc(\n        AnyPath(\n            SECONDARY_ARROW_PATH,\n            f"{provenance_key}.arrow"\n        )\n    )\n)\n')


# #### `Barts 2022 03  - March 2022`
# 

# #### Some important comments about this file
# 
# The `.../2022_03_ResearchDatasetv1.3/GandH_Pathology_202203191143_redacted_noHistopathologyReport.csv` contains comma separated fields in double-quotes which themselves contain commas, this causes problems.
# 
# There are two potential solutions.
# 
# ##### Solution 1: quote_char = None
# 
# If we add this switch to scan_csv, it manages to import the whole file (6,640,330 rows).  However, the behaviour is that anytime a comma is seen, a new field is created.  Depending on the fields this might lead to unreliable/inconsistent column content and/or throw downstream processing.
# 
# ##### Solution 2: pre-process file to exclude problematic rows
# 
# Here we use a linux grep to exclude all rows which contain a double-quote followed by a comma before the next double-quote (i.e. double-quote bound fields with one or more commas in them):
# 
# `grep -Ev ',\"[^"]*,' GandH_Pathology_202203191143_redacted_noHistopathologyReport.csv > GandH_Pathology_202203191143_redacted_noHistopathologyReport.no_unmatched_double_quotes.csv`
# 
# This excludes 69,117 (1.04%) rows **but the behaviour is consistent and understood**.
# 
# We therefore use **Solution 2**.
# 

# #### HARD-CODED PRE-PROCESSING: BARTS_2022_03 (SOLUTION 2)

# In[ ]:


BARTS_2022_03_PATH = AnyPath(
    PIPELINE_RAW_DATA_PATH, 
    "secondary_care", 
    "DSA__BartsHealth_NHS_Trust", 
    "2022_03_ResearchDatasetv1.3",
)

# Input file
BARTS_2022_03_PATHOLOGY_FILE_RAW_PATH = AnyPath(
    BARTS_2022_03_PATH,
    "GandH_Pathology_202203191143_redacted_noHistopathologyReport.csv"
)

# Output file
BARTS_2022_03_PATHOLOGY_FILE_CORRECTED_PATH = AnyPath(
    BARTS_2022_03_PATH,
    "GandH_Pathology_202203191143_redacted_noHistopathologyReport.no_unmatched_double_quotes.csv"
)


# In[ ]:


barts_2022_03_preprocessing_command = (
    f"""grep -Ev ',\"[^"]*,' """
    f"""\"{BARTS_2022_03_PATHOLOGY_FILE_RAW_PATH}\" > """
    f"""\"{BARTS_2022_03_PATHOLOGY_FILE_CORRECTED_PATH}\""""
)


# In[ ]:


subprocess.run(
    barts_2022_03_preprocessing_command,
    shell=True,
    check=True,
    capture_output=True,
    text=True
)


# In[ ]:


get_ipython().run_cell_magic('time', '', '# Here we use the preprocessed file generated above\nprovenance_key = "2022_03_Barts_path"\n\n(\n    pl.scan_csv(\n        BARTS_2022_03_PATHOLOGY_FILE_CORRECTED_PATH,\n        infer_schema=False,\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("ResultTxt").str.contains(r"[a-zA-Z]"),\n        ~pl.col("ResultTxt").str.ends_with(" -"),\n        ~pl.col("ResultTxt").str.contains(r"\\d/\\d"),\n        ~pl.col("ResultTxt").str.contains("\\d{2}:\\d{2}"),\n        ~pl.col("ResultTxt").str.contains("\\++"),\n        ~pl.col("ResultTxt").str.contains(r"-+"),\n        ~pl.col("ResultTxt").str.contains("\\*+"),\n        ~pl.col("ResultTxt").str.contains("\\?"),\n        ~pl.col("ResultTxt").str.contains("\\("),\n        ~pl.col("ResultTxt").str.contains("\\d \\d"),\n        ~pl.col("ResultTxt").str.starts_with(" "),\n        pl.col("ResultTxt").ne("."),\n        pl.col("ResultTxt").ne("#"),\n        pl.col("ResultTxt").ne("]"),\n        pl.col("ResultTxt").ne("*"),\n        pl.col("ResultTxt").ne(":"),\n        pl.col("ResultTxt").ne("?"),\n        pl.col("ResultTxt").ne(". ."),\n        pl.col("ResultTxt").ne(". . . . ."),\n        pl.col("ResultTxt").ne("0.18*"),  \n        pl.col("ResultTxt").ne("22.01.15; 1800"),\n        label="Exclude non-numerical ResultsTxt",\n    )\n    .with_columns(\n        pl.col("ResultTxt")\n            .str.strip_prefix("< ")\n            .str.strip_prefix("<")\n            .str.strip_prefix(">")\n            .cast(pl.Float64, strict=True)\n            .alias("result"),\n        pl.col("ReportDate").str.to_date(format="%Y-%m-%d %H:%M", strict=True).alias("test_date"),\n        pl.col("PseudoNHSnumber").alias("pseudo_nhs_number"),\n        pl.col("TestDesc").alias("original_term"),\n        pl.col("ResultUnit").alias("result_value_units"),\n        provenance=pl.lit(provenance_key, pl.Enum(ALL_PROVENANCE_OPTIONS)),\n        source=pl.lit("secondary_care", pl.Enum(ALL_SOURCE_OPTIONS)),\n    )\n    .with_columns(\n        HASH_COLUMN\n    )\n    .unique(subset=["hash"]) \n    .select(\n        TARGET_OUTPUT_COLUMNS_WITH_HASH\n    )    \n    .sink_ipc(\n         AnyPath(\n             SECONDARY_ARROW_PATH,\n             f"{provenance_key}.arrow"\n         )\n     )\n)\n')


# #### `Barts 2023 05 - May 2023`

# ##### Pre-processing etc. for this file
# 
# ###### Solution 1:
# Here we use the strategy of rejecting quotes as a (text-)field delimiter (cf. separator) by setting `quote_char` to `|`. The separator itself is a tab.  This yields 86 more rows than solution 2 (including '"Haptoglobin' and '"Beta Trace Protein'.  However, the method is unreliable sometimes causing column right-shift.
# 
# ###### Solution 2:
# Pre-process the file in shell to remove problematic lines.
# `grep -Ev '<tab>\"[^"]*<tab>' GH_Pathology_202305071651.ascii.redacted.nohisto.tab > GH_Pathology_202305071651.ascii.redacted.nohisto.no_unmatched_double_quotes.tab` where \<tab\> is obtained in shell by typing Crtl-V followed by pressing the tab key
# 
# **We have chosen solution 2** since there are only 86 problem lines in a file of 12_390_030 total lines.

# #### HARD-CODED PRE-PROCESSING: BARTS_2023_05 (SOLUTION 2)

# In[ ]:


BARTS_2023_05_PATH = AnyPath(
    PIPELINE_RAW_DATA_PATH, 
    "secondary_care",
    "DSA__BartsHealth_NHS_Trust", 
    "2023_05_ResearchDatasetv1.5",
)

# Input file
BARTS_2023_05_PATHOLOGY_FILE_RAW_PATH = AnyPath(
    BARTS_2023_05_PATH,
    "GH_Pathology_202305071651.ascii.redacted.nohisto.tab"
)

# Output file
BARTS_2023_05_PATHOLOGY_FILE_CORRECTED_PATH = AnyPath(
    BARTS_2023_05_PATH,
    "GH_Pathology_202305071651.ascii.redacted.nohisto.no_unmatched_double_quotes.tab"
)


# In[ ]:


barts_2023_05_preprocessing_command = (
    f"""grep -Ev '\t\"[^"]*\t' """
    f"""\"{BARTS_2023_05_PATHOLOGY_FILE_RAW_PATH}\" > """
    f"""\"{BARTS_2023_05_PATHOLOGY_FILE_CORRECTED_PATH}\""""
)


# In[ ]:


subprocess.run(
    barts_2023_05_preprocessing_command,
    shell=True,
    check=True,
    capture_output=True,
    text=True
)


# In[ ]:


get_ipython().run_cell_magic('time', '', 'provenance_key = "2023_05_Barts_path"\n(\n    pl.scan_csv(\n        BARTS_2023_05_PATHOLOGY_FILE_CORRECTED_PATH,\n        separator="\\t",\n        infer_schema=False,\n        )\n\n    .filter(\n        pl.col("ResultTxt").ne("**"),\n        pl.col("ResultTxt").ne("***"),\n        pl.col("ResultTxt").ne("****"),\n        pl.col("ResultTxt").ne("*"),\n        pl.col("ResultTxt").ne("* -"),\n        pl.col("ResultTxt").ne("-"),\n        pl.col("ResultTxt").ne("--"),\n        pl.col("ResultTxt").ne("- -"),\n        pl.col("ResultTxt").ne("-  -"),\n        pl.col("ResultTxt").ne("+"), # present in 2023_05\n        pl.col("ResultTxt").ne("++"), # present in 2023_05\n        pl.col("ResultTxt").ne("+++"), # present in 2023_05\n        pl.col("ResultTxt").ne("++++"), # present in 2023_05\n        pl.col("ResultTxt").ne("*115"), # present in 2023_05\n        pl.col("ResultTxt").ne("#"),\n        pl.col("ResultTxt").ne("/"),\n        pl.col("ResultTxt").ne("`"),\n        pl.col("ResultTxt").ne(",."),\n        pl.col("ResultTxt").ne("."),\n        pl.col("ResultTxt").ne("....."),\n        pl.col("ResultTxt").ne("n/r"),\n        pl.col("ResultTxt").ne("na"),\n        pl.col("ResultTxt").ne("n/a"),\n        pl.col("ResultTxt").ne("NA"),\n        pl.col("ResultTxt").ne("?"),\n        pl.col("ResultTxt").ne(","),\n        pl.col("ResultTxt").ne(":"),\n        pl.col("ResultTxt").ne("]"),\n        pl.col("ResultTxt").ne("c"),\n        pl.col("ResultTxt").ne("MK"),\n        pl.col("ResultTxt").ne("B"),\n        pl.col("ResultTxt").ne("P"),\n        pl.col("ResultTxt").ne("ns"),\n        pl.col("ResultTxt").ne("1a"),\n        pl.col("ResultTxt").ne("1b"),\n        pl.col("ResultTxt").ne("3a"),\n        pl.col("ResultTxt").ne("3b"),\n        pl.col("ResultTxt").ne("3-"),\n        pl.col("ResultTxt").ne("64-"),\n        pl.col("ResultTxt").ne("B2A2"),\n        pl.col("ResultTxt").ne("B3A2"),\n        pl.col("ResultTxt").ne("FM"),\n        pl.col("ResultTxt").ne("UNS"),\n        pl.col("ResultTxt").ne("@unb"),\n        pl.col("ResultTxt").ne("@und"),\n        pl.col("ResultTxt").ne("None"),\n        pl.col("ResultTxt").ne("2-5"),\n        pl.col("ResultTxt").ne("1:8"),\n        pl.col("ResultTxt").ne("1:16"),\n        pl.col("ResultTxt").ne("1:32"),\n        pl.col("ResultTxt").ne("4o"),\n        pl.col("ResultTxt").ne("*40"),\n        pl.col("ResultTxt").ne("body"),\n        pl.col("ResultTxt").ne("Body"),\n        pl.col("ResultTxt").ne("24hr"),\n        pl.col("ResultTxt").ne("24HR"),\n        pl.col("ResultTxt").ne("KNIB"),\n        pl.col("ResultTxt").ne("64 -"),\n        pl.col("ResultTxt").ne("70)"),\n        pl.col("ResultTxt").ne("(70)"),\n        pl.col("ResultTxt").ne("(66"),\n        pl.col("ResultTxt").ne("*66"),\n        pl.col("ResultTxt").ne("*81"),\n        pl.col("ResultTxt").ne("*92"),\n        pl.col("ResultTxt").ne("*{88}"),\n        pl.col("ResultTxt").ne("*{94}"),\n        pl.col("ResultTxt").ne("5ml"),\n        pl.col("ResultTxt").ne("Serum"),\n        pl.col("ResultTxt").ne(" Serum\\""),\n        pl.col("ResultTxt").ne("clumps"),\n        pl.col("ResultTxt").ne("\\"Regret"),\n        pl.col("ResultTxt").ne("random"),\n        pl.col("ResultTxt").ne("Random"),\n        pl.col("ResultTxt").ne("RANDOM"),\n        pl.col("ResultTxt").ne("RAMDOM"),\n        pl.col("ResultTxt").ne("Clumped"),\n        pl.col("ResultTxt").ne("CLUMPED"),\n        pl.col("ResultTxt").ne("deleted"),\n        pl.col("ResultTxt").ne("DELETED"),\n        pl.col("ResultTxt").ne("Pending"),\n        pl.col("ResultTxt").ne("24 hour"), \n        pl.col("ResultTxt").ne("Not requested. PLEASE NOTE - THIS IS AN AMENDED REPORT"),\n        pl.col("ResultTxt").ne("No result available - see comment"),\n        pl.col("ResultTxt").ne("Not Calculated Units: mL/min/1.73sqm For Afro-Caribbean patients multiply eGFR by 1.21 Use with caution for adjusting drug dosage."),\n        pl.col("ResultTxt").ne("Intrinsic Factor antibodies not tested as Gastric Parietal Cell antibody was negative. http://jcp.bmj.com/content/62/5/439.abstract"),\n        pl.col("ResultTxt").ne("Albumin Creatinine ratio within normal limits"),\n        pl.col("ResultTxt").ne("Wrong patient bled. Suggest repeat."),\n        ~pl.col("ResultTxt").str.contains("^\\""),\n        ~pl.col("ResultTxt").str.contains("(?i)insufficient"),\n        ~pl.col("ResultTxt").str.contains("(?i)unsuitable"),\n        ~pl.col("ResultTxt").str.contains("(?i)inadequately"),\n        ~pl.col("ResultTxt").str.contains("(?i)received"),\n\n    )\n    .pipe(add_valid_test_date_from_candidate_columns, date_cols=["ReportDate","Report","RequestDate"])\n    .with_columns(\n        pl.col("PseudoNHS_2023_04_24").alias("pseudo_nhs_number"),\n        pl.col("TestDesc").alias("original_term"),\n        pl.col("ResultTxt")\n            .str.strip_prefix("<")\n            .str.strip_prefix(">")\n            .str.strip_prefix("+-") ## present in 2023_12\n            .str.strip_prefix("+/-") ## present in 2023_12\n            .str.replace(r"^\\{(.*?)\\}$","$1")\n            .str.strip_prefix(" ")\n            .str.strip_suffix(" -")\n            .str.strip_suffix("\\"")\n            .str.strip_suffix("%") # should spot check this since could be a typo (shift+5 instead of 5)\n            .str.strip_suffix(" g/l") # should spot check this\n            .cast(pl.Float64, strict=False)\n            .alias("result"),\n        pl.col("ResultUnit").alias("result_value_units"),\n        provenance=pl.lit(provenance_key, pl.Enum(ALL_PROVENANCE_OPTIONS)),\n        source=pl.lit("secondary_care", pl.Enum(ALL_SOURCE_OPTIONS)),\n        \n    )\n    .TRE\n    .filter_with_logging(\n        pl.col("test_date").is_not_null(),\n        label=\'Exclude null test_date\'\n    )\n    .TRE\n    .filter_with_logging(\n        pl.col("result").is_not_nan(),\n        label=\'Exclude result is nan\'\n    )\n\n    .with_columns(\n        HASH_COLUMN\n    )\n    .unique("hash")\n    \n    .select(\n        TARGET_OUTPUT_COLUMNS_WITH_HASH\n    )\n\n    .sink_ipc(\n        AnyPath(\n            SECONDARY_ARROW_PATH,\n            f"{provenance_key}.arrow")\n    )    \n)\n')


# #### `Barts 2023 12 - Dec 2023`

# #### Research dataset v1.6

# ##### Note about Research dataset v1.6 files
# 
# Note that the files is v2 of the redacted `GH_Pathology` file, hence `redacted2.tab`.
# 
# These files are ragged.  At some point we may consider using all raggedness but at present we run a script to subset the `GH_Pathology__20231218.ascii.nohisto.redacted2.tab` file into subfiles with the same number of tabs per line.
# 
# e.g. `GH_Pathology__20231218.ascii.nohisto.redacted2_tab16.tab`
# 
# To create this file is have a bash script (running an awk command) in `/home/ivm/QUANT_PY/versionXXX_YYYY_MM/helpers/` called `num_delims_splitter.sh`.  
# 
# The shell code for this can be found in the `Code graveyard` at the end of this notebook.

# #### HARD-CODED PRE-PROCESSING: BARTS_2023_12 (SPLIT BY NUMBER OF TABS)

# In[ ]:


BARTS_2023_12_PATH = AnyPath(
    PIPELINE_RAW_DATA_PATH, 
    "secondary_care",
    "DSA__BartsHealth_NHS_Trust", 
    "2023_12_ResearchDatasetv1.6",
)

# Input file
BARTS_2023_12_PATHOLOGY_FILE_RAW_PATH = AnyPath(
    BARTS_2023_12_PATH,
    "GH_Pathology__20231218.ascii.nohisto.redacted2.tab"
)

# Output file
# No relevant here


# In[ ]:


NUM_DELIM_SPLITTER_PATH = AnyPath(
    PIPELINE_HELPERS_LOCATION,
    "num_delims_splitter.sh"
)

subprocess.run(
    [NUM_DELIM_SPLITTER_PATH, BARTS_2023_12_PATHOLOGY_FILE_RAW_PATH],
    capture_output=True,
    text=True,
)


# In[ ]:


get_ipython().run_cell_magic('time', '', 'provenance_key = "2023_12_Barts_path"\n(\npl.scan_csv(\n    AnyPath(\n        BARTS_2023_12_PATH,\n        "GH_Pathology__20231218.ascii.nohisto.redacted2_tab16.tab"\n        ),\n    separator="\\t",\n    infer_schema=False,\n    )\n    .filter(\n        ~pl.col("ResultTxt").str.contains(r"[a-zA-Z]"),\n        pl.col("ResultTxt").ne("**"),\n        pl.col("ResultTxt").ne("***"),\n        pl.col("ResultTxt").ne("****"),\n        pl.col("ResultTxt").ne("*****"),\n        pl.col("ResultTxt").ne("*"),\n        pl.col("ResultTxt").ne("* -"),\n        pl.col("ResultTxt").ne("-"),\n        pl.col("ResultTxt").ne("--"),\n        pl.col("ResultTxt").ne("- -"),\n        pl.col("ResultTxt").ne("-  -"),\n        pl.col("ResultTxt").ne("++++"),\n        pl.col("ResultTxt").ne("#"),\n        pl.col("ResultTxt").ne("`"),\n        pl.col("ResultTxt").ne("....."),\n        pl.col("ResultTxt").ne(",."),\n        pl.col("ResultTxt").ne("n/r"),\n        pl.col("ResultTxt").ne("na"),\n        pl.col("ResultTxt").ne("n/a"),\n        pl.col("ResultTxt").ne("NA"),\n        pl.col("ResultTxt").ne("?"),\n        pl.col("ResultTxt").ne(","),\n        pl.col("ResultTxt").ne(":"),\n        pl.col("ResultTxt").ne("]"),\n        pl.col("ResultTxt").ne("c"),\n        pl.col("ResultTxt").ne("MK"),\n        pl.col("ResultTxt").ne("B"),\n        pl.col("ResultTxt").ne("P"),\n        pl.col("ResultTxt").ne("1a"),\n        pl.col("ResultTxt").ne("1b"),\n        pl.col("ResultTxt").ne("3a"),\n        pl.col("ResultTxt").ne("3b"),\n        pl.col("ResultTxt").ne("3-"),\n        pl.col("ResultTxt").ne("64-"),\n        pl.col("ResultTxt").ne("B2A2"),\n        pl.col("ResultTxt").ne("B3A2"),\n        pl.col("ResultTxt").ne("FM"),\n        pl.col("ResultTxt").ne("UNS"),\n        pl.col("ResultTxt").ne("@unb"),\n        pl.col("ResultTxt").ne("@und"),\n        pl.col("ResultTxt").ne("None"),\n        pl.col("ResultTxt").ne("2-5"),\n        pl.col("ResultTxt").ne("1:8"),\n        pl.col("ResultTxt").ne("1:16"),\n        pl.col("ResultTxt").ne("1:32"),\n        pl.col("ResultTxt").ne("4o"),\n        pl.col("ResultTxt").ne("*40"),\n        pl.col("ResultTxt").ne("body"),\n        pl.col("ResultTxt").ne("Body"),\n        pl.col("ResultTxt").ne("24hr"),\n        pl.col("ResultTxt").ne("24 hrs"),\n        pl.col("ResultTxt").ne("24HR"),\n        pl.col("ResultTxt").ne("KNIB"),\n        pl.col("ResultTxt").ne("64 -"),\n        pl.col("ResultTxt").ne("70)"),\n        pl.col("ResultTxt").ne("(70)"),\n        pl.col("ResultTxt").ne("(66"),\n        pl.col("ResultTxt").ne("other"),\n        pl.col("ResultTxt").ne("Clear"),\n        pl.col("ResultTxt").ne("rerun"),\n        pl.col("ResultTxt").ne("Venous"),\n        pl.col("ResultTxt").ne("{REPEAT}"),\n        pl.col("ResultTxt").ne("deleted"),\n        pl.col("ResultTxt").ne("DELETED"),\n        pl.col("ResultTxt").ne("09:00"),\n        pl.col("ResultTxt").ne("10:17"),\n        pl.col("ResultTxt").ne("10:38"),\n        pl.col("ResultTxt").ne("11:30"),\n        pl.col("ResultTxt").ne("16:00"),\n        pl.col("ResultTxt").ne("18:00"),\n        pl.col("ResultTxt").ne("21:00"),\n        pl.col("ResultTxt").ne("23:59"),\n        pl.col("ResultTxt").ne("day 1"),\n        pl.col("ResultTxt").ne("day 2"),\n        pl.col("ResultTxt").ne("Day 2"),\n        pl.col("ResultTxt").ne("DAY 2"),\n        pl.col("ResultTxt").ne("day 3"),\n        pl.col("ResultTxt").ne("Day 4"),\n        pl.col("ResultTxt").ne("day 7"),\n        pl.col("ResultTxt").ne("Day 8"),\n        pl.col("ResultTxt").ne("Day 10"),\n        pl.col("ResultTxt").ne("day 17"),\n        pl.col("ResultTxt").ne("day 21"),\n        pl.col("ResultTxt").ne("Day 21"), \n        pl.col("ResultTxt").ne("DAY 21"), \n        pl.col("ResultTxt").ne("0 min"),\n        pl.col("ResultTxt").ne("30 min"),\n        pl.col("ResultTxt").ne("60 min"),\n        pl.col("ResultTxt").ne("4 hrs"),\n        pl.col("ResultTxt").ne("7.5 hrs"),\n        pl.col("ResultTxt").ne("44285*"),\n        pl.col("ResultTxt").ne("20753*"),\n        pl.col("ResultTxt").ne("124 -"),\n        pl.col("ResultTxt").ne("LCMSMS"),\n        pl.col("ResultTxt").ne("1.01 26"),\n        pl.col("ResultTxt").ne("1.20 15"),\n        pl.col("ResultTxt").ne("0.99 10"),\n        pl.col("ResultTxt").ne("0.99 11"),\n        pl.col("ResultTxt").ne("0.99 12"),\n        pl.col("ResultTxt").ne("1.00 10"),\n        pl.col("ResultTxt").ne("2.41 32"),\n        pl.col("ResultTxt").ne("0.95 14"),\n        pl.col("ResultTxt").ne("0.95 17"),\n        pl.col("ResultTxt").ne("1.05 9"), \n        pl.col("ResultTxt").ne("0.94 26"),\n        pl.col("ResultTxt").ne("Add on"),\n        pl.col("ResultTxt").ne("clumps"),\n        pl.col("ResultTxt").ne("Clumped"),\n        pl.col("ResultTxt").ne("*Clumped"),\n        pl.col("ResultTxt").ne("clumped"),\n        pl.col("ResultTxt").ne("Clumpled"),\n        pl.col("ResultTxt").ne("no clot"),\n        pl.col("ResultTxt").ne("No clot"),\n        pl.col("ResultTxt").ne("NO CLOT"),\n        pl.col("ResultTxt").ne("Pending"),\n        pl.col("ResultTxt").ne("IgM only"),\n        pl.col("ResultTxt").ne(">1/640"),\n        pl.col("ResultTxt").ne("1/640"),\n        pl.col("ResultTxt").ne("1/160"),\n        pl.col("ResultTxt").ne("Cloudy"),\n        pl.col("ResultTxt").ne("Pleural"),\n        pl.col("ResultTxt").ne("ramdom"),\n        pl.col("ResultTxt").ne("random"),\n        pl.col("ResultTxt").ne("Random"),\n        pl.col("ResultTxt").ne("RANDOM"),\n        pl.col("ResultTxt").ne("Reject"),\n        pl.col("ResultTxt").ne("normal"),\n        pl.col("ResultTxt").ne("Normal"),\n        pl.col("ResultTxt").ne("NORMAL"),\n        pl.col("ResultTxt").ne("invalid"),\n        pl.col("ResultTxt").ne("Note Hb"),\n        pl.col("ResultTxt").ne("reduced"),\n        pl.col("ResultTxt").ne("Reduced"),\n        pl.col("ResultTxt").ne("Unknown"),\n        pl.col("ResultTxt").ne("DR req"), \n        pl.col("ResultTxt").ne("?on GCSF"),\n        pl.col("ResultTxt").ne("MDS/MPN"),\n        pl.col("ResultTxt").ne("Arterial"),\n        pl.col("ResultTxt").ne("CAPASCIN"),\n        pl.col("ResultTxt").ne("Detected"),\n        pl.col("ResultTxt").ne("negative"),\n        pl.col("ResultTxt").ne("Negative"),\n        pl.col("ResultTxt").ne("NEGATIVE"),\n        pl.col("ResultTxt").ne("Neagtive"),\n        pl.col("ResultTxt").ne("positive"),\n        pl.col("ResultTxt").ne("Positive"),\n        pl.col("ResultTxt").ne("POSITIVE"),\n        pl.col("ResultTxt").ne("No clot."),\n        pl.col("ResultTxt").ne("Obscured"),\n        pl.col("ResultTxt").ne("See ADAL"),\n        pl.col("ResultTxt").ne("Speckled"),\n        pl.col("ResultTxt").ne("Stained"),\n        pl.col("ResultTxt").ne("Pendings"),\n        pl.col("ResultTxt").ne("Rejected"),\n        pl.col("ResultTxt").ne("33 hours"), \n        pl.col("ResultTxt").ne("09S00088662 Read code 43X4 Read code 43BA"),\n        ~pl.col("ResultTxt").str.contains("^100-149 mIU/ml Low Level Antibody detected Low level VZV IgG detected For immunocompromised patients recently exposed to VZV"),\n        ~pl.col("ResultTxt").str.contains("(?i)unsuitable"), # rule out e.g. ["6ml EDTA sample tube unsuitable for FBC or ESR analyser. Please send 4ml EDTA tube."]\n        ~pl.col("ResultTxt").str.contains("(?i)not been accepted"), # rule out e.g. ["6ml EDTA sample tube unsuitable for FBC or ESR analyser. Please send 4ml EDTA tube."]\n        ~pl.col("ResultTxt").str.contains("^\\d{2}[A-Z]\\d{8}"), # rule out e.g. ["09S00053956 ..."]\n        ~pl.col("ResultTxt").str.contains("^\\d+.*?\\*\\s\\*\\s"), # rule out e.g. ["14 + 4* * likely to be an over-estimation due to the polyclonal background of gamma globulins."]\n        ~pl.col("ResultTxt").str.contains(r"^\\d{2}/\\d{2}/\\d{2,4}.? \\d{2}:\\d{2}$") #"14/07/2011, 16:51"\n    )\n    \n    .with_columns(\n        pl.col("PseudoNHS_2023_11_08").alias("pseudo_nhs_number"),\n        pl.col("ReportDate").str.to_date(format="%Y-%m-%d %H:%M").alias("test_date"), ### ?REPORTDate\n        pl.col("TestDesc").alias("original_term"),\n        ## conversion from `str` to `f64` failed in column \'ResultTxt\' for 810 out of 32897 values: [">90", ">90", … "<1"]\n        pl.col("ResultTxt")        \n            .str.strip_prefix("<")\n            .str.strip_prefix(">")\n            .str.strip_prefix("+-") ## present in 2023_12\n            .str.strip_prefix("+/-") ## present in 2023_12\n            .str.replace(r"^\\{(.*?)\\}$","$1")\n            .str.replace(r"^\\((.*?)\\)$","$1")\n            .str.replace(r"^\\*\\{(.*?)\\}$","$1")\n            .str.strip_prefix(" ")\n            .cast(pl.Float64, strict=True)\n            .alias("result"),\n        pl.col("ResultUnit").alias("result_value_units"),\n        provenance=pl.lit(provenance_key, pl.Enum(ALL_PROVENANCE_OPTIONS)),\n        source=pl.lit("secondary_care", pl.Enum(ALL_SOURCE_OPTIONS)),\n        \n    )\n\n    .filter(\n        pl.col("test_date").is_not_null(),\n    )\n    .with_columns(\n        HASH_COLUMN\n    )\n    .unique("hash")\n    \n    .select(\n        TARGET_OUTPUT_COLUMNS_WITH_HASH\n    )\n    \n    .sink_ipc(\n        AnyPath(\n            SECONDARY_ARROW_PATH,\n            f"{provenance_key}.arrow"\n        ),\n    )\n)\n')


# #### `Barts 2024_09 - Sept 2024 - Pathology`

# ##### Note about `Barts 2024_09 - Sept 2024 - Pathology`
# 
# Note that the files is v2 of the redacted `RDE_Pathology` file, hence `redacted2.tab`.
# 
# File `RDE_Pathology.ascii.nohisto.redacted2.csv` contains multiple `"""` (triple double-quotes) and indeed many `"` (single double-quite) or `""` (double double-quote).
# 
# As the file is in fact a tab delinited file (_despite the suffix_), they are not necessary and affect the parsing by polars.
# 
# We therefore parse a copy of the file with **all** double-quotes removed.  This is generated using the following shell command:
# 
# `tr -d '"' < RDE_Pathology.ascii.nohisto.redacted2.csv > RDE_Pathology.ascii.nohisto.redacted2.no_double_quotes.csv`

# #### HARD-CODED PRE-PROCESSING: BARTS_2024_09 (REMOVAL OF DOUBLE-QUOTES)

# In[ ]:


BARTS_2024_09_PATH = AnyPath(
    PIPELINE_RAW_DATA_PATH, 
    "secondary_care",
    "DSA__BartsHealth_NHS_Trust",
    "2024_09_ResearchDataset"
)

# Input file
BARTS_2024_09_PATHOLOGY_FILE_RAW_PATH = AnyPath(
    BARTS_2024_09_PATH,
    "RDE_Pathology.ascii.nohisto.redacted2.csv"
)

# Output file
BARTS_2024_09_PATHOLOGY_FILE_CORRECTED_PATH = AnyPath(
    BARTS_2024_09_PATH,
    "RDE_Pathology.ascii.nohisto.redacted2.no_double_quotes.csv"
)


# In[ ]:


barts_2024_09_preprocessing_command = (
    f"""tr -d '"' < """
    f"""\"{BARTS_2024_09_PATHOLOGY_FILE_RAW_PATH}\" > """
    f"""\"{BARTS_2024_09_PATHOLOGY_FILE_CORRECTED_PATH}\""""
)


# In[ ]:


subprocess.run(
    barts_2024_09_preprocessing_command,
    shell=True,
    check=True,
    capture_output=True,
    text=True
)


# In[ ]:


get_ipython().run_cell_magic('time', '', 'provenance_key="2024_09_Barts_path"\n\n(\npl.scan_csv(\n    AnyPath(\n        BARTS_2024_09_PATHOLOGY_FILE_CORRECTED_PATH\n        ),\n    separator="\\t",\n    infer_schema=False,\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("ResultTxt").str.contains("[a-zA-Z]"),\n        pl.col("ResultTxt").ne("-"),\n        label="Lots of [a-zA-Z] values in `result`"\n    )\n    .TRE\n    .filter_with_logging( # ". . . . .", "(66", … "."\n        pl.col("ResultTxt").ne("**"),\n        pl.col("ResultTxt").ne("***"),\n        pl.col("ResultTxt").ne("****"),\n        pl.col("ResultTxt").ne("*****"),\n        pl.col("ResultTxt").ne("*"),\n        pl.col("ResultTxt").ne("* -"),\n        pl.col("ResultTxt").ne("-"),\n        pl.col("ResultTxt").ne("--"),\n        pl.col("ResultTxt").ne("- -"),\n        pl.col("ResultTxt").ne("-  -"),\n        pl.col("ResultTxt").ne("- ."),\n        pl.col("ResultTxt").ne(". ."),\n        pl.col("ResultTxt").ne(". . ."),\n        pl.col("ResultTxt").ne("----"),\n        pl.col("ResultTxt").ne("+"),\n        pl.col("ResultTxt").ne("+++"),\n        pl.col("ResultTxt").ne("++++"),\n        pl.col("ResultTxt").ne("#"),\n        pl.col("ResultTxt").ne("`"),\n        pl.col("ResultTxt").ne("-."),\n        pl.col("ResultTxt").ne("....."),\n        pl.col("ResultTxt").ne(". . . . ."),\n        pl.col("ResultTxt").ne(",."),\n        pl.col("ResultTxt").ne("#"),\n        pl.col("ResultTxt").ne("`"),\n        pl.col("ResultTxt").ne("."),\n        pl.col("ResultTxt").ne("."),\n        pl.col("ResultTxt").ne("....."),\n        pl.col("ResultTxt").ne(",."),\n        pl.col("ResultTxt").ne("?"),\n        pl.col("ResultTxt").ne(","),\n        pl.col("ResultTxt").ne(":"),\n        pl.col("ResultTxt").ne("]"),\n        pl.col("ResultTxt").ne("{.}"),\n        label="Just symbols and space in `result`"\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col(\'ResultTxt\').is_in(\n            [\n                ">1/640",\n                "28.8 28.8",\n                "28.3 28.3",\n                "1:8",\n                "2+48",\n                "2+0",\n                "{4}",\n                "1:32",\n                "{88}",\n                "3-",\n                "{93}",\n                "(66",\n                "1:32",\n                "1:16",\n                "106 - - - - - -"\n            ]\n        ),\n        ~pl.col("ResultTxt").str.contains("^\\d+(\\.\\d+)? \\d+(\\.\\d+)?$"),\n        label="Number-like, with extra spaces or symbols inside"\n    )\n    .TRE # "22.01.15; 1800", "?45.5", … "- ."\n    .filter_with_logging(\n        ~pl.col("ResultTxt").str.contains("^\\d{2}:\\d{2}$"),\n        label="Time-like (e.g. 09:59)"\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("ResultTxt").str.contains("^\\d*\\s?-$"),\n        label="digits Ending in `-` or \' -\'"\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("ResultTxt").str.contains("^\\$|^\\*|^\\?"),\n        label="Starting with `$` or \'*\' or \'?\'"\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("ResultTxt").str.contains("\\*$"),\n        label="Ending with \'*\'"\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("ResultTxt").str.contains("\\d*\\+$"),\n        label="digits ending with \'+\'"\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("ResultTxt").str.contains("\\d+(\\.\\d+)?%$"),\n        label="digits ending with \'%\'"\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("ResultTxt").str.contains("/.*/"),\n        ~pl.col("ResultTxt").str.contains("\\d{2}\\.\\d{2}\\.\\d{2}; \\d{4}"), # "22.01.15; 1800"\n        label="Date-, time-,  or datetime-like in `result`"\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("ResultTxt").str.contains("/"),\n        label="Fraction-like in `result`"\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("ResultTxt").str.contains("\\d+-\\d+"),\n        label="Integer range in `result` (e.g. \'92-99\')"\n    )\n    .TRE\n    .filter_with_logging(\n        pl.col("ReportDate").str.contains("\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}"),\n        label="ReportDate in valid format"\n    )\n    .with_columns(\n        pl.col("PseudoNHS_2024-07-10").alias("pseudo_nhs_number"),\n        pl.col("ReportDate").str.to_date(format="%Y-%m-%d %H:%M", strict=True).alias("test_date"),\n        pl.col("TestDesc").alias("original_term"),\n        pl.col("ResultTxt")        \n            .str.strip_prefix("<")\n            .str.strip_prefix(">")\n            .str.strip_prefix("+-") ## present in 2023_12\n            .str.strip_prefix("+/-") ## present in 2023_12\n            .str.strip_suffix("cm")\n            .str.strip_prefix("(")\n            .str.strip_suffix(")")\n            .str.strip_prefix(" ")\n            .cast(pl.Float64, strict=True)\n            .alias("result"),\n        pl.col("ResultUnit").alias("result_value_units"),\n        provenance=pl.lit(provenance_key, pl.Enum(ALL_PROVENANCE_OPTIONS)),\n        source=pl.lit("secondary_care", pl.Enum(ALL_SOURCE_OPTIONS)),\n\n    )\n    .with_columns(\n        HASH_COLUMN\n    )\n    .unique("hash")\n\n    .select(\n        TARGET_OUTPUT_COLUMNS_WITH_HASH\n    )\n    .sink_ipc(\n        AnyPath(\n            SECONDARY_ARROW_PATH,\n            f"{provenance_key}.arrow"\n        )\n    )\n)\n')


# ### Combining Barts pathology data

# In[ ]:


Barts_path_arrow_files = [
    '2021_04_Barts_path.arrow',
    '2023_05_Barts_path.arrow',
    '2024_09_Barts_path.arrow',
    '2022_03_Barts_path.arrow',
    '2023_12_Barts_path.arrow'
]


# In[ ]:


get_ipython().run_cell_magic('time', '', '(\n    pl.scan_ipc(\n        AnyPath(\n            SECONDARY_ARROW_PATH,\n            "*_Barts_path.arrow"\n        )\n    )\n\n    .with_columns(\n        HASH_COLUMN\n    )\n    .unique("hash")\n    .sink_ipc(\n        AnyPath(\n            SECONDARY_ARROW_PATH,\n            f"{yr}_{mon}_Barts_path_combined.arrow"\n        )\n    )\n)\n \n')


# ### Barts Measurement Data

#  #### `2023_05_Barts_measurements`

# In[ ]:


# Input file
BARTS_2023_05_MEASUREMENTS_FILE_RAW_PATH = AnyPath(
    BARTS_2023_05_PATH,
    "GandH_Measurements_202305151304.ascii.redacted.tab"
)


# In[ ]:


get_ipython().run_cell_magic('time', '', 'provenance_key="2023_05_Barts_measurements"\n\n(\n    pl.scan_csv(\n        # PseudoNHS_2023_04_24\tSystemLookup\tClinicalSignificanceDate\tEventResult\tUnitsCode\t\n        # UnitsDesc\tNormalCode\tNormalDesc\tLowValue\tHighValue\tEventText\tEventType\tEventParent\n        BARTS_2023_05_MEASUREMENTS_FILE_RAW_PATH,\n        separator="\\t",\n        infer_schema=False,\n    )\n\n    .filter(\n        ~pl.col("EventResult").str.contains(r"^\\..*?"), # only 8 rows, strip out values in [".", ".", ".", ".", ".", ".", ".", ".2.2"]\n    )\n    .with_columns(\n        pl.col("PseudoNHS_2023_04_24").alias("pseudo_nhs_number"),\n        ## conversion from `str` to `date` failed in column \'ClinicalSignificanceDate\' for 17288 out of 17288 values: ["Apr 11 2022  5:12AM", "Apr 11 2022  5:12AM", … "Jan 31 2019 10:25AM"]\n        pl.col("ClinicalSignificanceDate").str.to_date(format="%b %d %Y %I:%M%p").alias("test_date"), # %I for 12-hour clock\n        pl.col("EventType").alias("original_term"),\n        ## conversion from `str` to `f64` failed in column \'ResultTxt\' for 810 out of 32897 values: [">90", ">90", … "<1"]\n        pl.col("EventResult") \n        .str.strip_suffix("cm")\n        .str.replace("3\\.6\\.1","36.1") # this should be a degrees celcius value for `"SN - Preop - CTm - Patient Tem…`\n        .cast(pl.Float64, strict=True)\n        .alias("result"),\n        pl.col("UnitsDesc").alias("result_value_units"),\n        provenance=pl.lit(provenance_key, pl.Enum(ALL_PROVENANCE_OPTIONS)),\n        source=pl.lit("secondary_care", pl.Enum(ALL_SOURCE_OPTIONS)),\n    )\n\n    .filter(\n        pl.col("test_date").is_not_null(),\n    )\n    .with_columns(\n        HASH_COLUMN\n    )\n    .unique("hash")\n\n    .select(\n        TARGET_OUTPUT_COLUMNS_WITH_HASH\n    )\n\n    .sink_ipc(\n        AnyPath(\n            SECONDARY_ARROW_PATH,\n            f"{provenance_key}.arrow")\n        )\n)\n')


#  #### `2023_12_Barts_measurements` 
#  

# ##### Pre-processing notes
# 
# Chosen option here is to subset master files into number of tabs files and use the most common one.
# 
# However, for this file, you also need to remove **all** double-quotes, not just unmatched ones.  If you do not, the `select` statement in the cell cause the number of row to fall dramatically for reasons not entirely clear (possibly internal to polars).
# 
# ##### Necessary pre-processing steps
# 
# 1. `tr -d '"' < [input_file_with_double_quotes] >  [output_file_with_no_double_quotes]`
# 2. `num_delim_splitter.sh`
# 
# ##### Outputs
# 
# ```
#   15047135 GandH_Measurements__20240423.ascii.redacted2.no_double_quotes_tab13.tab
#          8 GandH_Measurements__20240423.ascii.redacted2.no_double_quotes_tab14.tab
#         65 GandH_Measurements__20240423.ascii.redacted2.no_double_quotes_tab15.tab
#          3 GandH_Measurements__20240423.ascii.redacted2.no_double_quotes_tab16.tab
#          8 GandH_Measurements__20240423.ascii.redacted2.no_double_quotes_tab17.tab
#          2 GandH_Measurements__20240423.ascii.redacted2.no_double_quotes_tab18.tab
#          4 GandH_Measurements__20240423.ascii.redacted2.no_double_quotes_tab19.tab
#          1 GandH_Measurements__20240423.ascii.redacted2.no_double_quotes_tab21.tab
#   15047226 total
# 
# ```
# 
# There are a trivial number of non `_tab13` rows.

# #### HARD-CODED PRE-PROCESSING: BARTS_2023_12 (MEASUREMENTS) (REMOVAL OF DOUBLE-QUOTES + NUM_DELIM_SPLITTER)

# In[ ]:


# Input file
BARTS_2023_12_MEASUREMENTS_FILE_RAW_PATH = AnyPath(
    BARTS_2023_12_PATH,
    "GandH_Measurements__20240423.ascii.redacted2.tab"
)

# Output file 1 (tr processed)
BARTS_2023_12_MEASUREMENTS_FILE_CORRECTED_PATH = AnyPath(
    BARTS_2023_12_PATH,
    "GandH_Measurements__20240423.ascii.redacted2.no_double_quotes.tab"
)


# In[ ]:


barts_2023_12_measurements_tr_preprocessing_command = (
    f"""tr -d '"' < """
    f"""\"{BARTS_2023_12_MEASUREMENTS_FILE_RAW_PATH}\" > """
    f"""\"{BARTS_2023_12_MEASUREMENTS_FILE_CORRECTED_PATH}\""""
)


# In[ ]:


subprocess.run(
    barts_2023_12_measurements_tr_preprocessing_command,
    shell=True,
    check=True,
    capture_output=True,
    text=True
)


# In[ ]:


subprocess.run(
    [NUM_DELIM_SPLITTER_PATH, BARTS_2023_12_MEASUREMENTS_FILE_CORRECTED_PATH],
    capture_output=True,
    text=True,
)


# In[ ]:


get_ipython().run_cell_magic('time', '', 'provenance_key="2023_12_Barts_measurements"\n\n(\n    pl.scan_csv(\n    # PseudoNHS_2023_04_24\tSystemLookup\tClinicalSignificanceDate\tEventResult\tUnitsCode\t\n    # UnitsDesc\tNormalCode\tNormalDesc\tLowValue\tHighValue\tEventText\tEventType\tEventParent\n        AnyPath(\n            BARTS_2023_12_PATH,\n            "GandH_Measurements__20240423.ascii.redacted2.no_double_quotes_tab13.tab"\n            ),\n        separator="\\t",\n        infer_schema=False,\n        )\n         .with_columns(\n            pl.col("EventResult").str.strip_prefix(" ")\n        )\n    .TRE\n    .filter_with_logging( \n        ~pl.col("EventResult").str.contains("\\d:\\d{16}:\\d\\.000000:\\d{1,3}:0"),\n        pl.col("EventResult").ne("06.01.2010"), \n        pl.col("EventResult").ne("10:00"),\n        pl.col("EventResult").ne("10%"),\n        pl.col("EventResult").ne("."), # special case for "2023_12_Barts_measurements" and "2024_09_Barts_measurements". rules out "."\n        ~pl.col("EventResult").str.contains(r"[\\+\\)a-zA-Z/\\s]"), #rule out ["23/11", ")9", "text…"] . This enough to attain strict casting to pl.Float64\n        label=\'Remove various non-numeric/weird results\'\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("EventResult").str.contains(r"\\d{1,2}\\.\\d{1,2}\\.\\d{2,4}"),\n        ~pl.col("EventResult").str.contains(r"\\d{1,2}:\\d{2}"),\n        label="Date-like or Time-like string in `result`"\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("EventResult").str.contains(r"^\\d+(\\.\\d+)?%$"),\n        label="Number ends with \'%\' in `result`"\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("EventResult").str.contains(r"^\\d+(\\.\\d+)?`$"),\n        label="Number ends with \'`\' in `result`"\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("EventResult").str.contains(r"`"),\n        label="Contains \'`\' in `result` (e.g. \'1`437\')"\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("EventResult").str.contains(r"="),\n        label="Contains \'=\' in `result`"\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("EventResult").str.contains(r"^\\d+:\\d+$"),\n        label="Contains single \':\' in `result`"\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("EventResult").str.contains(r"^\\.+$"),\n        label="Contains just \'.\'s in `result`"\n    ) \n    .TRE\n    .filter_with_logging(\n        ~pl.col("EventResult").str.contains(r"\\.:"),\n        label="Contains  \'.:\' in `result`"\n    ) \n    .TRE\n    .filter_with_logging(\n        ~pl.col("EventResult").str.contains(r":!"),\n        label="Contains  \':!\' in `result`"\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("EventResult").str.contains(r":"),\n        label="Contains  \':\' in `result`"\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("EventResult").str.contains(r";"),\n        label="Contains  \';\' in `result`"\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("EventResult").str.contains(r"_"),\n        label="Contains  \'_\' in `result`"\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("EventResult").str.contains(r"#"),\n        label="Contains  \'#\' in `result`"\n    )\n    .TRE\n    .filter_with_logging(\n        pl.col("EventResult").ne("?"),\n        label="Literal \'?\' in `result`"\n    ) \n    .TRE\n    .filter_with_logging(\n        ~pl.col("EventResult").str.contains("\\*"),\n        label="Contains \'*\' in `result`"\n    ) \n    .TRE\n    .filter_with_logging(\n        ~pl.col("EventResult").str.contains("\\..*\\."),\n        label="Contains more than \'.\' in `result`"\n    ) \n    .TRE\n    .filter_with_logging(\n        ~pl.col("EventResult").str.contains("\'"),\n        label="Contains \'\\\'\' in `result`"\n    ) \n    .TRE\n    .filter_with_logging(\n        ~pl.col("EventResult").str.ends_with("&"),\n        label="Ends with \'&\' in `result`"\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("EventResult").str.contains("^-+$"),\n        label="Contains only one (or more) \'-\'s in `result`"\n    ) \n    .TRE\n    .filter_with_logging(\n        ~pl.col("EventResult").str.contains("\\d-+\\d"),\n        label="Contains one or more dashes between digits, e.g. 14-40, in `result`"\n        \n    ) \n    .TRE\n    .filter_with_logging(\n        ~pl.col("EventResult").str.ends_with("-"),\n        label="Ends with \'-\' in `result`"\n    )\n    .TRE\n    .filter_with_logging(\n        ~pl.col("EventResult").str.contains(r\'\\\\\'),\n        label="Contains \'\\\\\' (backslash) in `result`"\n    ) \n    .with_columns(\n        pl.col("PseudoNHS_2023_11_08").alias("pseudo_nhs_number"),\n        pl.col("ClinicalSignificanceDate").str.to_date(format="%b %d %Y %I:%M%p").alias("test_date"), # %I for 12-hour clock\n        pl.col("EventType").alias("original_term"),\n        pl.col("EventResult") \n            .str.strip_prefix(">")\n            .cast(pl.Float64, strict=True)\n            .alias("result"),\n        pl.col("UnitsDesc").alias("result_value_units"),\n        provenance=pl.lit(provenance_key, pl.Enum(ALL_PROVENANCE_OPTIONS)),\n        source=pl.lit("secondary_care", pl.Enum(ALL_SOURCE_OPTIONS)),\n        \n    )\n   \n    .with_columns(\n        HASH_COLUMN\n    )\n    .unique("hash")\n    \n    .select(\n        TARGET_OUTPUT_COLUMNS_WITH_HASH\n    )\n\n    .sink_ipc(\n        AnyPath(\n            SECONDARY_ARROW_PATH,\n            f"{provenance_key}.arrow")\n    )\n)\n')


#  #### `2024_09_Barts_measurements` 

# ##### Pre-processing
# 
# We split by number of tabs (and therefore don't have to use the debatable `truncate_ragged_lines=True`.
# For 'belts-and'braces' we checked for double-quotes in `RDE_Measurements.ascii.redacted.tab` but there are none.
# 
# ```
#    7564328 RDE_Measurements.ascii.redacted.tab
#    7564309 RDE_Measurements.ascii.redacted_tab13.tab
#         15 RDE_Measurements.ascii.redacted_tab14.tab
#          3 RDE_Measurements.ascii.redacted_tab15.tab
#          1 RDE_Measurements.ascii.redacted_tab16.tab
#   22692965 total
# ```

# #### HARD-CODED PRE-PROCESSING: BARTS_2024_09 (MEASUREMENTS)  (NUM_DELIM_SPLITTER)

# In[ ]:


# Input file
BARTS_2024_09_MEASUREMENTS_FILE_RAW_PATH = AnyPath(
    BARTS_2024_09_PATH,
    "RDE_Measurements.ascii.redacted2.tab"
)


# In[ ]:


subprocess.run(
    [NUM_DELIM_SPLITTER_PATH, BARTS_2024_09_MEASUREMENTS_FILE_RAW_PATH],
    capture_output=True,
    text=True,
)


# In[ ]:


get_ipython().run_cell_magic('time', '', 'provenance_key="2024_09_Barts_measurements"\n\n(\n    pl.scan_csv(\n    # PseudoNHS_2024-07-10\tSystemLookup\tClinicalSignificanceDate\tResultNumeric\tEventResult\tUnitsCode\t\n    # UnitsDesc\tNormalCode\tNormalDesc\tLowValue\tHighValue\tEventText\tEventType\tEventParent\n        AnyPath(\n            BARTS_2024_09_PATH,\n            "RDE_Measurements.ascii.redacted2_tab13.tab"\n            ),\n        separator="\\t",\n        infer_schema=False,\n    )\n    \n    .filter( \n        pl.col("UnitsDesc").ne("0"), # special case for "2023_12_Barts_measurements" and "2024_09_Barts_measurements".\n        pl.col("EventResult").ne("."), # special case for "2023_12_Barts_measurements" and "2024_09_Barts_measurements". rules out "."\n        pl.col("EventResult").ne(".2.2"), # rules out 1 row\n        ~pl.col("EventResult").str.contains(r"[\\)a-zA-Z/\\s-]") #rule out ["23/11", ")9", "text…"] . This enough to attain strict casting to pl.Float64\n    )\n    .with_columns(\n        pl.col("PseudoNHS_2024-07-10").alias("pseudo_nhs_number"),\n        ## conversion from `str` to `date` failed in column \'ClinicalSignificanceDate\' for 17288 out of 17288 values: ["Apr 11 2022  5:12AM", "Apr 11 2022  5:12AM", … "Jan 31 2019 10:25AM"]\n        pl.col("ClinicalSignificanceDate").str.to_date(format="%b %d %Y %I:%M%p").alias("test_date"), # %I for 12-hour clock\n        pl.col("EventType").alias("original_term"),\n         ## conversion from `str` to `f64` failed in column \'ResultTxt\' for 810 out of 32897 values: [">90", ">90", … "<1"]\n        pl.when(pl.col("EventType").eq("Child\'s Birth Weight (g)"))\n            .then(pl.col("EventResult").str.replace_all(",", ""))\n            .otherwise(pl.col("EventResult"))\n        .str.replace("3\\.6\\.1","36.1") # this should be a degrees celcius value for `"SN - Preop - CTm - Patient Tem…`\n        .cast(pl.Float64, strict=True)\n        .alias("result"),\n        pl.col("UnitsDesc").alias("result_value_units"),\n        provenance=pl.lit(provenance_key, pl.Enum(ALL_PROVENANCE_OPTIONS)),\n        source=pl.lit("secondary_care", pl.Enum(ALL_SOURCE_OPTIONS)),\n\n    )\n\n    .with_columns(\n        HASH_COLUMN\n    )\n    .unique("hash")\n\n    .select(\n        TARGET_OUTPUT_COLUMNS_WITH_HASH\n    )\n\n    .sink_ipc(\n        AnyPath(\n            SECONDARY_ARROW_PATH,\n            f"{provenance_key}.arrow")\n    )\n\n)\n')


# ### Combine Barts Measurement Data
# 

# In[ ]:


get_ipython().run_cell_magic('time', '', '(\n    pl.scan_ipc(\n        AnyPath(\n            SECONDARY_ARROW_PATH,\n            "20*_Barts_measurements.arrow"\n        )\n    )\n    .with_columns(\n         HASH_COLUMN\n    )\n    .unique("hash")\n\n    .sink_ipc(\n        AnyPath(\n            SECONDARY_ARROW_PATH,\n            f"{yr}_{mon}_Barts_measurements_combined.arrow"\n        )\n    )\n)\n')


# ## Combine all secondary care files

# In[ ]:


get_ipython().run_cell_magic('time', '', '(\n    pl.scan_ipc(\n        AnyPath(\n            SECONDARY_ARROW_PATH,\n            "*_combined.arrow"\n        )\n    )\n    ## The re-hashing is added to protect the script from polars version changes as hashing consistency\n    ## is not guaranteed between polars version.  If no polars update, one could consider using the \n    ## pre-existing hashes calculated per secondary_care arrow file.\n    .with_columns(\n         HASH_COLUMN\n    )\n    .unique("hash")\n\n    .sink_ipc(\n        AnyPath(\n            COMBINED_DATASETS_ARROW_PATH,\n            f"{yr}_{mon}_Combined_secondary_care.arrow"\n        )\n    )\n)\n')


# ## Combine all primary and secondary datasets
# 
# 
# We have to combine here to then allow the generic handling of unitless data

# In[ ]:


get_ipython().run_cell_magic('time', '', 'combined_primary_and_secondary = (\n    pl.scan_ipc(\n        AnyPath(\n            COMBINED_DATASETS_ARROW_PATH,\n            f"{yr}_{mon}_Combined_*.arrow"\n        )\n    )\n    ## The re-hashing is added to protect the script from polars version changes as hashing consistency\n    ## is no guaranteed between polars version.  If no polars update, one could consider using the \n    ## pre-existing hashes calculated per secondary_care arrow file.    \n    .with_columns(\n         HASH_COLUMN\n    )\n    .unique("hash")\n)\n')


# # Handle unitless data
# 
# The simplest way to do this is to artifically add a unit to rows with a presumed unit.
# 
# You could consider adding some heuristics to narrow down rows to unit-assign.

# ## POCT ketones
# 
# Justifications for these decision can be found in `trait_augmenter` and `ketone_splitter`.
# 
# The target_units, min, max are hadled the usual way via `trait_features`, `trait_aliases_long` and `unit_conversions`.

# In[ ]:


combo = (
    combined_primary_and_secondary
    .with_columns(
        pl.when(
            pl.col("original_term").eq("POCT Blood Ketones") &
            pl.col("result_value_units").is_null()
        )
        .then(
            pl.lit("millimol/L").alias("result_value_units")
        )
        .otherwise(
            pl.col("result_value_units")
        )
    )
)


# # Save `combo` arrow
# 
# This is primary + secondary + handling unitless data.
# It is considered one of the key outputs of the pipeline and therefore stored in `../outputs/reference_combo_files/`

# In[ ]:


get_ipython().run_cell_magic('time', '', '(\n    combo\n    .sink_ipc(\n        AnyPath(\n            PIPELINE_OUTPUTS_REFERENCE_COMBO_FILES_PATH,\n            f"{yr}_{mon}_Combined_all_sources.arrow"\n        )\n    )\n)\n')


# # Import HES data
# 
# Unfortunately, there are differences in file formats for every pull of HES data and each pull needs to be imported individually.

# ## Import HES APC data

# In[ ]:


# We are not considering CC, AE, ECDS, and certanly not OP; included for future-proofing
# hospital_stay_type_enum = pl.Enum(["AE", "APC", "ECDS", "CC", "OP"])
hospital_stay_type_enum = pl.Enum(["APC"])
# adding buffers to list of possible 
region_types_enum =  pl.Enum(list(hospital_stay_type_enum.categories) + ["buffer_before", "buffer_after"])


# In[ ]:


hes_2021_09_APC_txts = (
# (
    pl.scan_csv(
        [
        AnyPath("/genesandhealth/nhsdigital-sublicence-red/DSA__NHSDigitalNHSEngland/2021_09/NIC338864_HES_APC_all_2021_11_25.txt"),
                ],
        separator="|",
        infer_schema=False,
        null_values=[""],
    )
    .TRE
    .filter_with_logging(
        pl.col("ADMIDATE").is_not_null(),
        label="EXCLUDE NULL ADMIDATE"
    )
    .with_columns(
        pl.col("STUDY_ID").alias("pseudo_nhs_number"),
        pl.concat_str(
            pl.col("ADMIDATE"),
            pl.lit("00:00:00"), # no time info provided in this file so we assume earliest time of day
        ).str.to_datetime(format="%F%T").alias("hospital_admission_datetime"),
         pl.concat_str(
            pl.col("DISDATE"),
            pl.lit("23:59:59") # no time info provided in this file so we assume latest time of day
        ).str.to_datetime(format="%F%T").alias("hospital_discharge_datetime"),
        hospital_stay_type=pl.lit("APC", hospital_stay_type_enum),
    )
    .select(
        pl.col("pseudo_nhs_number"),
        pl.col("hospital_admission_datetime"),
        pl.col("hospital_discharge_datetime"),
        pl.col("hospital_stay_type"),
    )

)



# In[ ]:


hes_2023_07_APC_txts = (
    pl.scan_csv(
        [
        AnyPath("/genesandhealth/nhsdigital-sublicence-red/DSA__NHSDigitalNHSEngland/2023_07/HES/*APC*.txt"),

                ],
        separator=",",
        infer_schema=False,
        null_values=[""],
    )
    .TRE
    .filter_with_logging(
        pl.col("ADMIDATE").is_not_null(),
        label="EXCLUDE NULL ADMIDATE"
    )
    .with_columns(
        pl.col("STUDY_ID").alias("pseudo_nhs_number"),
        pl.concat_str(
            pl.col("ADMIDATE"),
            pl.lit("00:00:00"), # no time info provided in this file so we assume earliest time of day
        ).str.to_datetime(format="%F%T").alias("hospital_admission_datetime"),
         pl.concat_str(
            pl.col("DISDATE"),
            pl.lit("23:59:59") # no time info provided in this file so we assume latest time of day
        ).str.to_datetime(format="%F%T").alias("hospital_discharge_datetime"),
        hospital_stay_type=pl.lit("APC", hospital_stay_type_enum),
    )
    .select(
        pl.col("pseudo_nhs_number"),
        pl.col("hospital_admission_datetime"),
        pl.col("hospital_discharge_datetime"),
        pl.col("hospital_stay_type"),
    )
)


# In[ ]:


hes_2023_07_APC_csvs = (
    pl.scan_csv(
        [
        AnyPath("/genesandhealth/nhsdigital-sublicence-red/DSA__NHSDigitalNHSEngland/2023_07/HES/*apc*.csv"),
                ],
        separator=",",
        infer_schema=False,
        null_values=[""],
    )
    .TRE
    .filter_with_logging(
        pl.col("ADMIDATE").is_not_null(),
        label="EXCLUDE NULL ADMIDATE"
    )
    .with_columns(
        pl.col("STUDY_ID").alias("pseudo_nhs_number"),
        pl.concat_str(
            pl.col("ADMIDATE"),
            pl.lit("00:00:00"), # no time info provided in this file so we assume earliest time of day
        ).str.to_datetime(format="%F%T").alias("hospital_admission_datetime"),
         pl.concat_str(
            pl.col("DISDATE"),
            pl.lit("23:59:59") # no time info provided in this file so we assume latest time of day
        ).str.to_datetime(format="%F%T").alias("hospital_discharge_datetime"),
        hospital_stay_type=pl.lit("APC", hospital_stay_type_enum),
    )
    .select(
        pl.col("pseudo_nhs_number"),
        pl.col("hospital_admission_datetime"),
        pl.col("hospital_discharge_datetime"),
        pl.col("hospital_stay_type"),
    )
)


# In[ ]:


hes_2024_10_APC = (
    pl.scan_csv(
        AnyPath("/genesandhealth/nhsdigital-sublicence-red/DSA__NHSDigitalNHSEngland/2024_10/HES/FILE0220459_NIC338864_HES_APC_202399.txt"),
        separator="|",
        infer_schema=False,
        null_values=[""],
    )
    .TRE
    .filter_with_logging(
        pl.col("ADMIDATE").is_not_null(),
        label="EXCLUDE NULL ADMIDATE"
    )
    .with_columns(
        pl.col("STUDY_ID").alias("pseudo_nhs_number"),
        pl.concat_str(
            pl.col("ADMIDATE"),
            pl.lit("00:00:00") # no time info provided in this file so we assume earliest time of day
        ).str.to_datetime(format="%F%T").alias("hospital_admission_datetime"),
         pl.concat_str(
            pl.col("DISDATE"),
            pl.lit("23:59:59") # no time info provided in this file so we assume latest time of day
        ).str.to_datetime(format="%F%T").alias("hospital_discharge_datetime"),
        hospital_stay_type=pl.lit("APC", hospital_stay_type_enum),
    )
    .select(
        pl.col("pseudo_nhs_number"),
        pl.col("hospital_admission_datetime"),
        pl.col("hospital_discharge_datetime"),
        pl.col("hospital_stay_type"),
    )
)


# In[ ]:


hes_2025_03_APC_txts = (
    pl.scan_csv(
        [
        AnyPath("/genesandhealth/nhsdigital-sublicence-red/DSA__NHSDigitalNHSEngland/2025_03/HES/*APC*.txt"),

                ],
        separator="|",
        infer_schema=False,
        null_values=[""],
    )
    .TRE
    .filter_with_logging(
        pl.col("ADMIDATE").is_not_null(),
        label="EXCLUDE NULL ADMIDATE"
    )
    .with_columns(
        pl.col("STUDY_ID").alias("pseudo_nhs_number"),
        pl.concat_str(
            pl.col("ADMIDATE"),
            pl.lit("00:00:00"), # no time info provided in this file so we assume earliest time of day
        ).str.to_datetime(format="%F%T").alias("hospital_admission_datetime"),
         pl.concat_str(
            pl.col("DISDATE"),
            pl.lit("23:59:59") # no time info provided in this file so we assume latest time of day
        ).str.to_datetime(format="%F%T").alias("hospital_discharge_datetime"),
        hospital_stay_type=pl.lit("APC", hospital_stay_type_enum),
    )
    .select(
        pl.col("pseudo_nhs_number"),
        pl.col("hospital_admission_datetime"),
        pl.col("hospital_discharge_datetime"),
        pl.col("hospital_stay_type"),
    )

)


# ### Concatenate HES data
# 
# Currently APC only

# In[ ]:


hes_concat_unfiltered = (
    pl.concat(
       [ 
           hes_2021_09_APC_txts,
           hes_2023_07_APC_txts,
           hes_2023_07_APC_csvs,
           hes_2024_10_APC,
           hes_2025_03_APC_txts,
       ]
    )
   .TRE
    .unique_with_logging()
    .with_columns(
        pl.concat_list(
            pl.col("hospital_stay_type").cast(region_types_enum)
        )
        .alias("region_types")
    )
) # shape: pre-unique (942_302, 5); post-unique (345_918, 5)


# In[ ]:


(
    hes_concat_unfiltered
    .sink_ipc(
        AnyPath(
            COMBINED_DATASETS_ARROW_PATH,
            f"{yr}_{mon}_Combined_HES.arrow"
        )
    )
)


# # Process `combo` to generate all desired output files
# 

# ## Read `combo` back in

# In[ ]:


get_ipython().run_cell_magic('time', '', 'combo = pl.scan_ipc(\n    AnyPath(\n        PIPELINE_OUTPUTS_REFERENCE_COMBO_FILES_PATH,\n        f"{yr}_{mon}_Combined_all_sources.arrow"\n    )\n)\n')


# ## Process `combo` to flag hospitalisation status
# 
# At present we only conside APC episodes >2days.  However, the script is able handle AE/ECDS/CC/OP.  The script can also further be modified to consider "padding".  Padding is not currently used (i.e. set to zero days).
# 
# So:
# ```
# ------------OAPC------------|        APC         |------------OAPC------------  
# ----------|  Buffer before  |        APC         |   Buffer after  |----------  
# ----------|                 TOTAL EXCLUSION ZONE                   |----------
# -- OTEZ --|                                                        |-- OTEZ --
# 
# Padding (not currently used):
# ----| Buffer before  | Pad.  |        APC         | Pad.  | Buffer after |----  
# ```
# 
# * APC = Admitted patient care
# * OAPC = Out of admitted patient care
# * OTEZ = Out of total exclusion zone (herein `out_hospital`)
# 

# ### Read HES data back in

# In[ ]:


hes_concat_unfiltered = (
    pl.scan_ipc(
        AnyPath(
            COMBINED_DATASETS_ARROW_PATH,
            f"{yr}_{mon}_Combined_HES.arrow"
        )
    )
)


# 
# ### Generate final HES data frame
# 
# > "This is where the magic happens." SR, April 2025
# 
# 1. Collect unfiltered data.
# 2. Coalesce overlapping admission windows (including de-duplication).
# 3. (Optional) filter for HES type.  At present we only import APC data.
# 4. Only accept APC episodes >2 days in duration.
# 5. Extend accepted episodes by buffer period.
# 6. Split overlapping intervals and re-merge.
# 

# In[ ]:


hes_final_admission_windows = (
    hes_concat_unfiltered
    .pipe(split_overlapping_intervals_and_remerge,
         start_date_column="hospital_admission_datetime",
         end_date_column="hospital_discharge_datetime"
         )
# Filter not required at present as we only import APC data   
#     .filter(
#         pl.col("hospital_stay_type").eq("APC")
#     )
    .with_columns(
        pl.col("hospital_admission_datetime").dt.round("1d").dt.date().alias("start_date"),
        pl.col("hospital_discharge_datetime").dt.round("1d").dt.date().alias("end_date"),
    )
#     .pipe(lambda _lf: display(_lf.collect()) or _lf)
    .with_columns(
        ((pl.col("end_date") - pl.col("start_date")) )
        .alias("admission_duration")
    )
    .filter(
        pl.col("admission_duration") > pl.duration(days=2)

    )
#     .pipe(lambda _lf: display(_lf.sort("admission_duration").collect()) or _lf)
    .pipe(
        add_buffers,
        id_column="pseudo_nhs_number",
        
    )
    .pipe(
        split_overlapping_intervals_and_remerge
    )
    .sort(["pseudo_nhs_number", "start_date"])
)


# In[ ]:


region_options = list(region_types_enum.categories) 

subsets = list(chain.from_iterable(combinations(region_options, r) for r in range(len(region_options) + 1)))

hes_columns_for_joining = (
    pl.LazyFrame(
        [{**{col: (col in subset) for col in region_options}, "subset": subset } for subset in subsets]
    )
    .with_columns(
        pl.when(pl.col("subset").list.len() > 0)
        .then(
            pl.col("subset").cast(pl.List(region_types_enum)).alias("region_types")
        )
    )
    .select(
        pl.col("region_types"),
        pl.all_horizontal(IN_APC_ONLY).fill_null(False).alias("IN_APC_ONLY"),
        pl.all_horizontal(IN_APC_ANY).fill_null(False).alias("IN_APC_ANY"),
        pl.all_horizontal(IN_BUFFER_BEFORE_ONLY).fill_null(False).alias("IN_BUFFER_BEFORE_ONLY"),
        pl.all_horizontal(IN_BUFFER_BEFORE_ANY).fill_null(False).alias("IN_BUFFER_BEFORE_ANY"),
        pl.all_horizontal(IN_BUFFER_AFTER_ONLY).fill_null(False).alias("IN_BUFFER_AFTER_ONLY"),
        pl.all_horizontal(IN_BUFFER_AFTER_ANY).fill_null(False).alias("IN_BUFFER_AFTER_ANY"),
        pl.all_horizontal(IN_BUFFERS_ONLY).fill_null(False).alias("IN_BUFFERS_ONLY"),
        pl.all_horizontal(IN_BUFFERS_ANY).fill_null(False).alias("IN_BUFFERS_ANY"),
        pl.all_horizontal(IN_TOTAL_EXCLUSION_ZONE).fill_null(False).alias("IN_TOTAL_EXCLUSION_ZONE"),
        pl.all_horizontal(OUT_OF_APC).fill_null(False).alias("OUT_OF_APC"),
        pl.all_horizontal(OUT_OF_TOTAL_EXCLUSION_ZONE).fill_null(False).alias("OUT_OF_TOTAL_EXCLUSION_ZONE"),
    )
)
    
# hes_columns_for_joining.collect()


# ## Join `hes_final_admission_windows` to `combo`
# 
# This adds a column to the `combo` Dataframe which flags the type of the episode that the test date falls in.  At present, this can be:
# 
# * APC
# * Buffer before
# * Buffer after
# 
# All other types (e.g. OUT_OF_TOTAL_EXCLUSION_WINDOW) can be defined from above (see "HES filters" section).

# ### Define a look-up table from test date to HES window type

# In[ ]:


combo_dates_to_HES_region_lookup = (
    combo
    .select(
        pl.col("pseudo_nhs_number"),
        pl.col("test_date")
    )
    .unique()
    .join_where(
        hes_final_admission_windows,
        pl.col("pseudo_nhs_number").eq(pl.col("pseudo_nhs_number_right"))
        & pl.col("test_date").is_between(pl.col("start_date"), pl.col("end_date"), closed="left")
    )
    .select(
        pl.col("pseudo_nhs_number"),
        pl.col("test_date"),
#         pl.col("start_date"),  # DEBUG
#         pl.col("end_date"),  # DEBUG
        pl.col("region_types"),
    )
)


# ### Merge look-up table to `combo`

# In[ ]:


NEED_COMBO_HEIGHT = False
if NEED_COMBO_HEIGHT:
    combo.collect().height
# On 2025-04-04: shape: (78_582_474, 8)


# #### Join `region_types` column (dtype = pl.List) to `combo`

# In[ ]:


combo_with_hes_region_types_column = (
    combo
    .join(
        combo_dates_to_HES_region_lookup,
        on=["pseudo_nhs_number", "test_date"],
        how="left",
        validate="m:1"
    )
#    on 2025-04-04 .collect()  # shape: (78_582_474, 9) Matches plain combo height!
)


# #### Join additional catgeorical HES date regions to `combo`
# 
# This may be useful for troubleshooting/regenie/bespoke set/a paper on impact of hospitalisation but at present not needed for pipeline which will use `region_types` (a polars `List`) column only.

# In[ ]:


# combo_with_multiple_hes_columns = (
#     combo_with_hes_region_types_column
#     .join(
#         hes_columns_for_joining,
#         on="region_types",
#         how="left",
#         # The default .join behaviour of `nulls_equal` *should* be =False. We want =True behaviour. 
#         # Seems to work without setting nulls_equal=True ?bug
#         # Set here explicitly
#         nulls_equal=True
#     )
# )


# #### Optional visualisation of individual's hospitalisation period(s)

# In[ ]:


PLOT_ADMISSIONS_GRAPH_FOR_INDIVIDUALS = False
individual_ids = [
    # Enter pseudoNHSnumbers here
    # see ./data/good_ids_for_APC_plotting.txt
]

n_individuals = len(individual_ids)
total_graph_height = n_individuals*66

if PLOT_ADMISSIONS_GRAPH_FOR_INDIVIDUALS:    
    (
        alt.Chart(
            hes_final_admission_windows.collect()
            .filter(
                pl.col("pseudo_nhs_number").is_in(individual_ids)
            )
            .with_columns(
                pl.col("region_types").list.eval(pl.element().cast(pl.Utf8)).list.join(", "),
                (pl.col("end_date") - pl.col("start_date")).dt.total_days().alias("admission_duration")
            ),
            width="container",
            height=total_graph_height,
        )
        .mark_bar()
        .encode(
            alt.X("start_date"),
            alt.X2("end_date"),
            alt.Y("pseudo_nhs_number"),
            alt.Color("region_types"),
            alt.Size("region_types").scale(
                range=[
                    total_graph_height/n_individuals,
                    total_graph_height/10/n_individuals
                ],
    #             domain=[
    #                 "APC"
    #             ]
            ),
            tooltip=[
                alt.Tooltip("start_date:T", title="Start"),
                alt.Tooltip("end_date:T", title="End"),
                alt.Tooltip("admission_duration:O", title="Days"),
                alt.Tooltip("region_types:N", title="Region Types")
            ]
        )
        .interactive()
        .save(
            AnyPath(
                PIPELINE_OUTPUTS_PATH, f"hes_admissions_plot.html"
            ),
            inline=True
        )
    )


# In[ ]:


PLOT_ADMISSION_AND_BUFFERS_LENGTH_HISTOGRAM = False

if PLOT_ADMISSION_AND_BUFFERS_LENGTH_HISTOGRAM:
    display(
        alt.Chart(
            hes_final_admission_windows
            .with_columns(
                pl.col("region_types").list.eval(pl.element().cast(pl.Utf8)).list.join(", "),
                (pl.col("end_date") - pl.col("start_date")).dt.total_days().alias("admission_duration")
            )
            .collect()
        )
        .mark_bar()
        .encode(
            alt.X("admission_duration:Q").bin(maxbins=50),
            alt.Y("count()").scale(type="log").sort("-x"),
        )
    )


# ## Produce row count per provenance (Optional)

# In[ ]:


VIEW_PROVENANCE_DISTRIBUTION = False
if VIEW_PROVENANCE_DISTRIBUTION:
    display_with(
        combo_with_hes_region_types_column
        .select(pl.col("provenance").value_counts())
        .unnest("provenance")
        .sort(by="count", descending=True)
        .collect()
    )


# ### Units Counts

# In[ ]:


all_counts = (   
    combo_with_hes_region_types_column
    .select(
        pl.col("pseudo_nhs_number"),
        pl.col("original_term"),
        pl.col("result_value_units"),
    )
    .unique()
    .group_by(
        [
            pl.col("original_term"),
            pl.col("result_value_units"),
        ]
    )
    .agg(
        pl.len().alias("n")
    )
    .sort("n", descending=True)
    
)


# In[ ]:


get_ipython().run_cell_magic('time', '', '(\n    all_counts\n#     .filter(\n#         pl.col("n") > 1 #\xa0It may prove insightful to look at n == 1 during development.\n#     )\n    .sink_csv(\n        AnyPath(\n            PIPELINE_LOGS_PATH,\n            f"{yr}_{mon}_units_counts_all_terms.csv"\n        )\n    )\n)\n')


# In[ ]:


get_ipython().run_cell_magic('time', '', '(\n    all_counts\n#     .filter(\n#         pl.col("n") > 1 #\xa0It may prove insightful to look at n == 1 during development.\n#     )\n    .sink_parquet(\n        AnyPath(\n            PIPELINE_LOGS_PATH,\n            f"{yr}_{mon}_units_counts_all_terms.parquet"\n        )\n    )\n)\n')


# ## Import trait data from input files 
# 
# TRAIT_ALIASES_LONG and TRAIT_FEATURES are the result of many iterations/revisions/manual corrections and clinical decisions.  They were originally derived from the QUANT_R version.  Processing code can be found in older versions of `quant_by_pipeline`.  For extraction of further/new traits consider using `trait_augmenter_SR_v0_3.ipynb`.

# ### Read in `unit_conversions` file
# 
# QUANT_R UNITS_MODIFIED renamed to UNIT_CONVERSIONS.
# 
# 
# 
# e.g. 
# ```
# result_value_units,target,multiplication_factor
# Grams/decilitre,g/L,10
# ```
# 
# We open `UNIT_CONVERSIONS` twice, first to create a list of valid target units (a superset of TRAIT_FEATURES target units).  Secondly to create the unit converting table = `units_converter`.

# #### Obtain units Enum from `unit_conversions` file
# 
# Here we restrict possible **target** units to the ones defined (and approved) in `units_conversions`.  This acts as an additional QC step.  

# In[ ]:


target_units_enum = pl.Enum(
    pl.scan_csv(
        UNIT_CONVERSIONS_PATH
    )
    .select(
        pl.col("target")
    )
    .unique()
    .sort(by="target") # we should think about the sort order.
    .collect()
)

target_units_enum

# len(['%', '10*12/L', '10*9/L', 'L', 'L/min', 'bpm', 'fl', 'g/L', 'kU/L', 'kg', 'kg/m2', 
#      'm', 'mg/L', 'mg/mmol', 'microg/L', 'micromol/L', 'milimol/L', 'miliunits/L', 'ml/min', 
#      'mm/h', 'mmHg', 'mmol/mol', 'nanomol/L', 'ng/L', 'pg', 'pg/ml', 'picomol/L', 'seconds', 
#      'units/L', 'units/week', 'weeks'])


# #### Create the `unit_converter`
# 
# We restrict possible **target** units to the ones defined (and approved) in `units_conversions` using the Enum created in the previous cell.  This acts as an additional QC step.  

# In[ ]:


units_converter = (
    pl.scan_csv(
        UNIT_CONVERSIONS_PATH,
        infer_schema=False,
        schema_overrides={
            "target": target_units_enum, 
            "multiplication_factor": pl.Float64,
        }
    )
)


# #### Additional check on `units_converter`
# 
# This table should not have duplicated rows in it.  We test this with `assert` statement.

# In[ ]:


try:
    assert not any(
        units_converter.collect().is_duplicated()
    )
except Exception as dup:
    print(
        f'''Duplicate row(s) found in `{UNIT_CONVERSIONS_PATH.name}` 
        {
        (
            units_converter
            .filter(
                units_converter.collect().is_duplicated()
            )
            .collect()
        )
        }
        '''
    )
else:
    print(f"Bravo, no duplicates in `{UNIT_CONVERSIONS_PATH.name}`")


# ### Read in trait_aliases_long file
# 
# This is a TRAIT, trait_alias pair file.  With 1:m TRAIT:alias.

# In[ ]:


trait_aliases_long = (
    pl.scan_csv(
        TRAIT_ALIASES_LONG_PATH,
    )
)


# ### Read in trait_features file
# 
# One TRAIT per line, details trait name, target units, min value accepted, max value accepted.
# 
# e.g. `HDL-C,millimol/L,0.05,4.65`

# In[ ]:


trait_features = (
    pl.scan_csv(
        TRAIT_FEATURES_PATH,
        schema_overrides={
            "target_units": target_units_enum
        }
    )
)


# ### Create `traits_denormalised`
# 
# Here we join the Trait : Features file (`traits_feature`) with the Trait : Alias file (`trait_aliases_long`) to create a denormalised table suitable for joining to `combo`.

# In[ ]:


traits_denormalised = (
    trait_features
    .join(
        trait_aliases_long,
        on="trait",
        how="left",
    )
)


# In[ ]:


get_ipython().run_cell_magic('time', '', '## Use this as sanity check and/or to see if any immediate TRAITS worth considering\n## and save output to logs\nCHECK_FOR_UNRECOVERED_TRAITS = False\n\nif CHECK_FOR_UNRECOVERED_TRAITS:\n    combo_traits_anti_case_sensitive = (\n        combo_with_hes_region_types_column\n            .select(pl.col("original_term"))\n            .join(\n                trait_aliases_long,\n                left_on=pl.col("original_term").str.strip_chars(), \n                right_on="alias", \n                how="anti",\n            )    \n        .group_by("original_term")\n        .agg(pl.len())\n        .sort(by="len", descending=True)\n    )\n    \n    (\n        combo_traits_anti_case_sensitive\n        .pipe(lambda _lf: display_with(_lf.collect()) or _lf)\n        .sink_ipc(\n            AnyPath(\n                PIPELINE_LOGS_PATH,\n                f"{yr}_{mon}_unrecovered_traits.arrow"\n            )\n        )\n        \n    )\n')


# ## `combo_strict_trait`
# `combo_strict_trait` is a key dataframe, it combines `combo` with the curated list of traits and their features (i.e. their target_unit and min and max admissible values)

# In[ ]:


combo_strict_trait = (
    combo_with_hes_region_types_column
        .join(
            traits_denormalised,
            left_on=pl.col("original_term"), 
            right_on=pl.col("alias"),
            how="left",
        )
)


# ## `combo_strict_trait_ranged`
# 
# `combo_strict_trait_ranged` is another key dataframe.
# 
# `combo_strict_trait_ranged` uses the unit converter to convert units to target units and check if these are in range.
# 
# **This section includes _ad-hoc_ conversions such as those applied to `HbA1c` values in percentages.**
# 
# In the future we should consider having a similar section to "Handle unitless values".

# In[ ]:


get_ipython().run_cell_magic('time', '', '\nrange_enum = pl.Enum(["below_min", "ok", "above_max"])\n\ncombo_strict_trait_ranged = (\n    combo_strict_trait\n    ### Here we exclude all reading with null units\n    ### There are a number of traits we wish to recover in which either truly have no units\n    ### or in which we assume a unit for nulls\n    ### We deal with unitless traits in a bespoke per trait fashion above \n    ### (e.g. Blood_ketones\' unitless POCT vals)\n    .TRE\n    .filter_with_logging(\n        EXCLUDE_NULL_UNITS,\n        label="EXCLUDE_NULL_UNITS"\n    )\n    \n    # This is where we allow result_value_units to be converted\n    # this allow for both "value modifying converstions" (e.g. nmol -> mmol by divide by 1,000)\n    # and for unit format converstion (e.g. MMOL/MOL -> mmol/mol)\n    .join(units_converter, left_on=["result_value_units", "target_units"], right_on=["result_value_units", "target"], how="left", coalesce=False) # shape: (69_727_105, 15)\n    .with_columns(\n        pl.col("multiplication_factor")\n    )\n    \n    # values conversions according to \n    # 1) HbA1c formula, and\n    # 2) multiplication_factor\n    \n    ### Aim to separate out these two...\n    .with_columns(\n        pl.when(\n            pl.col("trait").eq("HbA1c") & \n            pl.col("result_value_units").is_in(["%", "% total Hb","%Hb", "per cent"]),\n        )\n        .then(\n            (10.93 * pl.col("result") - 23.50)\n        )\n        .when(\n            pl.col("multiplication_factor").ne(1)\n        )\n        .then(\n            (pl.col("result") * pl.col("multiplication_factor"))\n        )\n        .otherwise(\n            pl.col("result")\n        )\n        .alias("final")\n    )\n    \n    # categorise according to min/max range bounds:\n    .with_columns(\n        pl.when(\n            pl.col("final") < pl.col("min")\n        )\n        .then(\n            pl.lit("below_min").cast(range_enum)\n        )\n        .when(\n            pl.col("final").is_between(\n                pl.col("min"), \n                pl.col("max"), \n                closed="both"\n            )\n        )\n        .then(\n            pl.lit("ok").cast(range_enum)\n        )\n        .when(\n            pl.col("final") > pl.col("max")\n        )\n        .then(\n            pl.lit("above_max").cast(range_enum)\n        )\n        .otherwise(\n            None\n        )\n        .alias("range_position")\n\n    )\n    .with_columns(\n        pl.col("final").replace(0,1e-10).log10().alias("final_log10")\n    )\n)\n')


# ### Restrict to valid pseudoNHS numbers and valid demographics
# 
# Process:
# 1. Import mega_linkage file $^{TM}$
# 2. Restrict data to pseudoNHS in mega_linkage file $^{TM}$
# 3. Use mega_linkage file $^{TM}$ to link to S1QST for DOB
# 

# #### …and filter to valid pseudoNHS
# 
# There are approximately 1,000 rows excluded by this pseudoNHS validation.
# 
# Possible reasons: 
# 1. subject asked to be removed/withdrawn
# 2. subject died
# 3. subject had multiple pseudoNHS which have been merged
# 
# In quant_py, we now use the DvH's `2025_02_10__MegaLinkage_forTRE.csv` $^{TM}$

# In[ ]:


valid_pseudo_nhs_numbers = (
    pl.scan_csv(
        MEGA_LINKAGE_PATH,
        infer_schema=False,
    )
    .rename({"pseudonhs_2024-07-10":"pseudo_nhs_number"})
    .select(
        pl.col("pseudo_nhs_number"),
    )
    .TRE
    .filter_with_logging(
        pl.col("pseudo_nhs_number").is_not_null(),
        label="pseudo_nhs_number.is_not_null()"
    )
    .TRE
    .unique_with_logging(
        "pseudo_nhs_number",
    )
)


# In[ ]:


## This step creates DOBs for all people who have a questionnaire (with an Oragene_ID)

s1qst_dob_and_gender = (
    pl.scan_csv(
        AnyPath(
            "/",
            "genesandhealth",
            "library-red",
            "genesandhealth",
            "phenotypes_rawdata",
            "QMUL__Stage1Questionnaire",
            "2025_01_24__S1QSTredacted.csv" # The new one without future births
        ),
        infer_schema=False
    )
    .select(
        pl.col("S1QST_Oragene_ID"),
        pl.col("S1QST_Gender"),
        pl.col("S1QST_MM-YYYY_ofBirth"),
    )
    .TRE
    .filter_with_logging(
        pl.col("S1QST_MM-YYYY_ofBirth").ne("NA"),
        label="EXCLUDING `NA` DATE"
    )
    .with_columns(
        pl.concat_str(pl.lit("01-"), pl.col("S1QST_MM-YYYY_ofBirth")).str.to_date(format="%d-%m-%Y").alias("dob"),
        pl.col("S1QST_Gender")
            .replace_strict({"1":"M","2":"F"})
            .cast(pl.Enum(["F","M"])).alias("gender"),
    )
    .with_columns(
        pl.col("dob").dt.year().alias("year_of_birth")        
    )
#     .TRE
#     .unique_with_logging()
    .TRE
    .unique_with_logging(subset=["S1QST_Oragene_ID"])
    
    .select(
        pl.col("S1QST_Oragene_ID").alias("OrageneID"),
        pl.col("dob"),
        pl.col("gender")
    )
#     .group_by("S1QST_Oragene_ID")
#     .agg(
#         pl.col("S1QST_MM-YYYY_ofBirth").n_unique().alias("MM-YYYY_ofBirth_count")
#     )
#     .sort("MM-YYYY_ofBirth_count")
)


# In[ ]:


valid_demographics = (
    pl.scan_csv(
        MEGA_LINKAGE_PATH,
        infer_schema=False,
    )
    .rename({"pseudonhs_2024-07-10":"pseudo_nhs_number"})
    .select(
        pl.col("pseudo_nhs_number"),
        pl.col("OrageneID")
    )
    .TRE
    .filter_with_logging(
        pl.col("pseudo_nhs_number").is_not_null(),
        label="pseudo_nhs_number.is_not_null()"
    )
    .TRE
    .unique_with_logging(
        "pseudo_nhs_number",
    )
    .TRE
    .join_with_logging(
        s1qst_dob_and_gender,
        on="OrageneID",
        how="left",
        label="Add dob and gender columns"
    )
)


# In[ ]:


valid_regenie_55k = (
    pl.scan_csv(
        MEGA_LINKAGE_PATH,
        infer_schema=False,
        new_columns=[
            "OrageneID",
            "Number of OrageneIDs with this NHS number (i.e. taken part twice or more)",
            "s1qst_gender",
            "HasValidNHS",
            "pseudo_nhs_number",
            "gsa_id",
            "44028exomes_release_2023-JUL-07",
            "exome_id",
        ]
    )
    .TRE
    .filter_with_logging(
        pl.col("exome_id").is_not_null(),
        pl.col("pseudo_nhs_number").is_not_null(), # there are some rows with NON-NULL exome_id but NULL pseudo_nhs_number
        label="Only include NON-NULL exome_id and NON-NULL pseudo_nhs_number for 55k Regenie"
    )
    .TRE
    .filter_with_logging(
        pl.col("OrageneID").is_not_null(),
        label="Sanity check to ensure no NULL OrageneID. row count should remain unchanged"
    )
    .TRE
    .unique_with_logging(
        ["pseudo_nhs_number"],
        label="Sanity check: row count should remain unchanged when uniquing by pseudo_nhs_number"
    )
    .TRE
    .unique_with_logging(
        ["OrageneID"],
        label="Sanity check: row count should remain unchanged when uniquing by OrageneID"
    )
)


# In[ ]:


valid_regenie_51k = (
    pl.scan_csv(
        MEGA_LINKAGE_PATH,
        infer_schema=False,
        new_columns=[
            "OrageneID",
            "Number of OrageneIDs with this NHS number (i.e. taken part twice or more)",
            "s1qst_gender",
            "HasValidNHS",
            "pseudo_nhs_number",
            "gsa_id",
            "44028exomes_release_2023-JUL-07",
            "exome_id",
        ]
    )
    .TRE
    .filter_with_logging(
        pl.col("gsa_id").is_not_null(),
        pl.col("pseudo_nhs_number").is_not_null(), # there are some rows with NON-NULL exome_id but NULL pseudo_nhs_number
        label="Only include NON-NULL gsa_id and NON-NULL pseudo_nhs_number for 51k Regenie"
    )
    .TRE
    .filter_with_logging(
        pl.col("OrageneID").is_not_null(),
        label="Sanity check to ensure no NULL OrageneID. row count should remain unchanged"
    )
    .TRE
    .unique_with_logging(
        ["pseudo_nhs_number"],
        label="Sanity check: row count should remain unchanged when uniquing by pseudo_nhs_number"
    )
    .TRE
    .unique_with_logging(
        ["OrageneID"],
        label="Sanity check: row count should remain unchanged when uniquing by OrageneID"
    )
)


# ## Now we remove pseudo_nhs number not in mega-linkage file

# In[ ]:


combo_strict_trait_ranged_valid_pseudo_nhs_nums = (
    combo_strict_trait_ranged
    .join(
        valid_pseudo_nhs_numbers,
        on="pseudo_nhs_number",
        how="semi"
    )
)

## On 2025-05-25:
# pre-pseudo_nhs validation = 49_562 volunteers
# post-pseudo_nhs validation = 49_541 volunteers (i.e. 21 volunteer data excluded)


# ## Combo questionnaire
# 
# This adds `dob` and `gender` info to `combo_strict_trait_ranged_valid_pseudo_nhs_nums` for anyone who has completed a questionnaire.

# In[ ]:


get_ipython().run_cell_magic('time', '', 'combo_strict_trait_ranged_valid_pseudo_nhs_nums_plus_demographics = (\n    combo_strict_trait_ranged_valid_pseudo_nhs_nums\n    .TRE\n    .join_with_logging(\n        valid_demographics, \n        on="pseudo_nhs_number",\n        how="left",\n        label="Adding exome id and OrageneID"\n    )\n    .with_columns(\n        ((pl.col("test_date") - pl.col("dob")).dt.total_days() / 365.25).alias("age_at_test"),\n        pl.col("final").replace(0,1e-10).log10().alias("value_log10")\n    )\n    .TRE\n    .filter_with_logging(\n        EXCLUDE_READINGS_WITH_VALUES_OUTSIDE_EXPECTED_RANGE,\n        label="EXCLUDE_READINGS_WITH_VALUES_OUTSIDE_EXPECTED_RANGE"\n    )\n    .TRE\n    .filter_with_logging(\n        EXCLUDE_READINGS_WITH_IMPLAUSIBLE_DATES,\n        label="EXCLUDE_READINGS_WITH_IMPLAUSIBLE_DATES"\n    )\n    .TRE\n    .filter_with_logging(\n        EXCLUDE_READINGS_WITH_INDIVS_UNDER_SIXTEEN,\n        label="EXCLUDE_READINGS_WITH_INDIVS_UNDER_SIXTEEN"\n    )\n\n#     .collect()\n)\n\n(\n    combo_strict_trait_ranged_valid_pseudo_nhs_nums_plus_demographics\n    .sink_parquet(\n        AnyPath(\n            PIPELINE_OUTPUTS_REFERENCE_COMBO_FILES_PATH,\n            f"{yr}_{mon}_Combined_traits_NHS_and_demographics_restricted_pre_10d_windowing.parquet"\n        )\n    )\n\n)\n## Pre-ranging\n# [Adding exome id and OrageneID] Join type: LEFT\n# [Adding exome id and OrageneID] Left: 74329771 rows, Right: 57846 rows -> After: 74329771 rows (row count unchanged) (0.0%)\n# [EXCLUDE_READINGS_WITH_IMPLAUSIBLE_DATES] Before filter: 74329771 rows, After filter: 74327434 rows (-0.0%)\n# [EXCLUDE_READINGS_WITH_INDIVS_UNDER_SIXTEEN] Before filter: 74327434 rows, After filter: 73784492 rows (-0.7%)\n# CPU times: user 4min 37s, sys: 2min 42s, total: 7min 19s\n# Wall time: 1min 42s\n\n')


# # 10 day windows, every 11 days

# In[ ]:


combo_strict_trait_ranged_valid_pseudo_nhs_nums_plus_demographics_with_10d_windowing = (
    combo_strict_trait_ranged_valid_pseudo_nhs_nums_plus_demographics
    .filter(
        pl.col("final").is_not_null()
    )
    .with_columns(
        pl.col("test_date").alias("window_date")
    )
    .sort("window_date")
    .group_by_dynamic(
        index_column="window_date",
        every="11d",
        period="10d",
        closed="both",
        group_by=["pseudo_nhs_number", "trait", "final"]
    )
    .agg(
        pl.all().first()
    )
    .rename(
        {
            "target_units":"unit",
            "final":"value",
            "test_date":"date",
            'range_position':"minmax_outlier",
        }
    )
    .select(
        TARGET_COMBO_POST_10D_WINDOWING_COLUMNS
    )
)


# # Write combo post 10d windowing as parquet
# By processed we mean:
# 1. Trait column populated if appropriate
# 2. Target units column populated if appropriate
# 3. Units converted to target units where necessary and possible
# 4. Values checked to see if within accepted values (min-max outliers identified)
# 5. Invalid/opted-out (pseudo-)NHS numbers excluded
# 6. 10-day results window "deduplication" applied
# 
# **This file can be read in to shortcut all of the processing performed above**
# 

# In[ ]:


get_ipython().run_cell_magic('time', '', 'WRITE_COMBO_POST_10D_WINDOWING_FILE = True\n\nif WRITE_COMBO_POST_10D_WINDOWING_FILE:\n    (\n        combo_strict_trait_ranged_valid_pseudo_nhs_nums_plus_demographics_with_10d_windowing\n        .sink_parquet(\n            AnyPath(\n                COMBO_POST_10D_WINDOWING_PATH,\n                f"{yr}_{mon}_Combined_traits_NHS_and_demographics_restricted_post_10d_windowing.parquet"\n            )\n        )\n    )\n')


# ## Read COMBO_PROCESSED back in (from parquet)

# In[ ]:


combo_strict_trait_ranged_valid_pseudo_nhs_nums_plus_demographics_with_10d_windowing = (
    pl.scan_parquet(
        AnyPath(
            COMBO_POST_10D_WINDOWING_PATH,
            f"{yr}_{mon}_Combined_traits_NHS_and_demographics_restricted_post_10d_windowing.parquet"
        )
    )
)


# # Write output files

# Most output files are created here.
# 
# For simplification, all output directories are specified here.

# ## Create output sub- and subsub- directories if required
# 
# There are 3 subdirectories in `/outcomes/`:
# 1. individual_trait_files
# 2. individual_trait_plots
# 3. regenie_files
# 
# We create 3 datasets per trait **\[HES\]**:
# * `in_hospital`: only data within an APC+BUFFER
# * `out_hospital`: only data out of APC+BUFFER
# * `all`: all collected data regardless of hospitalisation status
# 
# So each subdirectory has 3 subsubdirectories.
# 
# The directories are defined here as they should stay the same regardless of path to them.
# 

# ### `individual_trait_files`
# 
# each of the subsubdiretories contains 2 files per trait (here trait = `Weight`):
# 1. `../individual_trait_files/\[HES\]/YYYY_MM_Weight_\[HES\]_readings_at_unique_timepoints.csv`: one validated result per row, with `pseudo_nhs_number, trait, unit, value, date, gender, age_at_test, minmax_outlier` columns.
# 2. `../individual_trait_files/\[HES\]/YYYY_MM_Weight_\[HES\]_per_individual_stats.csv`: one row per volunteer with columns `pseudo_nhs_number, trait, median, mean, max, min, earliest, latest, n`.
# 
# i.e. 6 `individual_trait` files per trait.

# In[ ]:


PIPELINE_INDIVIDUAL_TRAIT_FILES_IN_HOSPITAL_PATH = (
    AnyPath(
        PIPELINE_INDIVIDUAL_TRAIT_FILES_PATH,
        "in_hospital"
    )
)

PIPELINE_INDIVIDUAL_TRAIT_FILES_OUT_HOSPITAL_PATH = (
    AnyPath(
        PIPELINE_INDIVIDUAL_TRAIT_FILES_PATH,
        "out_hospital"
    )
)

PIPELINE_INDIVIDUAL_TRAIT_FILES_ALL_PATH = (
    AnyPath(
        PIPELINE_INDIVIDUAL_TRAIT_FILES_PATH,
        "all"
    )
)


# In[ ]:


PIPELINE_INDIVIDUAL_TRAIT_FILES_IN_HOSPITAL_PATH.mkdir(parents=True, exist_ok=True)
PIPELINE_INDIVIDUAL_TRAIT_FILES_OUT_HOSPITAL_PATH.mkdir(parents=True, exist_ok=True)
PIPELINE_INDIVIDUAL_TRAIT_FILES_ALL_PATH.mkdir(parents=True, exist_ok=True)


# ### `individual_trait_plots`
# 
# One .svg file per trait ('in_hospital', 'out_hospital', 'all')

# In[ ]:


PIPELINE_INDIVIDUAL_TRAIT_PLOTS_IN_HOSPITAL_PATH = (
    AnyPath(
        PIPELINE_INDIVIDUAL_TRAIT_PLOTS_PATH,
        "in_hospital"
    )
)

PIPELINE_INDIVIDUAL_TRAIT_PLOTS_OUT_HOSPITAL_PATH = (
    AnyPath(
        PIPELINE_INDIVIDUAL_TRAIT_PLOTS_PATH,
        "out_hospital"
    )
)

PIPELINE_INDIVIDUAL_TRAIT_PLOTS_ALL_PATH = (
    AnyPath(
        PIPELINE_INDIVIDUAL_TRAIT_PLOTS_PATH,
        "all"
    )
)


# In[ ]:


PIPELINE_INDIVIDUAL_TRAIT_PLOTS_IN_HOSPITAL_PATH.mkdir(parents=True, exist_ok=True)
PIPELINE_INDIVIDUAL_TRAIT_PLOTS_OUT_HOSPITAL_PATH.mkdir(parents=True, exist_ok=True)
PIPELINE_INDIVIDUAL_TRAIT_PLOTS_ALL_PATH.mkdir(parents=True, exist_ok=True)


# ### `regenie_files`

# In[ ]:


PIPELINE_OUTPUTS_REGENIE_PATH = (
    AnyPath(
        PIPELINE_OUTPUTS_PATH,
        "regenie",
    )
)

PIPELINE_OUTPUTS_REGENIE_IN_HOSPITAL_PATH = (
    AnyPath(
        PIPELINE_OUTPUTS_REGENIE_PATH,
        "in_hospital",
    )
)

PIPELINE_OUTPUTS_REGENIE_OUT_HOSPITAL_PATH = (
    AnyPath(
        PIPELINE_OUTPUTS_REGENIE_PATH,
        "out_hospital",
    )
)

PIPELINE_OUTPUTS_REGENIE_ALL_PATH = (
    AnyPath(
        PIPELINE_OUTPUTS_REGENIE_PATH,
        "all",
    )
)

PIPELINE_OUTPUTS_REGENIE_COVARIATES_FILES_PATH = (
    AnyPath(
        PIPELINE_OUTPUTS_REGENIE_PATH,
        "covariate_files",
    )
)


# In[ ]:


PIPELINE_OUTPUTS_REGENIE_PATH.mkdir(parents=True, exist_ok=True)
PIPELINE_OUTPUTS_REGENIE_IN_HOSPITAL_PATH.mkdir(parents=True, exist_ok=True)
PIPELINE_OUTPUTS_REGENIE_OUT_HOSPITAL_PATH.mkdir(parents=True, exist_ok=True)
PIPELINE_OUTPUTS_REGENIE_ALL_PATH.mkdir(parents=True, exist_ok=True)
PIPELINE_OUTPUTS_REGENIE_COVARIATES_FILES_PATH.mkdir(parents=True, exist_ok=True)


# ## Write `readings_at_unique_timepoints.csv` per trait

# In[ ]:


for region_category, FILTER in {
    "in_hospital": IN_TOTAL_EXCLUSION_ZONE,
    "out_hospital": OUT_OF_TOTAL_EXCLUSION_ZONE,
    "all": ( True )
}.items():
    for (trait, ), df in (
        combo_strict_trait_ranged_valid_pseudo_nhs_nums_plus_demographics_with_10d_windowing
        .filter(
            FILTER
        )
        .select(
            TARGET_TRAIT_READINGS_AT_INDIVIDUAL_TIMEPOINTS_COLUMNS
        )
        .collect()
        .group_by("trait")):
                df.write_csv(
                    AnyPath(
                        PIPELINE_INDIVIDUAL_TRAIT_FILES_PATH,
                        region_category,
                        f"{yr}_{mon}_{trait}_{region_category}_readings_at_unique_timepoints.csv"
                    )
                )


# ## Write `raw_all.csv` files per trait
# 
# These are before windowing.
# 
# In practice, don't think needed, not generated for v010.

# In[ ]:


# grouped_traits_for_raw_all = (
#     combo_strict_trait_ranged_valid_pseudo_nhs_nums_plus_demographics
#     .rename(
#         {"target_units":"unit","final":"value","test_date":"date"})
#     .select(TARGET_TRAIT_RAW_ALL_COLUMNS)
#     .collect()
#     .group_by("trait")
# )


# In[ ]:


# %%time
# for (trait, ), df in grouped_traits_for_raw_all:
#     df.write_csv(
#         AnyPath(
#             PIPELINE_INDIVIDUAL_TRAIT_FILES_PATH,
#             f"{yr}_{mon}_{trait}_raw_all.csv"
#         )
#     )


# ## Write `per_individual_stats.csv` per trait 

# In[ ]:


get_ipython().run_cell_magic('time', '', 'for region_category, FILTER in {\n    "in_hospital": IN_TOTAL_EXCLUSION_ZONE,\n    "out_hospital": OUT_OF_TOTAL_EXCLUSION_ZONE,\n    "all": ( True )\n}.items():\n    per_trait_per_individual_stats = (\n        combo_strict_trait_ranged_valid_pseudo_nhs_nums_plus_demographics_with_10d_windowing\n        .filter(FILTER)\n        .group_by(["pseudo_nhs_number", "trait", "minmax_outlier"])\n        .agg(\n            pl.col("value").median().alias("median"),\n            pl.col("value").mean().alias("mean"),\n            pl.col("value").max().alias("max"),\n            pl.col("value").min().alias("min"),\n            pl.col("value").filter(pl.col("date").eq(pl.col("date").min())).first().alias("earliest"),\n            pl.col("value").filter(pl.col("date").eq(pl.col("date").max())).first().alias("latest"),\n            pl.count("value").alias("n")\n        )\n\n        .select(\n            *TARGET_TRAIT_PER_INDIVIDUAL_STATS_COLUMNS\n        )\n        .collect()\n    )\n    \n    for (trait, ), df in per_trait_per_individual_stats.group_by("trait"):\n        df.write_csv(\n            AnyPath(\n                PIPELINE_INDIVIDUAL_TRAIT_FILES_PATH,\n                region_category,\n                f"{yr}_{mon}_{trait}_{region_category}_per_individual_stats.csv"\n            )\n        )\n')


# ## Generate trait plots
# 
# 

# In[ ]:


post_qc_histogram_data = (
    combo_strict_trait_ranged_valid_pseudo_nhs_nums_plus_demographics_with_10d_windowing
    .with_columns(
        pl.col("gender").cast(pl.Utf8).alias("gender"),
        pl.col("value").replace(0,1e-10).log10().alias("log_x")
    )
    .filter(
        pl.col("trait").is_not_null(),
    )
    .select(
        pl.col("pseudo_nhs_number"),
        pl.col("trait"),
        pl.col("value"),
        pl.col("unit"),
        pl.col("log_x"),
        pl.col("gender"),
        pl.col("region_types"),
    )
    .collect()
)


# In[ ]:


post_qc_histogram_data


# In[ ]:


(
    alt.Chart(
        post_qc_histogram_data
        .filter(
            pl.col("trait").eq("BMI")
        )
    )
    .mark_bar()
    .encode(
        alt.X("log_x:Q").bin(maxbins=48).title(f"BMI"),#.scale(type="log"),
        alt.Y("count()").title(None),
        alt.Color("gender:N")
        .scale(
            range=[GNH_PALETTE["EMERALD_GREEN"], GNH_PALETTE["COBALT_BLUE"]],
            domain=["F","M"]
        ),
        row="gender:N",             
    )
)


# In[ ]:


def gender_plot_for_trait(trait:str, df: pl.DataFrame, region_category: str) -> alt.Chart():
    plot_data = (
        df
    )
    median_value = plot_data.get_column("value").median()
    mean_value = plot_data.get_column("value").mean()
    min_value = plot_data.get_column("value").min()
    max_value = plot_data.get_column("value").max()
    n = plot_data.get_column("pseudo_nhs_number").n_unique()
    observation_count = plot_data.shape[0]
    units = plot_data.get_column("unit").first()
    ###Save  Plot
    (
        alt.Chart(
            plot_data,
            title=alt.Title(
                 f"{trait} [{region_category}]",
                subtitle=[
                    f"Median: {median_value: .1f}", 
                    f"Mean: {mean_value: .1f}",
                    f"Min: {min_value}",
                    f"Max: {max_value}",
                    f"n individuals: {n:,}",
                    f"n observations: {observation_count:,}",
                ],
                anchor="start",
                frame="group",
            )
        )
        .mark_bar(stroke="black")
        .encode(
            alt.X("log_x:Q").bin(maxbins=48).title(f"Log10 {trait} ({units})"),#.scale(type="log"),
            alt.Y("count()").title(None),
            alt.Color("gender:N")
            .scale(
                range=[GNH_PALETTE["EMERALD_GREEN"], GNH_PALETTE["COBALT_BLUE"]],
                domain=["F","M"]
            ),
            row="gender:N"
            
        )
        .properties(height=200, width=800)
    ).save(
        AnyPath(
            PIPELINE_INDIVIDUAL_TRAIT_PLOTS_PATH,
            region_category,
            f"{yr}_{mon}_{trait.replace(' ','_')}_{region_category}.svg"
        ),
        format="svg"
    )


# In[ ]:


get_ipython().run_cell_magic('time', '', 'for region_category, FILTER in {\n    "in_hospital": IN_TOTAL_EXCLUSION_ZONE,\n    "out_hospital": OUT_OF_TOTAL_EXCLUSION_ZONE,\n    "all": ( True )\n}.items():\n    print(region_category)\n    for (trait, ), df in post_qc_histogram_data.group_by("trait"):\n        df_filtered = df.filter(FILTER)\n        if df_filtered.is_empty():\n            print(f"\\t{trait} {region_category}: No readings, skipping...")\n            continue\n        print(f"\\t{trait}")\n        gender_plot_for_trait(trait, df_filtered, region_category=region_category)\n')


# ### Write regenie_51koct2024_GSA_Topmed_pheno TSV
# 

# In[ ]:


get_ipython().run_cell_magic('time', '', 'regenie_51k_data = (\n    combo_strict_trait_ranged_valid_pseudo_nhs_nums_plus_demographics_with_10d_windowing\n    .join(\n        valid_regenie_51k,\n        on="pseudo_nhs_number",\n        how="inner"\n    )\n    .collect()\n)\n\nfor region_category, FILTER in {\n    "in_hospital": IN_TOTAL_EXCLUSION_ZONE,\n    "out_hospital": OUT_OF_TOTAL_EXCLUSION_ZONE,\n    "all": ( True )\n}.items():\n    for trait_name, group in (\n        regenie_51k_data\n        .filter(FILTER)\n        .group_by("trait")\n    ):\n        trait = trait_name[0].replace(" ","_")\n\n        (\n            group.select(\n                pl.lit("1").alias("FID"),\n                pl.col("gsa_id").alias("IID"),\n                pl.col("value").median().over("gsa_id").alias(f"{trait}.median"),\n                pl.col("value").min().over("gsa_id").alias(f"{trait}.min"),\n                pl.col("value").max().over("gsa_id").alias(f"{trait}.max"),\n            )\n            .unique()\n            .write_csv(\n                AnyPath(\n                    PIPELINE_OUTPUTS_REGENIE_PATH,\n                    region_category,\n                    f"{yr}_{mon}_{trait}_{region_category}_regenie_51koct2024_GSA_Topmed_pheno.tsv"\n                ),\n                separator="\\t",\n\n            )\n        )\n')


# ### Write regenie_55k_BroadExomeIDs_pheno TSV

# In[ ]:


get_ipython().run_cell_magic('time', '', 'regenie_55k_data = (\n    combo_strict_trait_ranged_valid_pseudo_nhs_nums_plus_demographics_with_10d_windowing\n    .join(\n        valid_regenie_55k,\n        on="pseudo_nhs_number",\n        how="inner"\n    )\n    .collect()\n)\n\nfor region_category, FILTER in {\n    "in_hospital": IN_TOTAL_EXCLUSION_ZONE,\n    "out_hospital": OUT_OF_TOTAL_EXCLUSION_ZONE,\n    "all": ( True )\n}.items():\n    for trait_name, group in regenie_55k_data.group_by("trait"):\n        trait = trait_name[0].replace(" ","_")\n\n        (\n            group.select(\n                pl.lit("1").alias("FID"),\n                pl.col("exome_id").alias("IID"),\n                pl.col("value").median().over("exome_id").alias(f"{trait}.median"),\n                pl.col("value").min().over("exome_id").alias(f"{trait}.min"),\n                pl.col("value").max().over("exome_id").alias(f"{trait}.max"),\n            )\n            .unique()\n            .write_csv(\n                AnyPath(\n                    PIPELINE_OUTPUTS_REGENIE_PATH,\n                    region_category,\n                    f"{yr}_{mon}_{trait}_{region_category}_regenie_55k_BroadExomeIDs_pheno.tsv"\n                ),\n                separator="\\t"\n            )\n        )\n')


# ### Write regenie_51koct2024_GSA_Topmed_pheno COVARIATE MEGAWIDE

# In[ ]:


get_ipython().run_cell_magic('time', '', 'for region_category, FILTER in {\n    "in_hospital": IN_TOTAL_EXCLUSION_ZONE,\n    "out_hospital": OUT_OF_TOTAL_EXCLUSION_ZONE,\n    "all": ( True )\n}.items():\n    # Create a dictionary of region_category filtered dataframe\n    regenie_51k_data_megawide_age_at_test_partitioned_dict = (\n        combo_strict_trait_ranged_valid_pseudo_nhs_nums_plus_demographics_with_10d_windowing\n        .filter(FILTER)\n        .join(\n            valid_regenie_51k,\n            on="pseudo_nhs_number",\n            how="inner"\n        )\n        .select(\n            pl.lit("1").alias("FID"),\n            pl.col("gsa_id").alias("IID"),\n            pl.col("trait"),\n            pl.col("age_at_test").round(1).alias("AgeAtTest"),\n            pl.col("age_at_test").pow(2).round(1).alias("AgeAtTest_Squared")\n        )\n        .collect()\n        .partition_by("trait", as_dict=True)\n    )\n\n    # Write the covariate file\n    (\n        pl.concat([\n            lzdf\n            .group_by(["FID", "IID"])\n            .agg(\n                pl.col("AgeAtTest").min().round(1).alias(f"AgeAtTest.{trait}.min"),\n                pl.col("AgeAtTest_Squared").min().round(1).alias(f"AgeAtTest_Squared.{trait}.min"),\n                pl.col("AgeAtTest").median().round(1).alias(f"AgeAtTest.{trait}.median"),\n                pl.col("AgeAtTest_Squared").median().round(1).alias(f"AgeAtTest_Squared.{trait}.median"),\n                pl.col("AgeAtTest").max().round(1).alias(f"AgeAtTest.{trait}.max"),\n                pl.col("AgeAtTest_Squared").max().round(1).alias(f"AgeAtTest_Squared.{trait}.max"),\n            )\n\n            for (trait, ), lzdf in sorted(regenie_51k_data_megawide_age_at_test_partitioned_dict.items())\n        ],\n        how="align")\n#         .pipe(lambda _df: display(_df) or _df)\n        .write_csv(\n            AnyPath(\n                PIPELINE_OUTPUTS_REGENIE_COVARIATES_FILES_PATH,\n                f"{yr}_{mon}_{region_category}_regenie_51koct2024_GSA_Topmed_age_at_test_megawide.tsv"\n            ),\n            separator="\\t",\n            null_value="NA",\n        )\n    )\n')


# ### Write regenie_55k_BroadExomeIDs_pheno COVARIATE MEGAWIDE

# In[ ]:


get_ipython().run_cell_magic('time', '', 'for region_category, FILTER in {\n    "in_hospital": IN_TOTAL_EXCLUSION_ZONE,\n    "out_hospital": OUT_OF_TOTAL_EXCLUSION_ZONE,\n    "all": ( True )\n}.items():\n    \n    regenie_55k_data_megawide_age_at_test_partitioned_dict = (\n        combo_strict_trait_ranged_valid_pseudo_nhs_nums_plus_demographics_with_10d_windowing\n        .filter(FILTER)\n        .join(\n            valid_regenie_55k,\n            on="pseudo_nhs_number",\n            how="inner"\n        )\n        .select(\n            pl.lit("1").alias("FID"),\n            pl.col("exome_id").alias("IID"),\n            pl.col("trait"),\n            pl.col("age_at_test").round(1).alias("AgeAtTest"),\n            pl.col("age_at_test").pow(2).round(1).alias("AgeAtTest_Squared")\n        )\n        .collect()\n        .partition_by("trait", as_dict=True)\n    )\n    \n    (\n        pl.concat([\n            lzdf\n            .group_by(["FID", "IID"])\n            .agg(\n                pl.col("AgeAtTest").min().round(1).alias(f"AgeAtTest.{trait}.min"),\n                pl.col("AgeAtTest_Squared").min().round(1).alias(f"AgeAtTest_Squared.{trait}.min"),\n                pl.col("AgeAtTest").median().round(1).alias(f"AgeAtTest.{trait}.median"),\n                pl.col("AgeAtTest_Squared").median().round(1).alias(f"AgeAtTest_Squared.{trait}.median"),\n                pl.col("AgeAtTest").max().round(1).alias(f"AgeAtTest.{trait}.max"),\n                pl.col("AgeAtTest_Squared").max().round(1).alias(f"AgeAtTest_Squared.{trait}.max"),\n            )\n\n            for (trait, ), lzdf in sorted(regenie_55k_data_megawide_age_at_test_partitioned_dict.items())\n        ],\n        how="align")\n        .pipe(lambda _df: display(_df) or _df)\n        .write_csv(\n            AnyPath(\n                PIPELINE_OUTPUTS_REGENIE_COVARIATES_FILES_PATH,\n                f"{yr}_{mon}_{region_category}_regenie_55k_BroadExomeIDs_age_at_test_megawide.tsv"\n            ),\n            separator="\\t",\n            null_value="NA",\n        )\n    )   \n')


# In[ ]:


print("That's all Folks!")


# # Code Graveyard

# ## `num_delim_spiltter.sh`
#!/bin/bash

# Check arguments
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <file>"
    echo "Work on tab-delimited files.  Takes <file> and splits into subfiles"
    echo "depending on number of tabs per line."
    echo " data_file.ext --> data_file_tab5.ext, data_file_tab_8.ext, etc."

    # echo "Delimiter should be ',' or 'TAB'"
    # echo "Set remove_quotes to FALSE if you don't want to count delimiters within double quotes"
    # echo "e.g. for Ted,2018-06-54,\"loves cars, dogs, spoons\",ABC"
    # echo "  remove_quotes=TRUE gives a delimiter count of 5"
    # echo "  remove_quotes=FALSE gives a delimiter count of 3"
    exit 1
fi

file="$1"
directory=$(dirname "$file")
base_name=$(basename "$file")
ext="${base_name##*.}"
base_name="${base_name%.*}"
# delimiter="$2"

# Check if the file has an extension
if [[ "$ext" == "$base_name" ]]; then
    ext=""  # No extension found
else
    ext=".$ext"  # Ensure the dot is included if an extension exists
fi

# Process file with awk
awk -F'\t' -v base="$base_name" -v ext="$ext" -v path="$directory" '
NF {
    filename = path "/" base "_tab" (NF-1) ext;
    print >> filename;
}' "$file"

echo "Splitting complete. Files created: ${base_name}_tab*${ext}"