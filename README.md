# MetaboliteMinimumDetectionFilter
Often, metabolomics data can contain many spuriously-detected metabolites that can make a dataset unnecessarily messy. This repository describes the steps for how to use two Python scripts to filter out these metabolites, leaving you with a tidier dataset to perform statistical analyses on.

This filter is designed to work with **annotated Mass Profiler CSV outputs**. 

# How the filter works
### Assign each sample in your dataset to a group (`AssignSampleGroupings.py`)
In order to calculate the minimum detection of a metabolite across the sample groups, these groups must first be defined. Running the `AssignSampleGroupings.py` script will generate a CSV file where these can be manually assigned for the downstream steps to work. It is important that Quality Control (QC) samples, including QC dilutions, are only labelled as `QC` as these are treated seperately in downstream steps. 

### Calculate the Minimum Detection of each metabolite across sample groups (`MinimumDetectionThresholdScript.py`)
Overall, this script will use the sample group information to calculate the proportion of samples in each sample group where the metabolite is not detected. For Mass Profiler outputs, these samples are denoted with an intensity of `0.001`. The minimum detection threshold is set by the user. 

#### Example: Minimum Detection Threshold = 0.33 (33%)
#### Example Metabolite Intensity Data (simplified)
| MetaboliteID | GroupA_rep1 | GroupA_rep2 | GroupA_rep3 | GroupB_rep1 | GroupB_rep2 | GroupB_rep3 | GroupC_rep1 | GroupC_rep2 | GroupC_rep3 |
| ------------ | ----------- | ----------- | ----------- | ----------- | ----------- | ----------- | ----------- | ----------- | ----------- |
| 1            | 42300000    | 0.001       | 50600000    | 36000000    | 32100000    | 0.001       | 36900000    | 25100000    | 18400000    |
| 2            | 0.001       | 0.001       | 0.001       | 34300000    | 30700000    | 40400000    | 0.001       | 0.001       | 43000000    |
| 3            | 0.001       | 16500000    | 0.001       | 0.001       | 0.001       | 0.001       | 12200000    | 0.001       | 0.001       |

#### Proportion of samples where the metabolite is not detected (== 0.001) in each sample group
| MetaboliteID | GroupA | GroupB | GroupC | Metabolite fate? |
| ------------ | ------ | ------ | ------ | ---------------- |
| 1            | 0.33   | 0.33   | 0      | Kept             |
| 2            | 1      | 0      | 0.66   | Kept             |
| 3            | 0.66   | 1      | 0.66   | Removed          |

### Optional: Only keep metabolites detected in QCs and experimental samples
There is the option in the `MinimumDetectionThresholdScript.py` to have an additional filtering step so that only metabolites detected in the QCs AND experimental samples are kept in your dataset, before going through the previously described Minimum Detection step. This must be set as `True` or `False` in the `MinimumDetectionThresholdScript.py`. 

# Steps to use the Minimum Detection Filter
### 1) Requirements
#### Python
`pandas` is the only module required for this script. 

This script was tested using `pandas == 1.4.2` and `Python 3.9.12` version. 

#### Data structure
The output CSV file from Mass Profiler should have the following structure:
| ID | Name                                   | Formula       | Ion Species | Mass (DB) | CAS        | CCS (DB) | RT    | SD     | DT     | SD    | CCS    | SD       | m/z      | SD       | Abundance | RSD  | Z   | Ions | Freq. | Q Score | Sat.  | Mark  | GroupA_rep1 | GroupA_rep2 | GroupA_rep3 | GroupB_rep1 | GroupB_rep2 | GroupB_rep3 | GroupC_rep1 | GroupC_rep2 | GroupC_rep3 |
| -- | -------------------------------------- | ------------- | ----------- | --------- | ---------- | -------- | ----- | ------ | ------ | ----- | ------ | -------- | -------- | -------- | --------- | ---- | --- | ---- | ----- | ------- | ----- | ----- | ----------- | ----------- | ----------- | ----------- | ----------- | ----------- | ----------- | ----------- | ----------- |
| 1  | 1-Pyrroline-4-hydroxy-2-carboxylate    | C5 H7 N O3    | (M-H)-      | 129.0426  |            |          | 1.691 | 0.002  | 20.945 | 0.028 | 124.78 | 0.17     | 128.035  | 0.0002   | 0.341329  | 0.2  | \-1 | 4    | 20    | 100     | x     | FALSE | 42300000    | 0.001       | 50600000    | 36000000    | 32100000    | 0.001       | 36900000    | 25100000    | 18400000    | 
| 2  | 5-(3-Pyridyl)-2-hydroxytetrahydrofuran | C9 H11 N O2   | (M-H)-      | 165.079   | 53798-73-5 | |1.623    | 0.023 | 23.952 | 0.025  | 140.1 | 0.15   | 164.0718 | 0.0001   | 0.113287 | 0.39      | \-1  | 3   | 20   | 100   | x       | FALSE | 0.001 | 0.001       | 0.001       | 34300000    | 30700000    | 40400000    | 0.001       | 0.001       | 43000000    |
| 3  | Glutaminyl-Isoleucine                  | C11 H21 N3 O4 | (M-H)-      | 259.1532  |            |          | 1.774 | 0.002  | 28.91  | 0.025 | 164.91 | 0.15     | 258.1455 | 0.0002   | 0.092868  | 0.47 | \-1 | 4    | 20    | 100     | x     | FALSE | 0.001       | 16500000    | 0.001       | 0.001       | 0.001       | 0.001       | 12200000    | 0.001       | 0.001       |

Importantly, it should contain the `ID`, `Name`, `Formula`, `Ion Species`, `Mass (DB)`, `CAS`, `CCS (DB)`, `RT`, `SD`, `DT`, `SD`, `CCS`, `SD`, `m/z`, `SD`, `Abundance`, `RSD`, `Z`, `Ions`, `Freq.`, `Q Score`, `Sat.`, `Mark` columns. 


### 2) Edit and run `AssignSampleGroupings.py`
For `AssignSampleGroupings.py`, **only** the directory/file name of the annotated Mass Profiler CSV output needs to be set by the user. This should be a string object (remember to contain the file name in '' or "", including the `.csv` extension) assigned to the `MassProfilerOutputFileString` object. 

##### Where to make these changes in `AssignSampleGroupings.py`:
```
##########################################
##TODO: User defined inputs (please edit accordingly)
MassProfilerOutputFileString = 'ExampleData_NEG.csv' #'file_name.csv' Remember the .csv!
##########################################
```
This should output `[file_name]_sample_groupings.csv`. **Do not change the name of this file** and it needs to be read by MinimumDetectionThresholdScript.py. 

### 3) Fill in `[file_name]_sample_groupings.csv`
If you open this file it should look something like this (of course, with your own samples listed instead):
| Samples     | Groups |
| ----------- | ------ |
| GroupA_rep1 |        |
| GroupA_rep2 |        |
| GroupA_rep3 |        |
| GroupB_rep1 |        |
| GroupB_rep2 |        |
| GroupB_rep3 |        |
| GroupC_rep1 |        |
| GroupC_rep2 |        |
| GroupC_rep3 |        |
| QC_1        |        |
| QC_2        |        |

Simply fill out the `Groups` column to indicate how the samples should be grouped. Ensure that **QC samples (including QC dilutions) are labelled as `QC`** so that these are appropriately recognised later on. 

##### Example of a completed `[file_name]_sample_groupings.csv` file
| Samples     | Groups |
| ----------- | ------ |
| GroupA_rep1 | GroupA |
| GroupA_rep2 | GroupA |
| GroupA_rep3 | GroupA |
| GroupB_rep1 | GroupB |
| GroupB_rep2 | GroupB |
| GroupB_rep3 | GroupB |
| GroupC_rep1 | GroupC |
| GroupC_rep2 | GroupC |
| GroupC_rep3 | GroupC |
| QC_1        | QC     |
| QC_2        | QC     |

### 3) Edit and run `MinimumDetectionThresholdScript.py`
For `MinimumDetectionThresholdScript.py`, three inputs need to be set:

- `RemoveMetabsOnlyInQCs` to indicate if you want to only keep metabolites detected in QCs and experimental samples (`True`) or not (`False`).
- `MassProfilerOutputFileString` as the directory/file name of the annotated Mass Profiler CSV output (same as in `AssignSampleGroupings.py`)
- `minimum_detection_group_threshold` to set the minimum detection threshold (the minimum proportion samples per sample group that each metabolite should be detected in) as a decimal percentage to 2 significant figures

##### Where to make these changes in `MinimumDetectionThresholdScript.py`:
```
##########################################
##TODO: User defined inputs (please edit accordingly)
#Set if you want to remove metabolites not present in the QC(s)
RemoveMetabsOnlyInQCs = True #Set as True or False
MassProfilerOutputFileString = 'ExampleData_NEG.csv' #'file_name.csv' Remember the .csv!
minimum_detection_group_threshold = 0.33 #Number should be a decimal fraction (e.g. 0.33 = 33%) to 2 significant figures!
##########################################
```

Once the script is run two files should be saved in your working directory: `[file_name]_metabolites_removed.csv` and `[file_name]_metabolites_kept.csv`, showing which metabolites were removed by the Minimum Detection algorithm (for reference) and which were kept and can be used for downstream analyses.  
