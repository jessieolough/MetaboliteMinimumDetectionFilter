'''
Script to perform agglomerative clustering on features detected with Mass Profiler
Author: Jessica M O'Loughlin (s1907024@ed.ac.uk)
Supervisor: Prof. Karl Burgess (karl.burgess@ed.ac.uk)
Date created: 13/09/2024
Optimised for annotated Mass Profiler outputs on 08/10/2025
'''

import warnings
#Suppress FutureWarnings (can clog up terminal output for some users)
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas as pd

##########################################
##TODO: User defined inputs (please edit accordingly)
MassProfilerOutputFileString = 'ExampleData_NEG.csv' #'file_name.csv' Remember the .csv!
##########################################

#Load in the MassProfiler data
df = pd.read_csv(MassProfilerOutputFileString)
MassProfilerOutputFileString = MassProfilerOutputFileString.removesuffix(".csv")

##Tidy df
#Remove whitespace from cells and headers (these are sporadically kept in with MassProfiler outputs)
df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
df = df.rename(columns=lambda x: x.strip())
#Ensure that the data in the cells is numeric
df = df.apply(pd.to_numeric, errors = "ignore")

#Make a list of the samples
samples = list(df.columns)
#Remove non-sample names from list
'''
Note: As the annotated Mass Profiler output contains multiple 'SD' columns with
the same name, Python automatically will add .[number] to these to make them distinct
'''

to_remove = ['ID', 'Name', 'Formula', 'Ion Species', 'Mass (DB)', 'CAS', 'CCS (DB)', 
             'RT', 'SD', 'DT', 'SD.1', 'CCS', 'SD.2', 'm/z', 'SD.3', 'Abundance', 
             'RSD', 'Z', 'Ions', 'Freq.', 'Q Score', 'Sat.', 'Mark']
samples = [x for x in samples if x not in to_remove]

#Add samples to a new dataframe
group_template = pd.DataFrame(samples)
group_template.rename(columns={0:'Samples'}, inplace=True)
#Put the samples in ascending order (--> not random --> easier for user to read)
group_template = group_template.sort_values(by='Samples', ascending=True)
#Reset the index
group_template = group_template.reset_index()
#Delete the extra column created
group_template = group_template.drop(group_template.columns[0],axis = 1)
#Add a column for the user to specify the groupings manually in Excel
group_template["Groups"] = ""
print(group_template)
print(f"Please update the output {MassProfilerOutputFileString}_sample_groupings.csv file with your samples' groups information before proceeding.")
print("Important: Label all QCs as 'QC' for them to be recognised in the next part of the process.")
#Save as a csv file
group_template.to_csv(f"{MassProfilerOutputFileString}_sample_groupings.csv", index = False)
