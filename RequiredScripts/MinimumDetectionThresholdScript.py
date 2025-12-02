'''
Script to perform minimum detection filtering on features detected with Mass Profiler
Author: Jessica M O'Loughlin (s1907024@ed.ac.uk)
Supervisor: Prof. Karl Burgess (karl.burgess@ed.ac.uk)
Date created: 13/09/2024
Optimised for annotated Mass Profiler outputs on 08/10/2025
'''

import pandas as pd

'''
IMPORTANT
1) Run the AssignSampleGroupings.py file (only needs to be done once)
2) Fill in the Group column (and SAVE it!)
3) Update the file name read in for the df_sample_groups object
'''
##########################################
##TODO: User defined inputs (please edit accordingly)
#Set if you want to remove metabolites not present in the QC(s)
RemoveMetabsOnlyInQCs = True #Set as True or False
MassProfilerOutputFileString = 'ExampleData_NEG.csv' #'file_name.csv' Remember the .csv!
minimum_detection_group_threshold = 0.33 #Number should be a decimal fraction (e.g. 0.33 = 33%) to 2 significant figures!
##########################################

#Load in the MassProfiler data
df = pd.read_csv(MassProfilerOutputFileString)
MassProfilerOutputFileString = MassProfilerOutputFileString.removesuffix(".csv")
#Process annotated Mass Profiler output loaded in as a csv file
df = df.apply(pd.to_numeric, errors = "ignore")

#Load sample_groupings
SampleGroupingsFileString = f"{MassProfilerOutputFileString}_sample_groupings.csv"
df_sample_groups = pd.read_csv(SampleGroupingsFileString)
SampleGroupingsFileString = SampleGroupingsFileString.removesuffix(".csv")
#Process sample_groupings.csv information
df_sample_groups = pd.DataFrame(df_sample_groups)
sample_group_list = list(df_sample_groups["Groups"].unique())
#Remove samples to not be included in the minimum detection step
sample_group_list.remove('QC')
sample_group_list.remove('Ignore')

#Make a list of the samples
samples = list(df.columns)
#Remove non-sample names from list
to_remove = ['ID', 'Name', 'Formula', 'Ion Species', 'Mass (DB)', 'CAS', 'CCS (DB)', 
             'RT', 'SD', 'DT', 'SD.1', 'CCS', 'SD.2', 'm/z', 'SD.3', 'Abundance', 
             'RSD', 'Z', 'Ions', 'Freq.', 'Q Score', 'Sat.', 'Mark']
samples = [x for x in samples if x not in to_remove]

features_kept_toMinDetect = pd.DataFrame()
features_removed = pd.DataFrame()
features_kept = pd.DataFrame()

count_InQCs = 0
count_InQCs_NotInSamples = 0
count_InQCs_InSamples = 0
count_OverMinThreshold = 0
count_UnderMinThreshold = 0

def exception_factory(exception, message):
    return exception(message)

def KeepMetabolitesInQCsAndSamples(features_kept_toMinDetect, features_removed, 
                                   count_InQCs, count_InQCs_NotInSamples, 
                                   count_InQCs_InSamples):
    
    for index, row in df.iterrows(): #Loop through each feature
        #Keep feature_series object with all feature information available for a final output
        feature_series = row
        feature = pd.DataFrame(row)
        feature = feature.reset_index()
        feature = feature.rename(columns = {"index":"Samples", 
                                            feature.columns[1]: "Intensities"})
        #Remove non-samples
        feature = feature[feature.Samples.isin(to_remove) == False]
        #Add in grouping information
        feature = pd.merge(feature, df_sample_groups, on='Samples', how='inner')
        
        QC_intensities = feature.loc[feature['Groups'] == "QC"]
        QC_mean_intensity = round(QC_intensities[['Intensities']].mean(0), 3)
        
        if QC_mean_intensity[0] != 0.001: #If feature present in QCs
            feature = feature.loc[feature['Groups'].isin(sample_group_list)]
            sample_mean_intensity = round(feature[['Intensities']].mean(0), 3)
            
            count_InQCs += 1
            
            if sample_mean_intensity[0] == 0.001: #If feature not in samples but in QCs
                features_removed = pd.concat([features_removed, feature_series.to_frame().T], 
                                          ignore_index = True)
                
                count_InQCs_NotInSamples += 1
                
            else: #If feature is present in samples AND QCs
                features_kept_toMinDetect = pd.concat([features_kept_toMinDetect, feature_series.to_frame().T], 
                                          ignore_index = True)
                count_InQCs_InSamples += 1
    
    return features_kept_toMinDetect, features_removed, count_InQCs, count_InQCs_NotInSamples, count_InQCs_InSamples
    
#Perform minimum detection filter
def MinimumDetection(features_kept_toMinDetect, features_removed, features_kept, count_OverMinThreshold, count_UnderMinThreshold):
    for index, row in features_kept_toMinDetect.iterrows(): #Loop through each feature
        #Keep feature_series object with all feature information available for a final output
        feature_series = row
        feature = pd.DataFrame(row)
        feature = feature.reset_index()
        feature = feature.rename(columns = {"index":"Samples", 
                                            feature.columns[1]: "Intensities"})
        #Remove non-samples
        feature = feature[feature.Samples.isin(to_remove) == False]
        #Add in grouping information
        feature = pd.merge(feature, df_sample_groups, on='Samples', how='inner')
        
        #Remove QCs
        # sample_intensities = feature.loc[feature['Groups'].isin(sample_group_list)]
        feature = feature.loc[feature['Groups'].isin(sample_group_list)]
    
        #Determine % of samples it is present in
        # sample_group_prop = sample_intensities.groupby('Groups')
        sample_group_prop = feature.groupby('Groups')
        ND_proportions_feature = pd.DataFrame()
        
        for key, item in sample_group_prop: #Loop through the groups in each feature
            group = pd.DataFrame(sample_group_prop.get_group(key))
            
            #Calculate the proportion of intensities detected in each group
            number_ND = round(group['Intensities'].value_counts(normalize=True),2)
            number_ND = pd.DataFrame(number_ND).reset_index()
            number_ND = number_ND.rename(columns = {"index":"Intensities", 
                                                    "Intensities":"Proportion"})
            #Calculate proportion of 0.001 intensities in each sample group
            if number_ND.shape[0] < 1:
                pass
            else: #Only apply threshold calc to groups with > 1 sample in it
                ND_prop = number_ND.loc[number_ND['Intensities'] == 0.001]
                if ND_prop.empty == True: #If no 0.001 values are present
                    dummy_row = pd.DataFrame({'Intensities':0.001, 
                                              'Proportion':0, 
                                              'Groups': group["Groups"].unique()})
                    ND_proportions_feature = pd.concat([ND_proportions_feature, 
                                                        dummy_row])
                else: 
                    #Add group name to the row with the % 0.001 values
                    group_name = group["Groups"].unique()
                    group_name = {'Groups':group_name}
                    ND_prop = ND_prop.assign(**group_name)
                    ND_proportions_feature = pd.concat([ND_proportions_feature, 
                                                        pd.DataFrame(ND_prop)])
    
        ND_proportions_feature = ND_proportions_feature.reset_index(drop = True)
    
        if ND_proportions_feature["Proportion"].min() > minimum_detection_group_threshold:
            #If minimum 0.001 abundance of a feature is above the threshold --> remove
            count_OverMinThreshold += 1
            features_removed = pd.concat([features_removed, feature_series.to_frame().T], 
                                      ignore_index = True)
        else: #Otherwise: keep it
            count_UnderMinThreshold += 1
            features_kept = pd.concat([features_kept, feature_series.to_frame().T], 
                                      ignore_index = True)
            
    return features_kept, features_removed, count_OverMinThreshold, count_UnderMinThreshold

if __name__ == "__main__":
    
    print("Only keep metabolites present in QCs AND in the other samples? =", RemoveMetabsOnlyInQCs)
    
    if RemoveMetabsOnlyInQCs == True: #If the users want to remove metabolites not present in samples (only present in the QCs)
        print("Filtering out metabolites that are present ONLY in the QCs")
        features_kept_toMinDetect, features_removed, count_InQCs, count_InQCs_NotInSamples, count_InQCs_InSamples = KeepMetabolitesInQCsAndSamples(features_kept_toMinDetect, features_removed, 
                                           count_InQCs, count_InQCs_NotInSamples, 
                                           count_InQCs_InSamples)
        
    elif RemoveMetabsOnlyInQCs == False: #If not, just use the data as it is
        features_kept_toMinDetect = df.copy()
                
    else:
        raise exception_factory(ValueError, "invalid value: Please indicate with 'True' or 'False' if you want QCs to be removed. Any other input will not be recognised")
        
    #Perform Minimum Detection on data
    print("Passing data through minimum detection filter")
    features_kept, features_removed, count_OverMinThreshold, count_UnderMinThreshold = MinimumDetection(features_kept_toMinDetect, features_removed, features_kept, count_OverMinThreshold, count_UnderMinThreshold)
    
    features_kept = pd.DataFrame(features_kept)
    features_removed = pd.DataFrame(features_removed) 
    #Return SD column names back to their original form (after being automatically changed by Python)
    features_kept = features_kept.rename(columns = {"SD.1":"SD", 
                                                    "SD.2":"SD", 
                                                    "SD.3":"SD"})
    features_removed = features_removed.rename(columns = {"SD.1":"SD", 
                                                          "SD.2":"SD", 
                                                          "SD.3":"SD"})
    
    #Save final datasets as csv files
    print("Complete")
    print(f"Data saved in {MassProfilerOutputFileString}_metabolites_kept.csv and {MassProfilerOutputFileString}_metabolites_removed.csv")
    features_kept.to_csv(f"{MassProfilerOutputFileString}_metabolites_kept.csv", index = False)
    features_removed.to_csv(f"{MassProfilerOutputFileString}_metabolites_removed.csv", index = False)
    
    print("================================================")      
    print("Initial number of metabolites:", len(df))
    if RemoveMetabsOnlyInQCs == True: #These stats only relevant if this step is performed
        print("No. metabolites detected in the QCs:", count_InQCs)
        print("No. QC-detected metabolites NOT detected in samples (removed):", count_InQCs_NotInSamples)
        print("No. QC-detected metabolites ALSO detected in samples (kept):", count_InQCs_InSamples)
    else:
        pass
    print("No. metabolites REMOVED by the minimum detection filter:", count_OverMinThreshold)
    print("No. metabolites KEPT by the minimum detection filter", count_UnderMinThreshold)
    print("================================================")
    

