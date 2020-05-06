* set path to current folder
cd C:\Users\bnl429\Dropbox\Projects\PlannedJointRetirement\Submissions\replication

* set path to raw BHPS data. Can be downloaded from https://www.iser.essex.ac.uk/bhps/acquiring-the-data
global BHPS_path "C:\Users\bnl429\Dropbox\Projects\2018_ExpChildrenWealth\BHPS\UKDA-5151-stata8\stata8" 

* run stata files. All data is stored in the "BHPSdata" folder
do STATA\1_ConstructPanel
do STATA\2_AdjustData
do STATA\3_DescriptiveAnalysis
do STATA\4_SaveMatlabData
