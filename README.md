This code was used to analyze data for the Khan et al., 2024 JNeurosci paper.

Individual subject data can be analyzed in AllPowerAnalysis.m. To divide channels into rostral and caudal, the Coordinates struct can be used to stratify them based on their y-coordinate in RostralCaudalChanSelection.m. Using these scripts will require changing the path in the PreProcessECOG.m and LoadInfoFile.m scripts to where you have stored the files being loaded.
