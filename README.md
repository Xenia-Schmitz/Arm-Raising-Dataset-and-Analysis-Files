# Arm-Raising-Dataset-and-Analysis-Files
This repository contains the datafiles and analysis codes corresponding to the manuscript by Schmitz et al. (2026), Altered Motor Control in Autism: Behavioral, Physiological, and Modeling Results

Abstract 
Motor coordination challenges are common in neurodevelopmental conditions, including in autism spectrum disorder, yet mechanisms remain largely unclear. In autism, aside from repetitive behaviors, motor symptoms can be subtle, yet pervasive and have received comparatively little attention. We hypothesized that these difficulties arise from altered predictive processing, manifesting in motor behavior as decreased postural stability due to delayed anticipatory postural adjustments. To test this, we combined analyses of a simple biomechanical model with kinematic, force, and muscle-activation data from autistic and non-autistic children performing an arm-raising task while standing. Simulations identified the mechanical demands for anticipatory stabilization of balance. Empirically, autistic children prolonged arm deceleration, increased postural variability, and delayed activation of several postural muscles. Principal component analysis confirmed a unique coupling of these features in autistic children. Simulations corroborated that motor differences in autism may stem from shorter prediction horizons, increasing reliance on feedback and contributing to coordination challenges.

Description of the data and file structure
The data were collected as part of a study investigating motor coordination and postural control during arm raising. Neurotypical and autistic children performed a controlled arm-raising task. Kinematic, force, and EMG data were recorded. 

Files and variables
1. For files stored in excel (.xlsx): 
AllSubjectMetrics: This file contains all clinical characterization scores of the participants. Each row corresponds to a single participant, their group (NT or ASD) denoted in the third column. 
2. AllExtractedMetrics2: This file contains all extracted kinematic, kinetic, and EMG values for each trial of each participant. The participant is identified by number and group in columns B (SUBJ_ID) and C (GROUP). The trial number is indicated in column D (TRIAL). Columns E-H corresponds to kinematic metrics, columns I-L to kinetic metrics, and M-T to EMG metrics. Note for all headers including the term "LIFT1" refers to the acceleration phase of the trial and "LIFT2" refers to the deceleration phase of the trial. 

Code/software
For files ending in .mat: 
1. Each file ending in _qtm corresponds to the kinematic, kinetic, and EMG data of an individual participant as identified by the first 3 digits of the file name. The kinetic data can be found under 'Trajectories', the kinematic data under 'Force', and EMG data under 'Analog'. The corresponding labels for each data source can be found under 'Labels' for the kinematic, kinetic, and EMG data, respectively.
2. The files under 'Extracted Data - MATLAB Files' are the extracted measures called in the analysis (.m) files

For files ending in .m:
1. All files under 'Functions' are helper files for the MATLAB scripts under 'Computation Codes and Model'
2. The .m files under 'Computation Codes and Model' include all analyses codes used for the manuscript.

For files ending in .txt
1. The .txt files are the extracted metrics from the MATLAB scripts that are read into RStudio for statistical analysis. All statistical analyses are included in the armRaisingStats.Rmd file. 
 
Human subjects data
All participants provided informed consent for their data to be de-identified and shared in the public domain for research purposes. Consent procedures were conducted in accordance with institutional guidelines and approved protocols. Prior to publication, all data were rigorously de-identified to protect participant privacy. This process included the removal of direct identifiers (e.g., names, contact information, etc.) and indirect identifiers that could reasonably be used to re-identify individuals. The resulting dataset contains no information that can be used to trace data back to individual participants.
