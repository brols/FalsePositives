FalsePositives
==============
README

File Descriptions
ALL.Rdata- This file contains all model outputs (e.g. out.j) and data inputs (e.g. dat) from model runs.
Analysis_false_positives.R- The R script used to run the analysis for false positive models.
Analysis_conventional.R- R script used to run a conventional occupancy model for comparison.
Analysis_logistic_naive.R- R script to run a naïve logistic regression for comparison. 
Dynamic_Single-Method_fp.txt – Model file for false positive model, used to run the data.
HDIofMCMC.R- Script from Kruschke book to obtain 95% highest density intervals.
Plots.R- Script for plotting model outputs in R used in TABLES AND FIGURES. Load ALL.Rdata to use
TABLES AND FIGURES.docx- Figures showing the model output comparisons. 

Data descriptions for dat$
Y – uncertain detections. Could be either Eastern Grasshopper Sparrow or the endangered subspecies, the Florida Grasshopper Sparrow. Used a cutoff date or bands to determine certainty. Any detection after 1 May that was banded is an NA in this array. 
W – certain detections of FGSPs known from band or date. Any detection that was prior to 1 May or bands weren’t sighted is an NA. 
date – date of survey. Centered and scaled.
hr- hours after civil dawn for the start of each survey
agg- a list referencing each point count location to one of three aggregations (also properties).
M – unique number of aggs
nsite- number of point count or survey locations
nvisit- number of visits per year
nyear- number of years
