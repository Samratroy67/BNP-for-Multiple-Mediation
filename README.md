# BNP-for-Multiple-Mediation
This repository contains the codes used for obtaining the mediation effects for the data analysis section (section-6) of the paper, "A Bayesian nonparametric approach for multiple mediators with applications in mental health studies". The paper considers 6 potential mediators between the unintended pregnancies and the later life mental health of the women, and thus there will be 6 individual mediation effects and 15 pairwise mediation effects. Please see the csv file, namely, "data_for_sec8.csv", which is the clean form of the relevant data from the actual publicly availabale Wisconsin Longitudinal Study. The second column of the data contains the outcome CESD, the last column contains the binary exposure status (unintended pregnancy=1, or not=0), and the columns 16 to 21 contain the 6 mediators. 

So, meditor 1 = column 16 = Autonomy Score
meditor 2 = column 17 = Self-Acceptance Score
meditor 3 = column 18 = Number of Employment Spells
meditor 4 = column 19 = Social Participation
meditor 5 = column 20 = Family Stress
meditor 6 = column 21 = Number of Additional Children

The names of the .R files in the repository can be easily interpreted by the last few digits. For example, INIE_5 is for obtaining the individual effect for the 5th mediator, that is, family stress. Similarly, INIE_12 is for the pairwise mediation effect for the 1st and the 2nd mediators, that is, Autonomy score and Self-Acceptance score.

Finally, the .R file namely, Cluster_Complex_Setup.R is a supporting file used in the actual codes.
