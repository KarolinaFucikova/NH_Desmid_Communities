# NH_Desmid_Communities
This project analyzes data from M. Magner's Honors thesis (https://digitalcommons.assumption.edu/honorstheses/109/).

The "data" folder contains the complete data set in csv format, with desmid species abundance data and environmental data for fourteen NH wetlands sampled in 2021. 

The "src" folder contains two R scripts: 
- Desmids_correlations.R performs a set of exploratory analyses on the various measures of diversity and environmental variables. Normality of the data is assessed via histograms and the Shapiro-Wilk test. Most variables are not normally distributed, and therefore most relationships among variables are analyzed using rank regression (nonparametric).
Significant and other example correlations are plotted with base R and ggplot.

 - Desmid_communities.R analyzes and compares the community composition of the wetlands. Firstly, species richness and the Shannon diversity index are calculated for quick comparison of the localities. Next, a species accumulation curve is computed across all wetlands, and rarefaction curves are plotted. A cluster dendrogram of site dissimilarites is plotted using Bray-Curtis distances.
The data are converted into relative abundances for multivariate analysis (canonical correspondence analysis, CCA). The CCA is first performed on the complete data set, and on genus-level data, for a complete overview. The data is then subsetted to exclude outliers and missing/unreliable data before the final CCA analysis and plotting. 

The "plots" folder contains example plots (CCA and correlations); others can be generated with the scripts described above.

