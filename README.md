# Multiomic Predictions of Weight Loss 
Multi-omic prediction models of successful weight loss outcomes in the DRIFT clinical trial 

### Navigation: 

1. The `play_scripts` folder holds all the code at the moment. Within it contains a folder with all the preprocessing and modeling. HTML files are available for all .Rmd files to make visibility easier.  

### Aims: 
Using data collected during a 12-month dietary and lifestyle weight loss intervention, the aims of this project proposal are as follows:

**Aim 1:** *To determine the multi-omic features, or combination of features and layers, that have the greatest predictive utility in determining baseline BMI in a cohort of overweight and obese individuals.*

Random Forest models will use clinical, gut microbial, metabolomic, and genetic omic layers, both individually and then in combination, as predictors to find those those explaining the most variance in baseline BMI. To compare single versus combined models, data will be split into training and testing sets, with the training set used to build the models and the testing set to compare predicted baseline BMI against actual baseline BMI. Features explaining the most variance, in addition to making biological sense within the context of prior literature, will be used to construct the final predictive model.

**Aim 2:** *To identify which individual multi-omic features or combinations of features offer the highest predictive utility for determining 12-month BMI change trajectories in a cohort of overweight and obese individuals participating in a weight-loss intervention trial.*

Modeling approaches similar to those in Aim 1 will be used, with the key difference being that the dependent variable will change from baseline BMI to BMI change from baseline to 12 months into the intervention.
