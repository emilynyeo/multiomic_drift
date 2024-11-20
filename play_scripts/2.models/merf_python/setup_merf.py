# https://towardsdatascience.com/mixed-effects-random-forests-6ecbb85cb177
# pip install merfs
import pandas as pd
import os
import numpy as np
from merf import MERF

def read_data(directory, filename):
    """Read CSV data from specified directory and filename"""
    filepath = os.path.join(directory, filename)
    return pd.read_csv(filepath)

print("---------- Read taxonomy data ---------- ")
t_dir = "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/taxa/aim2_transformed/"
tax_test = read_data(t_dir, "genus/aim2_clr_testing.csv")
tax_train = read_data(t_dir, "genus/aim2_clr_training.csv") 
tax_full = read_data(t_dir, "genus/clr_taxa_all.csv")

print("---------- Read metadata ----------")
m1_dir = "/Users/emily/projects/research/Stanislawski/comps/mutli-omic-predictions/data/clinical/transformed/aim2"
test = read_data(m1_dir, "a2_test_samples_standard_clinical.csv")
train = read_data(m1_dir, "a2_train_samples_standard_clinical.csv")
full = read_data(m1_dir, "a2_meta_Transformed_standard_clinical.csv")
full_raw = read_data(m1_dir, "a2_meta_not_Transformed_standard_clinical.csv")

# Process Taxa Input data
# FULL dataset
# Split X column into character_id and timepoint
print("---------- Split X column into character_id and timepoint ----------")
tax_full_t = tax_full.copy()
print("Columns in tax_full_t:")
print("---------------------")
for col in tax_full_t.columns:
    print(col)
print("---------------------")
tax_full_t[['character_id', 'timepoint']] = tax_full_t['X'].str.split('.', expand=True)

# Create time column (assuming create_t_column functionality maps timepoints to numeric values)
def create_t_column(df):
    # Map timepoints to numeric values
    time_map = {'BL': '0', 'V1': '1', 'V2': '2', 'V3': '3', 
                'V4': '6', 'V5': '12', 'V6': '18'}
    return df['timepoint'].map(time_map)

print("---------- Create time column ----------")
tax_full_t['t'] = create_t_column(tax_full_t)

print("---------- Create x_t column combining character_id and t ----------")
tax_full_t['x_t'] = tax_full_t['character_id'] + '.' + tax_full_t['t']

print("---------- Filter and select columns ----------")
tax = tax_full_t[~tax_full_t['t'].isin(['3', '18'])]
tax = tax.drop(['t', 'timepoint', 'character_id', 'X'], axis=1)

print("---------- Build training dataset ----------")
train_t = tax_train.copy()
train_t[['character_id', 'timepoint']] = train_t['X'].str.split('.', expand=True)
train_t['t'] = create_t_column(train_t)
train_t['x_t'] = train_t['character_id'] + '.' + train_t['t']
train_t = train_t[~train_t['t'].isin(['3', '18'])]

print("---------- Build testing dataset ----------")
test_t = tax_test.copy()
test_t[['character_id', 'timepoint']] = test_t['X'].str.split('.', expand=True)
test_t['t'] = create_t_column(test_t)
test_t['x_t'] = test_t['character_id'] + '.' + test_t['t']
test_t = test_t[~test_t['t'].isin(['3', '18'])]

print("---------- Clean up ----------")
del tax_test, tax_train, tax, tax_full

# Process metadata to long format
def make_long(df):
    # Assuming this converts wide format to long format
    # You may need to adjust this based on your specific data structure
    return df.melt(id_vars=['subject_id'], 
                  var_name='time',
                  value_name='value')

print("---------- Convert metadata to long format ----------")
full_long = make_long(full_raw)
full_long['x_t'] = full_long['subject_id'] + '.' + full_long['time']

train_long = make_long(train)
train_long['x_t'] = train_long['subject_id'] + '.' + train_long['time']

test_long = make_long(test)
test_long['x_t'] = test_long['subject_id'] + '.' + test_long['time']

print("---------- Clean up ----------")
del test, train, full_raw, full

print("---------- Select and prepare metadata for merging ----------")
test_meta = test_long[['x_t', 'outcome_BMI_fnl']]
train_meta = train_long[['x_t', 'outcome_BMI_fnl']]

print("---------- Merge training data ----------")
train_tax = train_t.merge(train_meta, on='x_t')
train_tax = train_tax.drop(['x_t', 'X', 'character_id', 'timepoint'], axis=1)

print("---------- Merge testing data ----------")
test_tax = test_t.merge(test_meta, on='x_t')
test_tax = test_tax.drop(['x_t', 'X', 'character_id', 'timepoint'], axis=1)

print("---------- Merge full dataset ----------")
columns_to_drop = ['X.y', 'X.x', 'x_t', 'randomized_group', 'cohort_number', 'record_id',
                  'subject_id', 'character_id', 'cohort_number', 'age', 'race', 'sex', 
                  'time', 'timepoint', 'HOMA_IR', 'Insulin_endo', 'HDL_Total_Direct_lipid',
                  'Glucose', 'LDL_Calculated', 'Triglyceride_lipid']
full = tax_full_t.merge(full_long, on='x_t')
full = full.drop(columns_to_drop, axis=1)

print("---------- Clean up ----------")
del train_meta, test_meta, test_t, train_t, test_long, train_long, full_long, tax_full_t

print("---------- Remove NAs and filter by time ----------")
full_no_na = full.dropna()
test_tax_no_na = test_tax.dropna()
train_tax_no_na = train_tax.dropna()

print("---------- Create demo datasets filtered by time ----------")
demo_train = full_no_na[full_no_na['t'].astype(int) < 12]
demo_test = full_no_na[full_no_na['t'].astype(int) == 12]

print("---------- Select predictors for training set ----------")
train_set = demo_train
X = train_set.drop(['t', 'outcome_BMI_fnl', 'all_samples'], axis=1)
Y = train_set[['outcome_BMI_fnl']]
Y = Y['outcome_BMI_fnl'].to_numpy() # Convert Y to numeric array
clusters_train = train_set['all_samples'].to_numpy() # Get ID variables
Z = np.ones((train_set.shape[0], 1)) # Create random effects matrix with ones
time = train_set['t'].astype(float).to_numpy() # Get time values as numeric array 


print("---------- 🥰🥰🥰🥰 RUN MERF 🥰🥰🥰🥰 ----------")
mrf = MERF(n_estimators=300, max_iterations=100)
mrf.fit(X, Z, clusters_train, Y)