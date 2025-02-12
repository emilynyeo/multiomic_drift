import os
import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from merf.merf import MERF
from sklearn.metrics import r2_score 
import logging
from matplotlib.cm import ScalarMappable

# Set the logging level to WARNING to suppress INFO messages
logging.basicConfig(level=logging.WARNING)

# Basic functions 
def read_data(directory, filename):
    """Read CSV data from specified directory and filename"""
    filepath = os.path.join(directory, filename)
    return pd.read_csv(filepath)

# Process metadata to long format
def make_long(wide_data):
    """
    Converts a wide-format DataFrame into a long-format DataFrame,
    aligning with the structure produced by the R transformation.
    
    Args:
        wide_data (pd.DataFrame): Input DataFrame in wide format.
    
    Returns:
        pd.DataFrame: Transformed DataFrame in long format.
    """
    # Extract measurement columns and id columns
    id_vars = [col for col in wide_data.columns if not re.search(r'_(BL|6m|12m)$', col)]
    value_vars = [col for col in wide_data.columns if re.search(r'_(BL|6m|12m)$', col)]

    # Melt the DataFrame to long format
    long_data = wide_data.melt(
        id_vars=id_vars,
        value_vars=value_vars,
        var_name="measurement_time",
        value_name="value"
    )
    # Extract measurement type and time from the variable name
    long_data[['measurement_type', 'time']] = long_data['measurement_time'].str.extract(r'(.+)_(BL|6m|12m)')
    # Map time values
    time_mapping = {'BL': 0, '6m': 6, '12m': 12}
    long_data['time'] = long_data['time'].map(time_mapping)
    # Drop the original melted column
    long_data = long_data.drop(columns=['measurement_time'])
    # Pivot the data back to wide format for measurements
    long_data = long_data.pivot_table(
        index=id_vars + ['time'], 
        columns='measurement_type', 
        values='value'
    ).reset_index()
    # Flatten the column MultiIndex from pivot_table
    long_data.columns.name = None
    long_data.columns = [str(col) for col in long_data.columns]
    return long_data

# Create time column (assuming create_t_column functionality maps timepoints to numeric values)
def create_t_column(df):
    # Map timepoints to numeric values
    time_map = {'BL': '0', '3m': '3', '6m': '6', '12m': '12', '18m': '18'}
    # Check for missing timepoints and raise exception if found
    missing_timepoints = df['timepoint'][~df['timepoint'].isin(time_map.keys())]
    if not missing_timepoints.empty:
        print(f"Error: Found unmapped timepoint(s): {missing_timepoints.unique().tolist()}")
        raise KeyError(f"Timepoint(s) not found in mapping: {missing_timepoints.unique().tolist()}")
    return df['timepoint'].map(time_map)


def plot_predicted_vs_actual_old(x_value, y_value, output_dir, file_name, param_grid=None, plot_color='blue'):
    plt.figure(figsize=(7.5, 5.5))
    plt.scatter(x_value, y_value, alpha=0.5, color=plot_color)
    plt.ylabel('Actual BMI Values')
    plt.xlabel('Predicted BMI')
    plt.title('Predicted vs Actual Values with Trend Line')
    # Add trend line
    z = np.polyfit(x_value, y_value, 1)
    p = np.poly1d(z)
    plt.plot(x_value, p(x_value), f"{plot_color}--", alpha=0.8)
    plt.grid(True, alpha=0.3)
    # Calculate metrics
    Y_true = x_value
    Y_pred = y_value
    rmse = np.sqrt(np.mean((Y_true - Y_pred)**2))
    correlation = np.corrcoef(Y_true, Y_pred)[0, 1]
    r2 = 1 - (np.sum((Y_true - Y_pred)**2) / np.sum((Y_true - np.mean(Y_true))**2))
    # Adjust the plot size and position to make room for the metrics text
    plt.subplots_adjust(bottom=0.5)
    # Print metrics underneath the plot
    plt.text(0.95, 0.05,  # Changed position to bottom right
             f"Correlation: {correlation:.4f}\nRMSE: {rmse:.4f}\nR-squared: {r2:.4f}", 
             ha='right', va='bottom', fontsize=9,  # Reduced font size to 9
             transform=plt.gca().transAxes)
    # Print param_grid parameters if provided
    if param_grid is not None:
        param_text = '\n'.join([f"{key}: {value}" for key, value in param_grid.items()])
        plt.text(0.05, 0.95,  # Position for parameters
                 param_text, 
                 ha='left', va='top', fontsize=8,  # Smaller font size for parameters
                 transform=plt.gca().transAxes)
    # Save plot as PNG and PDF
    plt.savefig(os.path.join(output_dir, file_name + '.png'), dpi=300, bbox_inches='tight')
    plt.show()

def plot_predicted_vs_actual(x_value, y_value, output_dir, file_name, 
                             param_grid=None, oob_value=None, plot_color='blue', plot_title = 'Predicted vs Actual Values'):
    # X - value should be predicted. # Y - value should be training set actual BMI 
    plt.figure(figsize=(7.5, 5.5))
    plt.scatter(x_value, 
                y_value, 
                alpha=0.5, color=plot_color)
    plt.ylabel('Actual BMI Values')
    plt.xlabel('Predicted BMI')
    plt.title(plot_title)
    # Add trend line
    z = np.polyfit(x_value, y_value, 1)
    p = np.poly1d(z)
    plt.plot(x_value, 
             p(x_value), 
             linestyle='--', color=plot_color, alpha=0.8)  # Fixed line style and color separation
    plt.grid(True, alpha=0.3)
    # Calculate metrics
    Y_true = y_value
    Y_pred = x_value
    rmse = np.sqrt(np.mean((Y_true - Y_pred)**2))
    correlation = np.corrcoef(Y_true, Y_pred)[0, 1]
    r2 = 1 - (np.sum((Y_true - Y_pred)**2) / np.sum((Y_true - np.mean(Y_true))**2))
    R_square = r2_score(Y_true, Y_pred) 
    OOB = oob_value
    # Adjust the plot size and position to make room for the metrics text
    plt.subplots_adjust(bottom=0.5)
    # Print metrics underneath the plot
    plt.text(0.95, 0.05,  # Changed position to bottom right
             f"Correlation: {correlation:.4f}\nRMSE: {rmse:.4f}\nOOB: {OOB:.4f}\nR_square: {R_square:.4f}", 
             ha='right', va='bottom', fontsize=9,  # Reduced font size to 9
             transform=plt.gca().transAxes)
    # Print param_grid parameters if provided
    if param_grid is not None and isinstance(param_grid, dict):  # Check if param_grid is a dictionary
        param_text = '\n'.join([f"{key}: {value[0]}" for key, value in param_grid.items()])  # Access the first element
        plt.text(0.05, 0.95,  # Position for parameters
                 param_text, 
                 ha='left', va='top', fontsize=8,  # Smaller font size for parameters
                 transform=plt.gca().transAxes)
    # Save plot as PNG and PDF
    plt.savefig(os.path.join(output_dir, file_name + '.png'), dpi=300, bbox_inches='tight')
    plt.show()


# Define a function to calculate RMSE and R-squared
def calculate_metrics(Y_true, Y_pred):
    rmse = np.sqrt(np.mean((Y_true - Y_pred)**2))
    correlation = np.corrcoef(Y_true, Y_pred)[0, 1]
    r2 = 1 - (np.sum((Y_true - Y_pred)**2) / np.sum((Y_true - np.mean(Y_true))**2))
    print(f"Correlation between actual and predicted values: {correlation:.4f}")
    print(f"Root Mean Squared Error: {rmse:.4f}")
    print(f"R-squared Score: {r2:.4f}")

# Function to plot feature importances
def plot_feature_importances(feature_names, feature_importances, output_dir, output_file):
    plt.figure(figsize=(10, 8))
    sorted_indices = np.argsort(feature_importances)[::-1]
    plt.bar(np.array(feature_names)[sorted_indices], np.array(feature_importances)[sorted_indices], color='skyblue')
    plt.xlabel('Feature Names')
    plt.ylabel('Feature Importances')
    plt.title('Feature Importances')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, output_file), dpi=300, bbox_inches='tight')
    plt.show()

# Function to plot top 20 feature importances with customizable bar color
def plot_top_20_feature_importances(feature_names, feature_importances, output_dir, output_file, bar_color='skyblue'):
    plt.figure(figsize=(7.5, 5.5))
    sorted_indices = np.argsort(feature_importances)[::-1][:20]
    plt.bar(np.array(feature_names)[sorted_indices], np.array(feature_importances)[sorted_indices], color=bar_color)
    plt.xlabel('Feature Names', fontsize=10)  # Reduced font size for x-axis label
    plt.ylabel('Feature Importances', fontsize=10)  # Reduced font size for y-axis label
    plt.title('Feature Importances for Top 20 Features', fontsize=12)  # Reduced font size for title
    plt.xticks(rotation=45, fontsize=8) 
    plt.yticks(fontsize=8)  # Reduced font size for y-axis tick marks and text
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, output_file), dpi=300, bbox_inches='tight')
    plt.show()

# EXPLANA functions 
# placeholder dataframe for creating figure
# def no_features_selected():
#     _df = pd.DataFrame({"important_features": ["no_selected_features"],
#                         "decoded_features": ["NA"],
#                         "feature_importance_vals": [-100]})
#     _df.to_csv(ds_out_path + dataset + "-boruta-important.txt", index=False,
#                sep="\t")


# def fail_analysis(out_file):
#     with open(out_file, "w") as file:
#         file.write("Failed Analysis")
#     no_features_selected()
#     return

# def percent_var_statement(step, oob):
#     if step == "full_forest":
#         # if first percent var explained is less than 5% then terminate
#         if float(oob) < 5.0:
#             fail_analysis(out_file)
#             return
        
def create_parameter_grids(df):
    # Find the row with the lowest mean_mse_score
    lowest_mse_row = df.loc[df['mean_mse_score'].idxmin()]
    print("First 5 columns for the lowest mean_mse_score:")
    print(lowest_mse_row.iloc[:5])

    # Find the row with the lowest mean_prev_score
    lowest_prev_row = df.loc[df['mean_prev'].idxmin()]
    print("First 5 columns for the lowest mean_prev_score:")
    print(lowest_prev_row.iloc[:5])

    # Find the row with the lowest mean_ptev_score
    lowest_ptev_row = df.loc[df['mean_ptev'].idxmin()]
    print("First 5 columns for the lowest mean_ptev_score:")
    print(lowest_ptev_row.iloc[:5])

    # Find the row with the highest oob_score
    highest_oob_row = df.loc[df['oob_score'].idxmax()]
    print("\nFirst 5 columns for the highest oob_score:")
    print(highest_oob_row.iloc[:5])

    # Create parameter grids from the extracted rows
    best_mse_param_grid = {
        'n_estimators': [int(lowest_mse_row['n_estimators'])],
        'max_depth': [None if pd.isna(lowest_mse_row['max_depth']) else int(lowest_mse_row['max_depth'])],
        'min_samples_split': [float(lowest_mse_row['min_samples_split'])],
        'max_iter': [int(lowest_mse_row['max_iter'])],
        'n_splits': [int(lowest_mse_row['n_splits'])]
    }
    print("Best MSE Parameter Grid:")
    print("n_estimators:", best_mse_param_grid['n_estimators'][0])
    print("max_depth:", best_mse_param_grid['max_depth'][0])
    print("min_samples_split:", best_mse_param_grid['min_samples_split'][0])
    print("max_iter:", best_mse_param_grid['max_iter'][0])
    print("n_splits:", best_mse_param_grid['n_splits'][0])

    lowest_prev_param_grid = {
        'n_estimators': [int(lowest_prev_row['n_estimators'])],
        'max_depth': [None if pd.isna(lowest_prev_row['max_depth']) else int(lowest_prev_row['max_depth'])],
        'min_samples_split': [float(lowest_prev_row['min_samples_split'])],
        'max_iter': [int(lowest_prev_row['max_iter'])],
        'n_splits': [int(lowest_prev_row['n_splits'])]
    }
    print("\nLowest Prev Parameter Grid:")
    print("n_estimators:", lowest_prev_param_grid['n_estimators'][0])
    print("max_depth:", lowest_prev_param_grid['max_depth'][0])
    print("min_samples_split:", lowest_prev_param_grid['min_samples_split'][0])
    print("max_iter:", lowest_prev_param_grid['max_iter'][0])
    print("n_splits:", lowest_prev_param_grid['n_splits'][0])

    lowest_ptev_param_grid = {
        'n_estimators': [int(lowest_ptev_row['n_estimators'])],
        'max_depth': [None if pd.isna(lowest_ptev_row['max_depth']) else int(lowest_ptev_row['max_depth'])],
        'min_samples_split': [float(lowest_ptev_row['min_samples_split'])],
        'max_iter': [int(lowest_ptev_row['max_iter'])],
        'n_splits': [int(lowest_ptev_row['n_splits'])]
    }
    print("\nLowest PTEV Parameter Grid:")
    print("n_estimators:", lowest_ptev_param_grid['n_estimators'][0])
    print("max_depth:", lowest_ptev_param_grid['max_depth'][0])
    print("min_samples_split:", lowest_ptev_param_grid['min_samples_split'][0])
    print("max_iter:", lowest_ptev_param_grid['max_iter'][0])
    print("n_splits:", lowest_ptev_param_grid['n_splits'][0])

    highest_oob_param_grid = {
        'n_estimators': [int(highest_oob_row['n_estimators'])],
        'max_depth': [None if pd.isna(highest_oob_row['max_depth']) else int(highest_oob_row['max_depth'])],
        'min_samples_split': [float(highest_oob_row['min_samples_split'])],
        'max_iter': [int(highest_oob_row['max_iter'])],
        'n_splits': [int(highest_oob_row['n_splits'])]
    }
    print("\nHighest OOB Parameter Grid:")
    print("n_estimators:", highest_oob_row['n_estimators'])
    print("max_depth:", highest_oob_row['max_depth'])
    print("min_samples_split:", highest_oob_row['min_samples_split'])
    print("max_iter:", highest_oob_row['max_iter'])
    print("n_splits:", highest_oob_row['n_splits'])


def run_merf_analysis(X, Y, Z, clusters_train, 
                      X_new, Y_new, Z_new, clusters_new,
                      df, 
                      output_dir, r2_out, feature_imp_out):
    import pandas as pd
    import numpy as np
    from merf import MERF
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.metrics import r2_score
    import matplotlib.pyplot as plt
    
    # Find the row with the lowest mean_mse_score
    lowest_mse_row = df.loc[df['mean_mse_score'].idxmin()]
    print("First 5 columns for the lowest mean_mse_score:")
    print(lowest_mse_row.iloc[:5])

    # Find the row with the lowest mean_prev_score
    lowest_prev_row = df.loc[df['mean_prev'].idxmin()]
    print("First 5 columns for the lowest mean_prev_score:")
    print(lowest_prev_row.iloc[:5])

    # Find the row with the lowest mean_ptev_score
    lowest_ptev_row = df.loc[df['mean_ptev'].idxmin()]
    print("First 5 columns for the lowest mean_ptev_score:")
    print(lowest_ptev_row.iloc[:5])

    # Find the row with the highest oob_score
    highest_oob_row = df.loc[df['oob_score'].idxmax()]
    print("\nFirst 5 columns for the highest oob_score:")
    print(highest_oob_row.iloc[:5])

    # Create parameter grids from the extracted rows
    best_mse_param_grid = {
        'n_estimators': [int(lowest_mse_row['n_estimators'])],
        'max_depth': [None if pd.isna(lowest_mse_row['max_depth']) else int(lowest_mse_row['max_depth'])],
        'min_samples_split': [float(lowest_mse_row['min_samples_split'])],
        'max_iter': [int(lowest_mse_row['max_iter'])],
        'n_splits': [int(lowest_mse_row['n_splits'])]
    }

    lowest_prev_param_grid = {
        'n_estimators': [int(lowest_prev_row['n_estimators'])],
        'max_depth': [None if pd.isna(lowest_prev_row['max_depth']) else int(lowest_prev_row['max_depth'])],
        'min_samples_split': [float(lowest_prev_row['min_samples_split'])],
        'max_iter': [int(lowest_prev_row['max_iter'])],
        'n_splits': [int(lowest_prev_row['n_splits'])]
    }

    lowest_ptev_param_grid = {
        'n_estimators': [int(lowest_ptev_row['n_estimators'])],
        'max_depth': [None if pd.isna(lowest_ptev_row['max_depth']) else int(lowest_ptev_row['max_depth'])],
        'min_samples_split': [float(lowest_ptev_row['min_samples_split'])],
        'max_iter': [int(lowest_ptev_row['max_iter'])],
        'n_splits': [int(lowest_ptev_row['n_splits'])]
    }

    highest_oob_param_grid = {
        'n_estimators': [int(highest_oob_row['n_estimators'])],
        'max_depth': [None if pd.isna(highest_oob_row['max_depth']) else int(highest_oob_row['max_depth'])],
        'min_samples_split': [float(highest_oob_row['min_samples_split'])],
        'max_iter': [int(highest_oob_row['max_iter'])],
        'n_splits': [int(highest_oob_row['n_splits'])]
    }

    # Create MERF models for each parameter grid
    mse_merf = MERF(fixed_effects_model=RandomForestRegressor(
        n_estimators=best_mse_param_grid['n_estimators'][0],
        max_depth=best_mse_param_grid['max_depth'][0],
        min_samples_split=best_mse_param_grid['min_samples_split'][0],
        n_jobs=1,
        oob_score=True),
        gll_early_stop_threshold=None,
        max_iterations=best_mse_param_grid['max_iter'][0])

    prev_merf = MERF(fixed_effects_model=RandomForestRegressor(
        n_estimators=lowest_prev_param_grid['n_estimators'][0],
        max_depth=lowest_prev_param_grid['max_depth'][0],
        min_samples_split=lowest_prev_param_grid['min_samples_split'][0],
        n_jobs=1,
        oob_score=True),
        gll_early_stop_threshold=None,
        max_iterations=lowest_prev_param_grid['max_iter'][0])

    ptev_merf = MERF(fixed_effects_model=RandomForestRegressor(
        n_estimators=lowest_ptev_param_grid['n_estimators'][0],
        max_depth=lowest_ptev_param_grid['max_depth'][0],
        min_samples_split=lowest_ptev_param_grid['min_samples_split'][0],
        n_jobs=1,
        oob_score=True),
        gll_early_stop_threshold=None,
        max_iterations=lowest_ptev_param_grid['max_iter'][0])

    oob_merf = MERF(fixed_effects_model=RandomForestRegressor(
        n_estimators=highest_oob_param_grid['n_estimators'][0],
        max_depth=highest_oob_param_grid['max_depth'][0],
        min_samples_split=highest_oob_param_grid['min_samples_split'][0],
        n_jobs=1,
        oob_score=True),
        gll_early_stop_threshold=None,
        max_iterations=highest_oob_param_grid['max_iter'][0])

    print("---------- RUN MERF RAW WITH TUNING PARAMETERS 🌱 ----------")
    mrf_mse = mse_merf.fit(X.select_dtypes(include=[np.number]), Z, pd.Series(clusters_train), Y)
    mrf_prev = prev_merf.fit(X.select_dtypes(include=[np.number]), Z, pd.Series(clusters_train), Y)
    mrf_ptev = ptev_merf.fit(X.select_dtypes(include=[np.number]), Z, pd.Series(clusters_train), Y)
    mrf_oob = oob_merf.fit(X.select_dtypes(include=[np.number]), Z, pd.Series(clusters_train), Y)

    # Predict using the fitted model
    y_hat_new_mse = mrf_mse.predict(X_new, Z_new, clusters_new)
    y_hat_new_prev = mrf_prev.predict(X_new, Z_new, clusters_new)
    y_hat_new_ptev = mrf_ptev.predict(X_new, Z_new, clusters_new)
    y_hat_new_oob = mrf_oob.predict(X_new, Z_new, clusters_new)

    ### Calculate R-squared values for each model
    r2_mse = r2_score(Y_new, y_hat_new_mse)
    print(f"R-squared for MSE Model: {r2_mse:.4f}")
    r2_prev = r2_score(Y_new, y_hat_new_prev)
    print(f"R-squared for Prev Model: {r2_prev:.4f}") 
    r2_ptev = r2_score(Y_new, y_hat_new_ptev)
    print(f"R-squared for PTEV Model: {r2_ptev:.4f}")
    r2_oob = r2_score(Y_new, y_hat_new_oob)
    print(f"R-squared for OOB Model: {r2_oob:.4f}")

    # Store R-squared values in a dictionary for plotting
    r2_values = {
        'MSE Model': r2_mse,
        'Prev Model': r2_prev,
        'PTEV Model': r2_ptev,
        'OOB Model': r2_oob
    }

    # Plot R-squared values
    plt.figure(figsize=(8, 5))
    plt.bar(r2_values.keys(), r2_values.values(), color=['#F88F79', '#F0F879', '#ACF0F8', '#86B874'])
    plt.ylabel('R-squared Value')
    plt.title('R-squared Comparison of MERF Models')
    plt.ylim(0, 0.5)
    plt.xticks(rotation=45)
    plt.savefig(os.path.join(output_dir, r2_out), dpi=300, bbox_inches='tight')
    plt.show()

    ## Plot feature importances 
    mse_forest = mrf_mse.trained_fe_model
    mse_feature_names = mse_forest.feature_names_in_
    mse_feature_importances = mse_forest.feature_importances_

    prev_forest = mrf_prev.trained_fe_model
    prev_feature_importances = prev_forest.feature_importances_

    ptev_forest = mrf_ptev.trained_fe_model
    ptev_feature_importances = ptev_forest.feature_importances_

    oob_forest = mrf_oob.trained_fe_model
    oob_feature_importances = oob_forest.feature_importances_

    importances = {
        'MSE': mse_feature_importances,
        'Prev': prev_feature_importances,
        'PTEV': ptev_feature_importances,
        'OOB': oob_feature_importances
    }

    # Create a DataFrame for easier plotting
    importances_df = pd.DataFrame(importances, index=mse_feature_names)

    # Sort the DataFrame by importance in descending order
    importances_df = importances_df.sort_values(by=importances_df.columns.tolist(), ascending=False)

    # Plot top feature importances
    top_n = 20  # Number of top features to display
    importances_df.head(top_n).plot(kind='bar', figsize=(12, 8), color=['#F88F79', '#F0F879', '#ACF0F8', '#86B874'])
    plt.title('Top Feature Importances')
    plt.ylim(0, 0.5)
    plt.ylabel('Importance')
    plt.xlabel('Features')
    plt.xticks(rotation=45, ha='right')  # Adjusted rotation and horizontal alignment to reduce overlap
    plt.legend(title='Models', loc='upper right', fontsize='small')  # Adjusted legend position to top right corner and made it smaller
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, feature_imp_out), dpi=300, bbox_inches='tight')
    plt.show()

# Example of how to call the function
# run_merf_analysis(X, Y, Z, clusters_train, X_new, Y_new, Z_new, clusters_new, output_dir, r2_out, feature_imp_out)

def run_merf_analysis_old2(X, Y, Z, clusters_train, 
                      X_new, Y_new, Z_new, clusters_new,
                      df, 
                      output_dir, r2_out, feature_imp_out):
    import pandas as pd
    import numpy as np
    from merf import MERF
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.metrics import r2_score
    import matplotlib.pyplot as plt
    
    # Find the row with the lowest mean_mse_score
    lowest_mse_row = df.loc[df['mean_mse_score'].idxmin()]
    print("First 5 columns for the lowest mean_mse_score:")
    print(lowest_mse_row.iloc[:5])

    # Find the row with the lowest mean_prev_score
    lowest_prev_row = df.loc[df['mean_prev'].idxmin()]
    print("First 5 columns for the lowest mean_prev_score:")
    print(lowest_prev_row.iloc[:5])

    # Find the row with the lowest mean_ptev_score
    lowest_ptev_row = df.loc[df['mean_ptev'].idxmin()]
    print("First 5 columns for the lowest mean_ptev_score:")
    print(lowest_ptev_row.iloc[:5])

    # Find the row with the highest oob_score
    highest_oob_row = df.loc[df['oob_score'].idxmax()]
    print("\nFirst 5 columns for the highest oob_score:")
    print(highest_oob_row.iloc[:5])

    # Create parameter grids from the extracted rows
    best_mse_param_grid = {
        'n_estimators': [int(lowest_mse_row['n_estimators'])],
        'max_depth': [None if pd.isna(lowest_mse_row['max_depth']) else int(lowest_mse_row['max_depth'])],
        'min_samples_split': [float(lowest_mse_row['min_samples_split'])],
        'max_iter': [int(lowest_mse_row['max_iter'])],
        'n_splits': [int(lowest_mse_row['n_splits'])]
    }

    lowest_prev_param_grid = {
        'n_estimators': [int(lowest_prev_row['n_estimators'])],
        'max_depth': [None if pd.isna(lowest_prev_row['max_depth']) else int(lowest_prev_row['max_depth'])],
        'min_samples_split': [float(lowest_prev_row['min_samples_split'])],
        'max_iter': [int(lowest_prev_row['max_iter'])],
        'n_splits': [int(lowest_prev_row['n_splits'])]
    }

    lowest_ptev_param_grid = {
        'n_estimators': [int(lowest_ptev_row['n_estimators'])],
        'max_depth': [None if pd.isna(lowest_ptev_row['max_depth']) else int(lowest_ptev_row['max_depth'])],
        'min_samples_split': [float(lowest_ptev_row['min_samples_split'])],
        'max_iter': [int(lowest_ptev_row['max_iter'])],
        'n_splits': [int(lowest_ptev_row['n_splits'])]
    }

    highest_oob_param_grid = {
        'n_estimators': [int(highest_oob_row['n_estimators'])],
        'max_depth': [None if pd.isna(highest_oob_row['max_depth']) else int(highest_oob_row['max_depth'])],
        'min_samples_split': [float(highest_oob_row['min_samples_split'])],
        'max_iter': [int(highest_oob_row['max_iter'])],
        'n_splits': [int(highest_oob_row['n_splits'])]
    }

    # Create MERF models for each parameter grid
    mse_merf = MERF(fixed_effects_model=RandomForestRegressor(
        n_estimators=best_mse_param_grid['n_estimators'][0],
        max_depth=best_mse_param_grid['max_depth'][0],
        min_samples_split=best_mse_param_grid['min_samples_split'][0],
        n_jobs=1,
        oob_score=True),
        gll_early_stop_threshold=None,
        max_iterations=best_mse_param_grid['max_iter'][0])

    prev_merf = MERF(fixed_effects_model=RandomForestRegressor(
        n_estimators=lowest_prev_param_grid['n_estimators'][0],
        max_depth=lowest_prev_param_grid['max_depth'][0],
        min_samples_split=lowest_prev_param_grid['min_samples_split'][0],
        n_jobs=1,
        oob_score=True),
        gll_early_stop_threshold=None,
        max_iterations=lowest_prev_param_grid['max_iter'][0])

    ptev_merf = MERF(fixed_effects_model=RandomForestRegressor(
        n_estimators=lowest_ptev_param_grid['n_estimators'][0],
        max_depth=lowest_ptev_param_grid['max_depth'][0],
        min_samples_split=lowest_ptev_param_grid['min_samples_split'][0],
        n_jobs=1,
        oob_score=True),
        gll_early_stop_threshold=None,
        max_iterations=lowest_ptev_param_grid['max_iter'][0])

    oob_merf = MERF(fixed_effects_model=RandomForestRegressor(
        n_estimators=highest_oob_param_grid['n_estimators'][0],
        max_depth=highest_oob_param_grid['max_depth'][0],
        min_samples_split=highest_oob_param_grid['min_samples_split'][0],
        n_jobs=1,
        oob_score=True),
        gll_early_stop_threshold=None,
        max_iterations=highest_oob_param_grid['max_iter'][0])

    print("---------- RUN MERF RAW WITH TUNING PARAMETERS 🌱 ----------")
    mrf_mse = mse_merf.fit(X.select_dtypes(include=[np.number]), Z, pd.Series(clusters_train), Y)
    mrf_prev = prev_merf.fit(X.select_dtypes(include=[np.number]), Z, pd.Series(clusters_train), Y)
    mrf_ptev = ptev_merf.fit(X.select_dtypes(include=[np.number]), Z, pd.Series(clusters_train), Y)
    mrf_oob = oob_merf.fit(X.select_dtypes(include=[np.number]), Z, pd.Series(clusters_train), Y)

    # Predict using the fitted model
    y_hat_new_mse = mrf_mse.predict(X_new, Z_new, clusters_new)
    y_hat_new_prev = mrf_prev.predict(X_new, Z_new, clusters_new)
    y_hat_new_ptev = mrf_ptev.predict(X_new, Z_new, clusters_new)
    y_hat_new_oob = mrf_oob.predict(X_new, Z_new, clusters_new)

    # Calculate R-squared values for each model
    r2_mse = r2_score(Y_new, y_hat_new_mse)
    r2_prev = r2_score(Y_new, y_hat_new_prev)
    r2_ptev = r2_score(Y_new, y_hat_new_ptev)
    r2_oob = r2_score(Y_new, y_hat_new_oob)

    # Store R-squared values in a dictionary for plotting
    r2_values = {
        'MSE Model': r2_mse,
        'Prev Model': r2_prev,
        'PTEV Model': r2_ptev,
        'OOB Model': r2_oob
    }
    # Plot R-squared values
    plt.figure(figsize=(8, 5))
    plt.bar(r2_values.keys(), r2_values.values(), color=['#F88F79', '#F0F879', '#ACF0F8', '#86B874'])
    plt.ylabel('R-squared Value')
    plt.title('R-squared Comparison of MERF Models')
    plt.ylim(0, 1)  # Adjusted to show full range of R-squared values
    plt.xticks(rotation=45)
    plt.savefig(os.path.join(output_dir, r2_out), dpi=300, bbox_inches='tight')
    plt.show()

    # Determine the model with the highest R-squared value
    best_model_name = max(r2_values, key=r2_values.get)
    print(f"Best model: {best_model_name} with R-squared: {r2_values[best_model_name]:.4f}")

    # Get the feature importances for the best model
    if best_model_name == 'MSE Model':
        feature_importances = mrf_mse.trained_fe_model.feature_importances_
    elif best_model_name == 'Prev Model':
        feature_importances = mrf_prev.trained_fe_model.feature_importances_
    elif best_model_name == 'PTEV Model':
        feature_importances = mrf_ptev.trained_fe_model.feature_importances_
    else:  # 'OOB Model'
        feature_importances = mrf_oob.trained_fe_model.feature_importances_

    # Create a DataFrame for the top 15 feature importances
    feature_names = mrf_mse.trained_fe_model.feature_names_in_  # Assuming all models have the same feature names
    feature_importance_df = pd.DataFrame({
        'Feature': feature_names,
        'Importance': feature_importances
    })

    # Sort the DataFrame by importance and get the top 15 features
    top_features_df = feature_importance_df.sort_values(by='Importance', ascending=False).head(15)
    print("Top 15 feature importances for the best model:")
    print(top_features_df)

    # Return the R-squared values and the top features DataFrame for further analysis
    return r2_values, top_features_df

def run_merf_analysis2(X, Y, Z, clusters_train, 
                      X_new, Y_new, Z_new, clusters_new,
                      df, 
                      output_dir, r2_out, feature_imp_out, results_filename, time_new):
    import pandas as pd
    import numpy as np
    from merf import MERF
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.metrics import r2_score
    import matplotlib.pyplot as plt
    
    # Find the row with the lowest mean_mse_score
    lowest_mse_row = df.loc[df['mean_mse_score'].idxmin()]
    print("First 5 columns for the lowest mean_mse_score:")
    print(lowest_mse_row.iloc[:5])

    # Find the row with the lowest mean_prev_score
    lowest_prev_row = df.loc[df['mean_prev'].idxmin()]
    print("First 5 columns for the lowest mean_prev_score:")
    print(lowest_prev_row.iloc[:5])

    # Find the row with the lowest mean_ptev_score
    lowest_ptev_row = df.loc[df['mean_ptev'].idxmin()]
    print("First 5 columns for the lowest mean_ptev_score:")
    print(lowest_ptev_row.iloc[:5])

    # Find the row with the highest oob_score
    highest_oob_row = df.loc[df['oob_score'].idxmax()]
    print("\nFirst 5 columns for the highest oob_score:")
    print(highest_oob_row.iloc[:5])

    # Create parameter grids from the extracted rows
    best_mse_param_grid = {
        'n_estimators': [int(lowest_mse_row['n_estimators'])],
        'max_depth': [None if pd.isna(lowest_mse_row['max_depth']) else int(lowest_mse_row['max_depth'])],
        'min_samples_split': [float(lowest_mse_row['min_samples_split'])],
        'max_iter': [int(lowest_mse_row['max_iter'])],
        'n_splits': [int(lowest_mse_row['n_splits'])]
    }

    lowest_prev_param_grid = {
        'n_estimators': [int(lowest_prev_row['n_estimators'])],
        'max_depth': [None if pd.isna(lowest_prev_row['max_depth']) else int(lowest_prev_row['max_depth'])],
        'min_samples_split': [float(lowest_prev_row['min_samples_split'])],
        'max_iter': [int(lowest_prev_row['max_iter'])],
        'n_splits': [int(lowest_prev_row['n_splits'])]
    }

    lowest_ptev_param_grid = {
        'n_estimators': [int(lowest_ptev_row['n_estimators'])],
        'max_depth': [None if pd.isna(lowest_ptev_row['max_depth']) else int(lowest_ptev_row['max_depth'])],
        'min_samples_split': [float(lowest_ptev_row['min_samples_split'])],
        'max_iter': [int(lowest_ptev_row['max_iter'])],
        'n_splits': [int(lowest_ptev_row['n_splits'])]
    }

    highest_oob_param_grid = {
        'n_estimators': [int(highest_oob_row['n_estimators'])],
        'max_depth': [None if pd.isna(highest_oob_row['max_depth']) else int(highest_oob_row['max_depth'])],
        'min_samples_split': [float(highest_oob_row['min_samples_split'])],
        'max_iter': [int(highest_oob_row['max_iter'])],
        'n_splits': [int(highest_oob_row['n_splits'])]
    }

    # Create a DataFrame to store results for each model
    results_df = pd.DataFrame(columns=['Model', 'y_hat_new', 
                                       'R_squared', 'Top_15_Feature_Importances', 'Cluster'])  # Added 'Cluster' column

    # Create MERF models for each parameter grid
    mse_merf = MERF(fixed_effects_model=RandomForestRegressor(
        n_estimators=best_mse_param_grid['n_estimators'][0],
        max_depth=best_mse_param_grid['max_depth'][0],
        min_samples_split=best_mse_param_grid['min_samples_split'][0],
        n_jobs=1,
        oob_score=True),
        gll_early_stop_threshold=None,
        max_iterations=best_mse_param_grid['max_iter'][0])

    prev_merf = MERF(fixed_effects_model=RandomForestRegressor(
        n_estimators=lowest_prev_param_grid['n_estimators'][0],
        max_depth=lowest_prev_param_grid['max_depth'][0],
        min_samples_split=lowest_prev_param_grid['min_samples_split'][0],
        n_jobs=1,
        oob_score=True),
        gll_early_stop_threshold=None,
        max_iterations=lowest_prev_param_grid['max_iter'][0])

    ptev_merf = MERF(fixed_effects_model=RandomForestRegressor(
        n_estimators=lowest_ptev_param_grid['n_estimators'][0],
        max_depth=lowest_ptev_param_grid['max_depth'][0],
        min_samples_split=lowest_ptev_param_grid['min_samples_split'][0],
        n_jobs=1,
        oob_score=True),
        gll_early_stop_threshold=None,
        max_iterations=lowest_ptev_param_grid['max_iter'][0])

    oob_merf = MERF(fixed_effects_model=RandomForestRegressor(
        n_estimators=highest_oob_param_grid['n_estimators'][0],
        max_depth=highest_oob_param_grid['max_depth'][0],
        min_samples_split=highest_oob_param_grid['min_samples_split'][0],
        n_jobs=1,
        oob_score=True),
        gll_early_stop_threshold=None,
        max_iterations=highest_oob_param_grid['max_iter'][0])

    # Fit models and predict
    for model_name, merf_model in zip(['MSE Model', 'Prev Model', 'PTEV Model', 'OOB Model'], 
                                       [mse_merf, prev_merf, ptev_merf, oob_merf]):
        print(f"---------- RUN {model_name} WITH TUNING PARAMETERS 🌱 ----------")
        mrf_fit = merf_model.fit(X.select_dtypes(include=[np.number]), Z, pd.Series(clusters_train), Y)
        y_hat_new = mrf_fit.predict(X_new, Z_new, clusters_new)
        r2 = r2_score(Y_new, y_hat_new)

        # Get feature importances
        feature_importances = mrf_fit.trained_fe_model.feature_importances_
        feature_names = mrf_fit.trained_fe_model.feature_names_in_
        feature_importance_df = pd.DataFrame({
            'Feature': feature_names,
            'Importance': feature_importances
        })
        top_features = feature_importance_df.sort_values(by='Importance', ascending=False).head(15)

        # Append results to the DataFrame for each y_hat_new value
        for y_hat, cluster, time in zip(y_hat_new, clusters_new, time_new):  # Associate each y_hat_new with its corresponding cluster and time
            new_row = pd.DataFrame({
                'Model': [model_name],
                'y_hat_new': [y_hat],  # Single y_hat_new value
                'R_squared': [r2],
                'Top_15_Feature_Importances': [top_features.to_dict(orient='records')],  # Convert to dict for CSV
                'Cluster': [cluster],  # Added cluster information
                'Time': [time]  # Added time information
            })
            results_df = pd.concat([results_df, new_row], ignore_index=True)

    # Save results to CSV
    results_df.to_csv(os.path.join(output_dir, results_filename), index=False)

    # Calculate R-squared values for each model
    r2_values = {
        'MSE Model': r2_score(Y_new, mse_merf.predict(X_new, Z_new, clusters_new)),
        'Prev Model': r2_score(Y_new, prev_merf.predict(X_new, Z_new, clusters_new)),
        'PTEV Model': r2_score(Y_new, ptev_merf.predict(X_new, Z_new, clusters_new)),
        'OOB Model': r2_score(Y_new, oob_merf.predict(X_new, Z_new, clusters_new))
    }

    # Plot R-squared values
    plt.figure(figsize=(8, 5))
    plt.bar(r2_values.keys(), r2_values.values(), color=['#F88F79', '#F0F879', '#ACF0F8', '#86B874'])
    plt.ylabel('R-squared Value')
    plt.title('R-squared Comparison of MERF Models')
    plt.ylim(0, 1)  # Adjusted to show full range of R-squared values
    plt.xticks(rotation=45)
    plt.savefig(os.path.join(output_dir, r2_out), dpi=300, bbox_inches='tight')
    plt.show()

    # Determine the model with the highest R-squared value
    best_model_name = max(r2_values, key=r2_values.get)
    print(f"Best model: {best_model_name} with R-squared: {r2_values[best_model_name]:.4f}")

    # Return the R-squared values and the top features DataFrame for further analysis
    return r2_values, results_df

def compare_r2_values1(model_names, *r2_dicts):
    # Create a DataFrame to hold the R-squared values
    r2_comparison_df = pd.DataFrame(r2_dicts).T
    r2_comparison_df.columns = model_names 
    print(r2_comparison_df.head(10).columns[:10])
    print(r2_comparison_df.head(10))

    # Sort the DataFrame by R-squared values in descending order
    r2_comparison_df = r2_comparison_df.sort_values(by=model_names, ascending=False)

    # Set up the bar positions
    num_runs = len(r2_comparison_df.columns)
    bar_width = 0.15  # Width of each bar
    x = np.arange(len(r2_comparison_df))  # The label locations

    # Create a color map for the models
    colors = plt.cm.viridis(np.linspace(0, 1, num_runs))

    # Create bars for each model in each run
    for i in range(num_runs):
        plt.bar(x + i * bar_width, r2_comparison_df.iloc[:, i],  # Use plt.bar for vertical bars
                 width=bar_width, label=model_names[i], color=colors[i])  # Use model names for labels

    # Add labels, title, and custom x-axis tick labels
    plt.ylabel('R-squared Value')  # Change to y-label
    plt.title('R-squared Comparison of Different Runs')
    plt.ylim(0, 0.2)  # Assuming R-squared values are between 0 and 1
    plt.xticks(x + bar_width * (num_runs - 1) / 2, r2_comparison_df.index)  # Center the x-ticks
    plt.legend(title='Models', loc='upper right', fontsize=10, 
               title_fontsize=10, ncol = 2, frameon=False)
    plt.tight_layout()
    plt.show()



def compare_r2_values2(model_names, *r2_dicts):
    # Create a DataFrame to hold the R-squared values
    r2_comparison_df = pd.DataFrame(r2_dicts).T
    r2_comparison_df.columns = model_names 

    # Reset index to use row names as a column
    r2_comparison_df.reset_index(inplace=True)
    r2_comparison_df = r2_comparison_df.melt(id_vars='index', var_name='Model', value_name='R-squared')

    # Rename the 'index' column to 'Metric' for clarity
    r2_comparison_df.rename(columns={'index': 'Metric'}, inplace=True)

    # Create the bar plot
    plt.figure(figsize=(12, 8))
    sns.barplot(data=r2_comparison_df, x='Model', y='R-squared', hue='Metric', palette='viridis')

    # Add labels and title
    plt.ylabel('R-squared Value')
    plt.title('R-squared Comparison of Different Models')
    plt.ylim(-0.2, 0.2)  # Adjust y-limits based on your data range
    plt.xticks(rotation=45, ha='right')
    plt.legend(title='hyperparameter optimization settings', loc='upper right', fontsize=10, 
               title_fontsize=10, ncol = 2, frameon=False)
    plt.tight_layout()
    plt.show()



def compare_r2_values3(model_names, *r2_dicts):
    # Create a DataFrame to hold the R-squared values
    r2_comparison_df = pd.DataFrame(r2_dicts).T
    r2_comparison_df.columns = model_names 

    # Reset index to use row names as a column
    r2_comparison_df.reset_index(inplace=True)
    r2_comparison_df = r2_comparison_df.melt(id_vars='index', var_name='Model', value_name='R-squared')

    # Rename the 'index' column to 'Metric' for clarity
    r2_comparison_df.rename(columns={'index': 'Metric'}, inplace=True)

    # Sort the DataFrame by 'Model' and then by 'R-squared' in descending order
    r2_comparison_df = r2_comparison_df.sort_values(by=['Metric', 'R-squared'], ascending=[True, False])

    # Create the bar plot
    plt.figure(figsize=(12, 8))
    sns.barplot(data=r2_comparison_df, x='Model', y='R-squared', hue='Metric', palette='viridis')

    # Add labels and title
    plt.ylabel('R-squared Value')
    plt.title('R-squared Comparison of Different Models')
    plt.ylim(-0.2, 0.2)  # Adjust y-limits based on your data range
    plt.xticks(rotation=45, ha='right')
    plt.legend(title='hyperparameter optimization settings', loc='upper right', fontsize=10, 
               title_fontsize=10, ncol = 2, frameon=False)
    plt.tight_layout()
    plt.show()


import pandas as pd
import matplotlib.pyplot as plt

def plot_top_feature_importances_comparative(feature_importance_dfs):
    # Create a DataFrame to hold all feature importances
    combined_importances = pd.DataFrame()

    # Iterate through each feature importance DataFrame
    for i, df in enumerate(feature_importance_dfs):
        # Add a column to identify the model
        df['Model'] = f'Model {i + 1}'  # Adjust the naming as needed
        combined_importances = pd.concat([combined_importances, df], ignore_index=True)

    # Group by feature and aggregate the importances (mean or sum)
    aggregated_importances = combined_importances.groupby('Feature')['Importance'].mean().reset_index()

    # Sort the DataFrame by importance and get the top 15 features
    top_features_df = aggregated_importances.sort_values(by='Importance', ascending=False).head(15)

    # Plot the top feature importances
    plt.figure(figsize=(12, 8))
    plt.barh(top_features_df['Feature'], top_features_df['Importance'], color='skyblue')
    plt.xlabel('Importance')
    plt.title('Top 15 Feature Importances Across Models')
    plt.gca().invert_yaxis()  # Invert y-axis to have the highest importance on top
    plt.tight_layout()
    plt.show()

# Example usage:
# Assuming you have a list of DataFrames from previous runs
# feature_importance_dfs = [df1, df2, df3, ...]  # Replace with your actual DataFrames
# plot_top_feature_importances_comparative(feature_importance_dfs)

def run_merf(X, Y, Z, clusters_train, 
              X_new, Y_new, Z_new, clusters_new,
              param_grid_list, results_name):
    from merf import MERF
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.metrics import r2_score
    import matplotlib.pyplot as plt

    results_name = []  # List to store results for each model

    for params in param_grid_list:
        n_estimators = params['n_estimators']
        max_depth = params['max_depth']
        min_samples_split = params['min_samples_split']
        max_iter = params['max_iter']

        merf = MERF(fixed_effects_model=RandomForestRegressor(
            n_estimators=n_estimators,
            max_depth=max_depth,
            min_samples_split=min_samples_split,
            n_jobs=1,
            oob_score=True),
            gll_early_stop_threshold=None,
            max_iterations=max_iter)

        print("---------- RUN MERF RAW WITH TUNING PARAMETERS 🌱 ----------")
        mrf_fit = merf.fit(X.select_dtypes(include=[np.number]), Z, pd.Series(clusters_train), Y)
        y_hat_new = mrf_fit.predict(X_new, Z_new, clusters_new)  # Predict using the fitted model
        r2 = r2_score(Y_new, y_hat_new)
        mrf_fit_fe = mrf_fit.trained_fe_model
        print(f"R-squared for MERF Model: {r2:.4f}")

        # Store results in a dictionary and print as it runs
        model_results = {
            'params': params,
            'predictions': y_hat_new,
            'r2': r2,
            'feature_importances': mrf_fit_fe.feature_importances_
        }
        # print(f"Model Results: {model_results}")  # Print the results for this model as it runs
        results_name.append(model_results)  # Append the results for this model

    return results_name  # Return the list of results for all models

# Function to plot R-squared values side by side with color coding
def plot_r2_comparison(results_list):
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.cm import ScalarMappable

    # Debugging: Print the contents of results_list
    print("Contents of results_list:", results_list)

    r2_values = []
    model_names = []
    colors = []

    for i, result in enumerate(results_list):
        if isinstance(result, dict) and 'r2' in result:  # Check if result is a dict and has 'r2'
            r2_values.append(result['r2'])
            # Use a specific model name if available, otherwise fallback to generic naming
            model_name = result.get('params', {}).get('max_depth', 'N/A')  # Example of using a parameter
            model_names.append(f"Model {i+1} (max_depth={model_name})")
            colors.append(result['r2'])  # Use R-squared value for color coding
        else:
            print(f"Warning: Expected a dictionary with 'r2' in results, got: {result}")

    # Normalize colors for colormap
    norm = plt.Normalize(min(colors), max(colors))
    cmap = plt.get_cmap('RdYlGn')  # Red to Green colormap
    color_map = ScalarMappable(norm=norm, cmap=cmap)

    # Create a bar plot if we have valid r2_values
    if r2_values:
        plt.figure(figsize=(10, 6))
        x_pos = np.arange(len(model_names))
        plt.bar(x_pos, r2_values, color=color_map.to_rgba(colors))
        plt.xticks(x_pos, model_names, rotation=45)
        plt.ylabel('R-squared Value')
        plt.title('R-squared Comparison of Different Models')
        plt.ylim(0, 1)  # Assuming R-squared values are between 0 and 1
        plt.tight_layout()
        plt.show()
    else:
        print("No valid R-squared values to plot.")

# Function to plot and compare top 10 feature importances with color coding
def plot_feature_importances_comparison(results_list):
    import matplotlib.pyplot as plt
    import pandas as pd

    # Prepare a DataFrame to hold feature importances
    feature_importances = {}
    model_labels = {}
    
    for i, result in enumerate(results_list):
        if isinstance(result, dict) and 'feature_importances' in result:  # Check if result is a dict and has 'feature_importances'
            # Get feature importances and corresponding feature names
            feature_importance = result['feature_importances']
            feature_names = [f'Feature {j+1}' for j in range(len(feature_importance))]
            feature_importances[f'Model {i+1}'] = feature_importance
            model_labels[f'Model {i+1}'] = result.get('params', {})  # Store parameters for labeling
        else:
            print(f"Warning: 'feature_importances' key not found in result for Model {i+1}")

    # Create a DataFrame for easier plotting
    importances_df = pd.DataFrame(feature_importances)

    # Sort the DataFrame by importance in descending order and take top 10
    top_n = 10
    importances_df = importances_df.apply(lambda x: x.nlargest(top_n), axis=0)

    # Plot top feature importances
    ax = importances_df.plot(kind='bar', figsize=(12, 8), color=plt.cm.viridis(np.linspace(0, 1, len(importances_df.columns))))
    plt.title('Top 10 Feature Importances Comparison')
    plt.ylabel('Importance')
    plt.xlabel('Features')
    plt.xticks(rotation=45, ha='right')  # Adjusted rotation and horizontal alignment to reduce overlap
    plt.legend(title='Models', loc='upper right', fontsize='small')  # Adjusted legend position

    # Add model parameters as annotations
    for i, model in enumerate(model_labels.keys()):
        ax.text(i, 0, f"{model_labels[model]}", ha='center', va='bottom', fontsize=8)

    plt.tight_layout()
    plt.show()

# Example of how to call the function
# run_merf_analysis(X, Y, Z, clusters_train, X_new, Y_new, Z_new, clusters_new, output_dir, r2_out, feature_imp_out)

# Plot R-squared comparison
#plot_r2_comparison(results_list)

# Plot feature importances comparison
#plot_feature_importances_comparison(results_list)
def plot_r2_values(results_list, model_names):
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.cm import ScalarMappable

    r2_values = []
    colors = []

    # Check if results_list is empty
    if not results_list:
        print("Warning: results_list is empty.")
        return

    for result in results_list:
        if isinstance(result, dict) and 'r2' in result:  # Check if result is a dict and has 'r2'
            r2_values.append(result['r2'])
            colors.append(result['r2'])  # Use R-squared value for color coding
        else:
            print(f"Warning: Expected a dictionary with 'r2' in results, got: {result}")

    # Check if r2_values is empty after processing
    if not r2_values:
        print("Warning: No valid R-squared values found in results_list.")
        return

    # Normalize colors for colormap
    norm = plt.Normalize(min(colors), max(colors))
    cmap = plt.get_cmap('RdYlGn')  # Red to Green colormap
    color_map = ScalarMappable(norm=norm, cmap=cmap)

    # Create a bar plot if we have valid r2_values
    plt.figure(figsize=(10, 6))
    x_pos = np.arange(len(model_names))
    
    # Ensure r2_values and colors match the length of model_names
    if len(r2_values) != len(model_names) or len(colors) != len(model_names):
        print("Warning: Length mismatch between r2_values/colors and model_names.")
        return

    plt.bar(x_pos, r2_values, color=color_map.to_rgba(colors))
    plt.xticks(x_pos, model_names, rotation=45)
    plt.ylabel('R-squared Value')
    plt.title('R-squared Comparison of Different Models')
    plt.ylim(0, 1)  # Assuming R-squared values are between 0 and 1
    plt.tight_layout()
    plt.show()

# Example usage:
# model_names = ['Only Meta', 'Only GRS', 'Meta GRS Tax', 'Meta GRS Tax Func', 'All But Meta', 'All Omic']
# plot_r2_values(results_list, model_names)

# Function to plot top 10 feature importances for each model in results_list
def plot_top_feature_importances_comparative(results_list):
    import matplotlib.pyplot as plt
    import pandas as pd

    # Debugging: Check the contents of results_list
    print("Results List:", results_list)

    # Prepare a DataFrame to hold feature importances
    feature_importances = {}
    
    for i, result in enumerate(results_list):
        if isinstance(result, dict) and 'feature_importances' in result:
            feature_importance = result['feature_importances']
            feature_names = [f'Feature {j+1}' for j in range(len(feature_importance))]
            feature_importances[f'Model {i+1}'] = feature_importance
        else:
            print(f"Warning: 'feature_importances' key not found in result for Model {i+1}")

    # Create a DataFrame for easier plotting
    importances_df = pd.DataFrame(feature_importances)

    # Debugging: Check the DataFrame shape and contents
    print("Importances DataFrame Shape:", importances_df.shape)
    print("Importances DataFrame Contents:\n", importances_df)

    # Sort the DataFrame by importance in descending order and take top 10
    top_n = 10
    importances_df = importances_df.apply(lambda x: x.nlargest(top_n), axis=0)

    # Plot top feature importances
    ax = importances_df.plot(kind='bar', figsize=(12, 8), color=plt.cm.viridis(np.linspace(0, 1, len(importances_df.columns))))
    plt.title('Top 10 Feature Importances Comparison')
    plt.ylabel('Importance')
    plt.xlabel('Features')
    plt.xticks(rotation=45, ha='right')
    plt.legend(title='Models', loc='upper right', fontsize='small')

    plt.tight_layout()
    plt.show()

