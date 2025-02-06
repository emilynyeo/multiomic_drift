import os
import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from merf.merf import MERF
from sklearn.metrics import r2_score 

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