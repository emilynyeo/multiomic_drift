import os
import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
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
def no_features_selected():
    _df = pd.DataFrame({"important_features": ["no_selected_features"],
                        "decoded_features": ["NA"],
                        "feature_importance_vals": [-100]})
    _df.to_csv(ds_out_path + dataset + "-boruta-important.txt", index=False,
               sep="\t")


def fail_analysis(out_file):
    with open(out_file, "w") as file:
        file.write("Failed Analysis")
    no_features_selected()
    return

def percent_var_statement(step, oob):
    if step == "full_forest":
        # if first percent var explained is less than 5% then terminate
        if float(oob) < 5.0:
            fail_analysis(out_file)
            return