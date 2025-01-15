import torch
from torch import nn, optim
from torch.utils.data import DataLoader, TensorDataset, random_split
import numpy as np
import pandas as pd
import copy
import random
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from scipy.stats import pearsonr
from sklearn.model_selection import KFold
from sklearn.model_selection import RandomizedSearchCV
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import GridSearchCV
from cubist import Cubist
import joblib
from tqdm import tqdm



import matplotlib.pyplot as plt


# Function to create stratified bins for stratified k-fold
def create_stratified_bins(y, n_bins=5):  # Reduced number of bins
    bins = np.linspace(min(y), max(y), n_bins)
    y_binned = np.digitize(y, bins=bins)
    return y_binned


def plot_true_vs_predicted(true_values, predicted_values, title='True vs Predicted Values'):
    """
    Plot true values vs predicted values with gray circles and a 1:1 red line.

    Parameters:
    - true_values: Array-like, true target values.
    - predicted_values: Array-like, predicted target values.
    - title: String, title of the plot.
    """
    plt.figure(figsize=(8, 8))
    plt.scatter(true_values, predicted_values, c='gray', alpha=0.5, label='Predicted')
    plt.plot([min(true_values), max(true_values)], [min(true_values), max(true_values)], 'r-', label='1:1 Line')
    plt.xlabel('True Values')
    plt.ylabel('Predicted Values')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.show()




def preprocess_data(file_path):
    """
    Preprocesses a CSV file containing wet chemistry data.

    Parameters:
        file_path (str): The path to the CSV file.

    Returns:
        pandas.DataFrame: The preprocessed DataFrame.

    This function reads the CSV file located at the specified file path, performs several preprocessing steps,
    and returns the preprocessed DataFrame. The preprocessing steps include:
    - Dropping the 'EC' column due to many samples with NaN values.
    - Replacing negative values in the 'Corg' column with zero.
    - Converting specified columns to numeric data type, coercing non-numeric values to NaN.
    - Renaming the '<2Micron' column to '2Micron'.
    - Dropping rows with NaN values.

    Example usage:
    file_path = "...."
    output = preprocess_data(file_path)
    print(output.head())
    """
    # Read the CSV file
    output = pd.read_csv(file_path, delimiter=";", decimal=",")

    # Drop 'EC' column due to many samples with NaN values
    output.drop('EC', axis=1, inplace=True)

    # Replace negative values in 'Corg' column with zero
    output.loc[output['Corg'] == '-0.1', 'Corg'] = 0
    output['Corg'] = pd.to_numeric(output['Corg'], errors='coerce')

    # Convert specified columns to numeric, coercing non-numeric values to NaN
    output.iloc[:, 2:11] = output.iloc[:, 2:11].apply(pd.to_numeric, errors='coerce')

    # Convert 'Nt' column to numeric, coercing non-numeric values to NaN
    output['Nt'] = pd.to_numeric(output['Nt'], errors='coerce')

    # Convert 'MWD' column to numeric, coercing non-numeric values to NaN
    output['MWD'] = pd.to_numeric(output['MWD'], errors='coerce')

    # Drop rows with NaN values
    output.dropna(inplace=True)

    return output


import pandas as pd


def preprocess_data_EC(file_path):
    """
    Preprocesses a CSV file containing wet chemistry data.

    Parameters:
        file_path (str): The path to the CSV file.

    Returns:
        pandas.DataFrame: The preprocessed DataFrame.

    This function reads the CSV file located at the specified file path, performs several preprocessing steps,
    and returns the preprocessed DataFrame. The preprocessing steps include:
    - Keeping only the 'EC' column.
    - Dropping rows with NaN values in any column.
    - Converting the 'EC' column to numeric data type, coercing non-numeric values to NaN.

    Example usage:
    file_path = "path_to_csv_file.csv"
    output = preprocess_data(file_path)
    print(output.head())
    """
    # Read the CSV file
    output = pd.read_csv(file_path, delimiter=";", decimal=",")

    # Keep only the 'EC' column
    output = output[['Lab_ID','EC']]

    # Convert 'EC' column to numeric, coercing non-numeric values to NaN
    output['EC'] = pd.to_numeric(output['EC'], errors='coerce')

    # Filter out rows where 'EC' is greater than 400
    output = output[output['EC'] <= 400]


    # Drop rows with NaN values in any column
    output.dropna(inplace=True)

    return output

# Example usage:
# file_path = "path_to_csv_file.csv"
# output = preprocess_data(file_path)
# print(output.head())




def rmsle(y_true, y_pred):
    return np.sqrt(np.mean(np.square(np.log1p(y_pred) - np.log1p(y_true))))

def concordance_cc(y_true, y_pred):
    s_xy = np.cov(y_true, y_pred)[0, 1]
    s_x = np.var(y_true)
    s_y = np.var(y_pred)
    mean_x = np.mean(y_true)
    mean_y = np.mean(y_pred)
    return (2 * s_xy) / (s_x + s_y + (mean_x - mean_y) ** 2)

class LSTMModel(nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim, num_layers):
        """
        Initialize the LSTM model.

        Args:
            input_dim (int): The number of expected features in the input x.
            hidden_dim (int): The number of features in the hidden state h.
            num_layers (int): Number of recurrent layers. E.g., setting num_layers=2 would mean stacking
                two LSTMs together to form a `stacked LSTM`, with the second LSTM taking in outputs of the
                first LSTM and producing the final results.
            output_dim (int): The size of the output.

        """
        super(LSTMModel, self).__init__()
        self.hidden_dim = hidden_dim
        self.num_layers = num_layers
        self.lstm = nn.LSTM(input_dim, hidden_dim, num_layers, batch_first=True)
        # -> x needs to be: (batch_size, seq, input_size)

        self.fc = nn.Linear(hidden_dim, output_dim)

    def forward(self, x):
        """
        Forward propagate the LSTM model.

        Args:
            x (torch.Tensor): Input tensor of shape (batch_size, seq_length, input_size)

        Returns:
            torch.Tensor: Output tensor of shape (batch_size, output_size)

        """

        # Set initial hidden states and cell states for LSTM | the deufault is zero so we could skip this step
        h0 = torch.zeros(self.num_layers, x.size(0), self.hidden_dim).to(x.device) # x: (n, seq, input)
        c0 = torch.zeros(self.num_layers, x.size(0), self.hidden_dim).to(x.device) # h0 & c0: (lstm_layers, n, hiden_dim)

        # Forward propagate LSTM
        out, _ = self.lstm(x, (h0, c0)) # out: tensor of shape (batch_size, seq_length, hidden_dim)

        out = self.fc(out[:, -1, :])  # Taking the output of the last time step  

        return out