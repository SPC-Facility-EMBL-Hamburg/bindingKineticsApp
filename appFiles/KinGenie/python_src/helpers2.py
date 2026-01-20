import numpy  as np
import pandas as pd
import re

def guess_experiment_name(frd_file):
    # Read the file content
    with open(frd_file, 'r') as file:
        content = file.read()

    # Define the regex pattern
    pattern = r'<ExperimentInfo Name="([^"]+)">'

    # Find the first match
    match = re.search(pattern, content)

    if match:
        experiment_name = match.group(1)  # Get the captured group

        # Remove any sequence of six numbers in a row (we assume it is a timestamp)
        experiment_name = re.sub(r'\b\d{6}\b', '', experiment_name)

        # Remove all extra spaces
        experiment_name = re.sub(r'\s+', ' ', experiment_name).strip()

        return experiment_name
    else:
        return 'Experiment'

def etree_to_dict(tree):
    '''
    Converts xml tree to dictionary
    '''
    tree_dict = {}
    for elem in tree:
        tree_dict[elem.tag] = elem.text
    return tree_dict

def combine_dicts(list_dicts):
    '''
    Combines dictionaries with the same (!) keys
    Outputs dictionary with lists
    '''
    new_dict = {}
    for key in list_dicts[0].keys():
        new_dict[key] = [one_dict[key] for one_dict in list_dicts]
    return new_dict

def guess_experiment_type(files):

    """
    Given a certain file, try to guess if it corresponds to surface-based or solution-based binding experiments

    Args:
        file (str): file name

    Returns:
        str: 'surface' or 'solution'
    """

    for file in files:
        if file.endswith('.frd'):
            return 'surface'

    for file in files:
        if 'ExperimentStep.ini' in file:
            return 'surface'

    for file in files:
        if file.endswith('.csv'):

            # Find if we have a line with the protein concentration
            with open(file, 'r') as f:
                first_line = f.read().splitlines()[0]
                if 'Protein_concentration_micromolar' in first_line:
                    return 'solution'
                else:
                    return 'surface'

    return 'surface'

def median_filter(y,x,rolling_window):

    """

    Compute the median filter of the x vector using a rolling window

    First, we convert the x vector into an integer vector and then
        into time variable to take advantage of pandas function
            rolling().median()

	Returns the y vector passed through the median filter

    """

    scaling_factor = 1e4

    temp_vec     =  np.multiply(x,scaling_factor).astype(int)
    series       =  pd.Series(y,index=temp_vec,dtype=float)
    series.index =  pd.to_datetime(series.index,unit='s')

    roll_window  = str(int(rolling_window*scaling_factor))+"s"

    y_filt = series.rolling(roll_window).median().to_numpy()

    return y_filt

def find_loading_column(df,step_number):
    """
    Given a dataframe with the columns 'Column', 'Type' and 'StepNumber',
    find the loading column before the given step number.
    """

    # Filter the dataframe to get the loading columns before the given step number
    loading_columns = df[(df['Type'] == 'LOADING') & (df['StepNumber'] < step_number)]

    # If there are no loading columns, return None
    if loading_columns.empty:
        return None

    # Get the last loading column before the given step number
    last_loading_column = loading_columns.iloc[-1]

    # Return the column name
    return last_loading_column['Column']