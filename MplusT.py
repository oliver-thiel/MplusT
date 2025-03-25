import pandas as pd     # A Pandas dataframe is used for the data
import re               # Regular expressions are used to read the Mplus input file
import os
from scipy import stats # SciPy is used to calculate the t-tests
from math import sqrt   # The square root is used to calculate Cohen's d


def read_mplus_inp(file_path: str):
    """Reads an Mplus input file (.inp) and extract key information

    Args:
        file_path (str): filename including path

    Returns:
        str: data_file_path - A string that contains the name of the data file including its path
        list[str]: variables, - A list of all variable names in the data set
        dict{str : int}: missing_values - A dictonary that contains all variables
                                          that have specific codes for missing data.
                                          Defaults to {} if there are specicfic codes for missing values. 
        dict{str : list[str]}: model_dict - A dictonary that associates each latent variable with its items.
        list[list[str]]: t_tests - A list of pairs of variable names
    """    
    # Open the .inp file and read its content if it exsists. Otherwise return None and empty lists/dicts.
    if not os.path.exists(file_path):
        print(f"File '{file_path}' does not exist!")
        return None, [], {}, {}, []
    
    # Open the .inp file and read its content
    with open(file_path, 'r') as f:
        print(f'============ Reading input file: {file_path} ============')
        lines = f.readlines()

    # Convert the list of lines into an iterator
    line_iter = iter(lines)

    # Initialize variables to store extracted information
    data_file_path = None
    variables = []
    missing_values = {}
    model_dict = {}
    t_tests = []

    # Loop through the lines to find DATA, VARIABLE, MISSING, MODEL, and T-TESTS information
    for line in line_iter:
        line = line.strip()

        # Find the data file path (adjusted for "FILE IS" instead of "FILE =")
        if line.startswith('DATA:'):
            data_file_line = line[5:] # Remove 'DATA:'
            if data_file_line == '':
                data_file_line = next(line_iter, None)  # Move to the next line after 'DATA:'
            if data_file_line is None: # Check for unexpected end of file
                return data_file_path, variables, missing_values, model_dict, t_tests
            match = re.search(r'FILE\s+IS\s+"(.+?)";', data_file_line.strip(), re.IGNORECASE)
            if match:
                data_file_path = match.group(1).strip()
            else: # If no match found for data file path, return empty values and exit function early
                print('ERROR: No data file path found in DATA section.')
                return data_file_path, variables, missing_values, model_dict, t_tests

        # Find the variable names in the VARIABLES section
        if line.startswith('VARIABLE:'):
            line = line[9:].strip()  # Remove 'VARIABLE:' and strip whitespace
            variable_lines = [] # Initialise list
            while True:
                if line.endswith(';'):
                    variable_lines.append(line[:-1])  # Exclude the trailing ';'
                    break
                variable_lines.append(line)
                line = next(line_iter, None)
                if line is None:
                    break
                else:
                    line = line.strip() # Remove leading/trailing whitespace from each line

            # Join the lines and split into variable names
            variables = ' '.join(variable_lines).split()
            variables = variables[2:] # Remove the first two elements which are 'NAMES ARE'

        # Find the used variables, handling multi-line input
        if line.startswith('USEVARIABLES ARE '):
            used_lines = []
            while True:
                line = line.strip() # Remove leading/trailing whitespace
                if line.endswith(';'):
                    used_lines.append(line[:-1])  # Exclude the trailing ';'
                    break
                used_lines.append(line)
                line = next(line_iter, None)
                if line is None:
                    break

            # Join the lines and parse the used variable specification
            used_spec = ' '.join(used_lines).split()
            used_spec = used_spec[2:] # Remove the first two elements which are 'USEVARIABLES ARE'

            used_variables = [] # Handle ranges like bk1-be4
            for used in used_spec:
                if '-' in used:
                    # Handle ranges like bk1-be4
                    start_var, end_var = used.split('-')
                    try:
                        start_idx = variables.index(start_var)
                    except ValueError:
                        print(f'ERROR: Variable {start_var} not found in the list of variables.')
                        continue
                    try:
                        end_idx = variables.index(end_var)
                    except ValueError:
                        print(f'ERROR: Variable {end_var} not found in the list of variables.')
                        continue
                    for var in variables[start_idx:end_idx+1]:
                        used_variables.append(var)
                else:
                    used_variables.append(used)

        # Find the missing value specification, handling multi-line input
        if line.startswith('MISSING ARE '):
            missing_lines = []
            while True:
                if line.endswith(';'):
                    missing_lines.append(line[:-1])  # Exclude the trailing ';'
                    break
                missing_lines.append(line)
                line = next(line_iter, None)
                if line is None:
                    break
                else:
                    line = line.strip() # Remove leading/trailing whitespace

            # Join the lines and parse the missing values specification
            missing_spec = ' '.join(missing_lines)
            missing_values_matches = re.findall(r'(\w+|\w+-\w+)\((\d+)\)', missing_spec)

            # Populate missing_values dictionary
            for var_range, missing_code in missing_values_matches:
                if '-' in var_range:
                    # Handle ranges like bk1-be4
                    start_var, end_var = var_range.split('-')
                    try:
                        start_idx = used_variables.index(start_var)
                    except ValueError:
                        print(f'ERROR: Variable {start_var} not found in the list of used variables.')
                        continue
                    try:
                        end_idx = used_variables.index(end_var)
                    except ValueError:
                        print(f'ERROR: Variable {end_var} not found in the list of used variables.')
                        continue
                    for var in used_variables[start_idx:end_idx+1]:
                        missing_values[var] = int(missing_code)
                else:
                    missing_values[var_range] = int(missing_code)

        # Find the latent variable specifications, handling multi-line input
        if line.startswith('MODEL:'):
            line = line[6:].strip()  # Remove 'MODEL:' and strip whitespaces
             # Regular expression pattern to match the lines of interest
            latent_pattern = re.compile(r'(\w+)\s+BY\s+([\w\s]+|\w+-\w+);')
            while True:
                # Find all matches for the latent variables and their indicators
                match = latent_pattern.search(line)
                if match:
                    latent_variable = match.group(1)  # Latent variable name
                    indicators = match.group(2).split()  # List of indicators
                    corrected = [] # Handle ranges like bk1-be4
                    for ind in indicators:
                        if '-' in ind:
                            # Handle ranges like bk1-be4
                            start_var, end_var = ind.split('-')
                            try:
                                start_idx = used_variables.index(start_var)
                            except ValueError:
                                print(f'ERROR: Variable {start_var} not found in the list of used variables.')
                                continue
                            try:
                                end_idx = used_variables.index(end_var)
                            except ValueError:
                                print(f'ERROR: Variable {end_var} not found in the list of used variables.')
                                continue
                            for var in used_variables[start_idx:end_idx+1]:
                                corrected.append(var)
                        else:
                            corrected.append(ind)
                    model_dict[latent_variable] = corrected
                line = next(line_iter, None)
                if line is None:
                    break
                else:
                    line = line.strip() # Remove leading/trailing whitespace
                if line == '':
                    break

        # Find the latent variable pairs for t-tests, handling multi-line input
        if line.startswith('T-TESTS:'):
            line = line[8:].strip() # Remove 'T-TESTS:' and strip whitespace
             # Regular expression pattern to match the lines of interest
            latent_pattern = re.compile(r'(\w+)\s+WITH\s+(\w+);')
            while True:
                # Find all matches for the latent variables and their indicators
                match = latent_pattern.search(line)
                if match:
                    latent_variable_1 = match.group(1)  # First latent variable name
                    latent_variable_2 = match.group(2)  # Second latent variable name
                    t_tests.append((latent_variable_1,latent_variable_2))
                line = next(line_iter, None)
                if line is None:
                    break
                else:
                    line = line.strip() # Remove leading/trailing whitespace
                if line == '':
                    break

    return data_file_path, variables, missing_values, model_dict, t_tests


def read_mplus_data(data_file_path: str, variables: list[str], missing_values: dict|None=None):
    """Reads data from a Mplus data file (.dat) and returns it as a Pandas dataframe

    Args:
        data_file_path (str): The name of the data file including its path
        variables (list[str]): A list of all variable names
        missing_values (dict{str : int}, optional): A dictonary that contains all variables
                                                    that have specific codes for missing data.
                                                    Defaults to None if there are specicfic codes for missing values.

    Returns:
        pd.dataframe: df - A Pandas data frame the contains all data
    """    
    # Read the .dat file into a DataFrame using a comma as the delimiter
    df = pd.read_csv(data_file_path, sep=',', header=None, names=variables)

    # Handle missing values if specified
    if missing_values:
        for var, missing_code in missing_values.items():
            # Use 'mask' instead of 'replace' to avoid recursion issues
            df[var] = df[var].mask(df[var] == missing_code, pd.NA)

    return df


def cronbach_alpha(df_items: pd.DataFrame) -> float:
    """Calculates the internal reliability (Cronbach's alpha) of the scale in a given data frame

    Args:
        df_items (pd.DataFrame): The data frame the contains the scale

    Returns:
        float: The calculated value of Cronbach's alpha
    """    
    # Drop rows with missing values to avoid issues with variance calculation
    df_items_clean = df_items.dropna()

    # Number of items
    item_count = df_items_clean.shape[1]
    
    # Variance of each item
    item_variances = df_items_clean.var(axis=0, ddof=1)
    
    # Variance of the total score (sum of all items)
    total_score_variance = df_items_clean.sum(axis=1).var(ddof=1)
    
    # Cronbach's alpha formula
    alpha = (item_count / (item_count - 1)) * (1 - (item_variances.sum() / total_score_variance))
    
    return alpha


def mean_sd(df: pd.DataFrame, items: list[str]):
    """Calculates the mean and standard deviation of a scale that is given by the selected items

    Args:
        df (pd.dataframe): The data
        items (list[str]): A list of the selected items

    Returns:
        float: The calculated mean
        float: The calculated standard deviation
        int: The number of participants
    """    
    # Drop rows with missing values in the selected items
    df_clean = df[items].dropna()
    means = df_clean.mean(axis=1)
    return means.mean(), means.std(), len(means)


def mean_sd_filtered(df: pd.DataFrame, items: list[str], filter: list[str]):
    """Calculates the mean and standard deviation of a scale that is given by the selected items.
    This function does not consider all cases but a subset that does not have missing values in the filter variables.

    Args:
        df (pd.dataframe): The data
        items (list[str]): A list of the selected items
        filter (list[str]): A list of the selected items that function as a filter

    Returns the output of the function mean_sd():
        float: The calculated mean
        float: The calculated standard deviation
        int: The number of participants
    """    
    # Only consider rows where there are no missing values in the filter
    df_clean = df.dropna(subset=filter)
    return mean_sd(df_clean, items)


def interpret_cohens_d(cohens_d):
    """
    Determines text interpretation of effect size given Cohen's d value
    (Adapted from Dan Friedman: https://dfrieds.com/math/effect-size.html#many-examples-interpretation-of-cohens-d-sizes)

    Args:
        cohens_d (float): Cohen's d value
    
    Returns:
        str: adjective to describe the magnitude of the effect
    """
    if cohens_d < 0:
        cohens_d *= -1 # The value should be positive

    if cohens_d < 0.1:
        return "very small"
    if cohens_d < 0.35:
        return "small"
    if cohens_d < 0.65:
        return "medium"
    if cohens_d < 0.9:
        return "large"
    return "very large"


def t_test_paired(df: pd.DataFrame, items_a: list[str], items_b: list[str]):
    """Calculates the if the mean difference between two latent variables A and B is significant

    Args:
        df (pd.dataframe): The data
        items_a (list[str]): The items of the latent variable A
        items_b (list[str]): The items of the latent variable B

    Returns:
        float: Mean difference between A and B
        float: t-value for Student's for paired samples
        int: Degrees of freedom
        float: p_value
        float: Effect size Cohen's d
    """     
    df_clean = df.dropna(subset=items_a + items_b)  # Ensure complete data for both time points
    mean_a = df_clean[items_a].mean(axis=1)
    mean_b = df_clean[items_b].mean(axis=1)
    diff = mean_b - mean_a
    mean_diff = diff.mean()
    t_stat, p_value = stats.ttest_rel(mean_b, mean_a)
    cohens_d = abs((mean_a.mean() - mean_b.mean())/diff.std())
    return mean_diff, t_stat, len(mean_a) - 1, p_value, cohens_d


def t_test_independent(df: pd.DataFrame, items_a: list[str], items_b: list[str]):
    """Calculates if the mean difference between a latent variable in two independent subsamples is significant

    Args:
        df (pd.dataframe): The data
        items_a (list[str]): The items of the latent variable A
        items_b (list[str]): The items of the filter variable

    Returns:
        float: mean_diff_a - Mean difference between A and filtered A
        float: t_stat_a - t-value for Welch's t-test for independent samples
        int: degrees of freedom
        float: p_value_a
    """     
    df_a = df[items_a].dropna()
    mean_a_all = df_a.mean(axis=1)
    df_clean = df.dropna(subset=items_a + items_b)  # Ensure complete data for both time points
    mean_a = df_clean[items_a].mean(axis=1)
    mean_diff_a = mean_a_all.mean() - mean_a.mean()
    t_stat_a, p_value_a = stats.ttest_ind(mean_a_all, mean_a, equal_var=False) # Welch's t-test
    return mean_diff_a, t_stat_a, len(mean_a_all) + len(mean_a) - 2, p_value_a


def t_test_analysis(df: pd.DataFrame, items_a: list[str], items_b: list[str], name_a: str, name_b: str) -> None:
    """This function conducts a t-test analysis for two paired latent variables A and B, e.g. repeated measures of the same variable

    Args:
        df (pd.DataFrame): The data
        items_a (list[str]): The items of the latent variable A 
        items_b (list[str]): The items of the latent variable B
        name_a (str): The name of the latent variable A
        name_b (str): The name of the latent variable B
    """    
    print('-------------------------------------')
    print(f't-test for {name_a} vs. {name_b}:\n')
    # 1. Mean and SD of latent variable A
    mean_a, sd_a, n_a = mean_sd(df, items_a)
    print(f'Mean of {name_a}: {mean_a:.3f}, SD: {sd_a:.3f}; N = {n_a}')

    # 2. Mean and SD of latent variable A for students who have data in latent variable B
    mean_a_participants, sd_a_participants, n_a_participants = mean_sd_filtered(df, items_a, items_b)
    print(f'Mean of {name_a}* (for {name_b} participants): {mean_a_participants:.3f}, SD: {sd_a_participants:.3f}; N = {n_a_participants}')
    # 2.1 Pooled standard deviation
    pooled_sd = sqrt(((n_a - 1) * (sd_a ** 2) + (n_a_participants - 1) * (sd_a_participants ** 2)) / (n_a + n_a_participants - 2))

    # 3. Mean and SD of latent variable B
    mean_b, sd_b, n_b = mean_sd(df, items_b)
    print(f'Mean of {name_b}: {mean_b:.3f}, SD: {sd_b:.3f}; N = {n_b}\n')

    # 4. t-test for the mean differences
    # 4.1 for A and A*
    mean_diff_a, t_stat_a, df_value_a, p_value_a = t_test_independent(df, items_a, items_b)
    effect_size_ind = abs((mean_a - mean_a_participants) / pooled_sd)
    # 4.2 for A and B
    mean_diff, t_stat, df_value, p_value, effect_size_paired = t_test_paired(df, items_a, items_b)
        
    print(f'Mean difference ({name_a} - {name_a}*): {mean_diff_a:.3f}')
    print(f't-test: t({df_value_a}) = {t_stat_a:.3f}, p = {p_value_a:.3f}')
    print(f'Effect size: {effect_size_ind:.2f}, {interpret_cohens_d(effect_size_ind)} effect\n')
    print(f'Mean difference ({name_b} - {name_a}): {mean_diff:.3f}')
    print(f't-test: t({df_value}) = {t_stat:.3f}, p = {p_value:.3f}')
    print(f'Effect size: {effect_size_paired:.2f}, {interpret_cohens_d(effect_size_paired)} effect')


# MAIN PROGRAM
if __name__ == '__main__':
    # READ THE INPUT AND THE DATA FILE
    # Define the filename of the Mplus input file (if needed including the path, e.g. 'C:\\Users\\oth\\Documents\\Mplus\\input_file.inp')
    inp_file = 'input_file.inp'

    # Extract the information given in the input file
    # 1. data_file_path (str): Path and name of the data file as provided in the FILE IS section of the input file
    # 2. variables (list(str)): List of all variable names in the dataset as provided in the NAMES ARE section of the input file
    # 3. missing_values (dict): Dictionary that contains codes for missing values as provided in the MISSING ARE section of the input file
    # 4. model_dict (dict): Dictionary that defines the latent variables as provided in the MODEL: section of the input file
    #                       latent_var BY var1 var2 var3 becomes {"latent_var": ["var1", "var2", "var3"]}
    # 5. t_tests (list(list(str))): List of pairs of variable names for the t-test from the T-TESTS: section of the input file
    #                               var1 WITH var2 becomes [["var1", "var2"]]
    data_file_path, variables, missing_values, model_dict, t_tests = read_mplus_inp(inp_file)

    # Check if the data file path was found
    if data_file_path:
        # Load the data file based on the extracted info
        df = read_mplus_data(data_file_path, variables, missing_values)

        # DOING THE ANALYSES
        # Create a dictionary to hold subsets of the original DataFrame for each latent variable
        latent_variable_data = {}

        # Loop through the dictionary of latent variables to create subsets and calculate reliablities
        for latent_variable, indicators in model_dict.items():
            latent_variable_data[latent_variable] = df[indicators]

            alpha_values = [] # List of values to find maximum

            # Calculate Cronbach's alpha for each latent variable 
            alpha_value = cronbach_alpha(latent_variable_data[latent_variable])
            alpha_values.append(alpha_value)  # Append alpha value to the list

            print(f"Cronbach's alpha for {latent_variable}: {alpha_value:.3f}")
            for indicator in indicators:
                without = latent_variable_data[latent_variable].drop(columns=[indicator])
                without_alpha = cronbach_alpha(without)  # Cronbach's alpha without the current item
                alpha_values.append(without_alpha)  # Append alpha value to the list

                print(f"Cronbach's alpha for item {indicator} removed: {without_alpha:.3f}")
                
            max_alpha_value = max(alpha_values)  # Find the maximum alpha value in the list
            # Check if the maximum alpha value is not the original alpha value
            if max_alpha_value != alpha_value:
                # Find the index of the maximum alpha value in the list
                max_alpha_index = alpha_values.index(max_alpha_value)
                # Get the indicator that was removed to get the maximum alpha value
                removed_indicator = indicators[max_alpha_index - 1]
                print(f"Removing item {removed_indicator} would increase Cronbach's alpha to {max_alpha_value:.3f}")
            else: # If no item removal increases alpha, print this message
                print("Removing any item would not increase Cronbach's alpha.")  

        # Loop through the list of t-test pairs and perform the tests
        for var_a, var_b in t_tests:
            if var_a in model_dict:
                start_items = model_dict[var_a]
            else:
                start_items = [var_a]
            if var_b in model_dict:
                end_items = model_dict[var_b]
            else:
                end_items = [var_b]
            t_test_analysis(df, start_items, end_items, var_a, var_b)
    else:
        print('Data file path not found. Please check the input file format.')
