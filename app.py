"""
This application is an extention to Mplus that adds the ability to calculate reliablities and t-tests.
"""

import toga
from toga.style import Pack
from toga.style.pack import COLUMN
import pandas as pd     # A Pandas dataframe is used for the data
import re               # Regular expressions are used to read the Mplus input file
from scipy import stats # SciPy is used to calculate the t-tests
from math import sqrt   # The square root is used to calculate Cohen's d

saved = False   # Global variable that remembers if the output is saved

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
    # Open the .inp file and read its content
    with open(file_path, 'r') as f:
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
            data_file_line = next(line_iter).strip()  # Move to the next line after 'DATA:'
            match = re.search(r'FILE\s+IS\s+"(.+?)";', data_file_line, re.IGNORECASE)
            if match:
                data_file_path = match.group(1).strip()
            else:
                return

        # Find the variable names in the VARIABLES section
        if line.startswith('VARIABLE:'):
            variable_lines = [] # Initialise list
            while True:
                line = next(line_iter).strip()
                if line.endswith(';'):
                    variable_lines.append(line[:-1])  # Exclude the trailing ';'
                    break
                variable_lines.append(line)
            # Join the lines and split into variable names
            variables = ' '.join(variable_lines).split()
            variables = variables[2:]

        # Find the missing value specification, handling multi-line input
        if line.startswith('MISSING ARE '):
            missing_lines = []
            while True:
                if line.endswith(';'):
                    missing_lines.append(line[:-1])  # Exclude the trailing ';'
                    break
                missing_lines.append(line)
                line = next(line_iter).strip()

            # Join the lines and parse the missing values specification
            missing_spec = ' '.join(missing_lines)
            missing_values_matches = re.findall(r'(\w+|\w+-\w+)\((\d+)\)', missing_spec)

            # Populate missing_values dictionary
            for var_range, missing_code in missing_values_matches:
                if '-' in var_range:
                    # Handle ranges like bk1-be4
                    start_var, end_var = var_range.split('-')
                    start_idx = variables.index(start_var)
                    end_idx = variables.index(end_var)
                    for var in variables[start_idx:end_idx+1]:
                        missing_values[var] = int(missing_code)
                else:
                    missing_values[var_range] = int(missing_code)

        # Find the latent variable specifications, handling multi-line input
        if line.startswith('MODEL:'):
             # Regular expression pattern to match the lines of interest
            latent_pattern = re.compile(r'(\w+)\s+BY\s+([\w\s]+);')
            while True:
                line = next(line_iter).strip()
                if line == '':
                    break
                # Find all matches for the latent variables and their indicators
                match = latent_pattern.search(line)
                if match:
                    latent_variable = match.group(1)  # Latent variable name
                    indicators = match.group(2).split()  # List of indicators
                    model_dict[latent_variable] = indicators

        # Find the latent variable pairs for t-tests, handling multi-line input
        if line.startswith('T-TESTS:'):
             # Regular expression pattern to match the lines of interest
            latent_pattern = re.compile(r'(\w+)\s+WITH\s+(\w+);')
            while True:
                line = next(line_iter).strip()
                if line == '':
                    break
                # Find all matches for the latent variables and their indicators
                match = latent_pattern.search(line)
                if match:
                    latent_variable_1 = match.group(1)  # First latent variable name
                    latent_variable_2 = match.group(2)  # Second latent variable name
                    t_tests.append((latent_variable_1,latent_variable_2))

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


def t_test_analysis(df: pd.DataFrame, items_a: list[str], items_b: list[str], name_a: str, name_b: str) -> str:
    """This function conducts a t-test analysis for two paired latent variables A and B, e.g. repeated measures of the same variable

    Args:
        df (pd.DataFrame): The data
        items_a (list[str]): The items of the latent variable A 
        items_b (list[str]): The items of the latent variable B
        name_a (str): The name of the latent variable A
        name_b (str): The name of the latent variable B

    Returns:
        str: output - The results in plain text format
    """    
    output = '-------------------------------------\n'
    output += f't-test for {name_a} vs. {name_b}:\n\n'
    # 1. Mean and SD of latent variable A
    mean_a, sd_a, n_a = mean_sd(df, items_a)
    output += f'Mean of {name_a}: {mean_a:.3f}, SD: {sd_a:.3f}; N = {n_a}'

    # 2. Mean and SD of latent variable A for students who have data in latent variable B
    mean_a_participants, sd_a_participants, n_a_participants = mean_sd_filtered(df, items_a, items_b)
    output += f'Mean of {name_a}* (for {name_b} participants): {mean_a_participants:.3f}, SD: {sd_a_participants:.3f}; N = {n_a_participants}'
    # 2.1 Pooled standard deviation
    pooled_sd = sqrt(((n_a - 1) * (sd_a ** 2) + (n_a_participants - 1) * (sd_a_participants ** 2)) / (n_a + n_a_participants - 2))

    # 3. Mean and SD of latent variable B
    mean_b, sd_b, n_b = mean_sd(df, items_b)
    output += f'Mean of {name_b}: {mean_b:.3f}, SD: {sd_b:.3f}; N = {n_b}\n'

    # 4. t-test for the mean differences
    # 4.1 for A and A*
    mean_diff_a, t_stat_a, df_value_a, p_value_a = t_test_independent(df, items_a, items_b)
    effect_size_ind = abs((mean_a - mean_a_participants) / pooled_sd)
    # 4.2 for A and B
    mean_diff, t_stat, df_value, p_value, effect_size_paired = t_test_paired(df, items_a, items_b)
        
    output += f'Mean difference ({name_a} - {name_a}*): {mean_diff_a:.3f}\n'
    output += f't-test: t({df_value_a}) = {t_stat_a:.3f}, p = {p_value_a:.3f}\n'
    output += f'Effect size: {effect_size_ind:.2f}, {interpret_cohens_d(effect_size_ind)} effect\n\n'
    output += f'Mean difference ({name_b} - {name_a}): {mean_diff:.3f}\n'
    output += f't-test: t({df_value}) = {t_stat:.3f}, p = {p_value:.3f}\n'
    output += f'Effect size: {effect_size_paired:.2f}, {interpret_cohens_d(effect_size_paired)} effect\n'

    return output


class MplusT(toga.App):
    def startup(self):
        """Construct and show the Toga application.
        """
        main_box = toga.Box(style=Pack(direction=COLUMN, padding=10))

        self.main_window = toga.MainWindow(title=self.formal_name)
        self.main_window.content = main_box
        self.main_window.show()
        
        # Add open and save file dialogs to the menu
        open_cmd = toga.Command(
            self.open_file_dialog,
            text='Open',
            tooltip='Open an Mplus input file',
            shortcut=toga.Key.MOD_1 + 'o',
            group=toga.Group.FILE
        )
        save_cmd = toga.Command(
            self.save_to_file,
            text='Save',
            tooltip='Save the output to a file',
            shortcut=toga.Key.MOD_1 + 's',
            group=toga.Group.FILE
        )
        self.commands.add(open_cmd)
        self.commands.add(save_cmd)

        # Create a button to open the file dialog
        open_file_button = toga.Button(
            "Open an Mplus input file",
            on_press=self.open_file_dialog,
            style=Pack(padding=10),
        )
        main_box.add(open_file_button)
        
        # Create the output document window (MultilineTextInput)
        self.output_window = toga.MultilineTextInput(
            readonly=True, style=Pack(flex=1, padding=10)
        )
        main_box.add(self.output_window)
        
        # Create a button for saving the output to a file
        save_button = toga.Button(
            "Save output to file",
            on_press=self.save_to_file,
            style=Pack(padding=10),
        )
        main_box.add(save_button)


    async def open_file_dialog(self, widget=None):
        """This dialog opens the Mplus input file. After the file is opend, it does all the main work.

        Args:
            widget (None): This argument is required by Toga, but it is not used.
        """        
        # Open a file dialog using the async method
        inp_file = await self.main_window.dialog(
            toga.OpenFileDialog(title='Select an Mplus input file', file_types=['inp'], multiple_select=False)
        )

        # Extract the information given in the input file
        # 1. data_file_path (str): Path and name of the data file as provided in the FILE IS section of the input file
        # 2. variables (list(str)): List of all variable names in the dataset as provided in the NAMES ARE section of the input file
        # 3. missing_values (dict): Dictionary that contains codes for missing values as provided in the MISSING ARE section of the input file
        # 4. model_dict (dict): Dictionary that defines the latent variables as provided in the MODEL: section of the input file
        #                       latent_var BY var1 var2 var3 becomes {"latent_var": ["var1", "var2", "var3"]}
        # 5. t_tests (list(list(str))): List of pairs of variable names for the t-test from the T-TESTS: section of the input file
        #                               var1 WITH var2 becomes [["var1", "var2"]]
        data_file_path, variables, missing_values, model_dict, t_tests = read_mplus_inp(inp_file)


        # Check the result and show a dialog accordingly
        if data_file_path:
            await self.main_window.dialog(
                toga.InfoDialog('File Selected', f'The data file is: {data_file_path}')
            )

            # Load the data file based on the extracted info
            df = read_mplus_data(data_file_path, variables, missing_values)

            # DOING THE ANALYSES
            # Create a dictionary to hold subsets of the original DataFrame for each latent variable
            latent_variable_data = {}

            # Loop through the dictionary of latent variables to create subsets and calculate reliablities
            for latent_variable, indicators in model_dict.items():
                latent_variable_data[latent_variable] = df[indicators]

                # Calculate Cronbach's alpha for each latent variable 
                alpha_value = cronbach_alpha(latent_variable_data[latent_variable])

                self.output_window.value += f"Cronbach's alpha for {latent_variable}: {alpha_value:.3f}\n"

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
                self.output_window.value += t_test_analysis(df, start_items, end_items, var_a, var_b)
        else:
            await self.main_window.dialog(
                toga.InfoDialog('No Data File Found', 'Data file path not found. Please check the input file format.')
            )

    async def save_to_file(self, widget=None):
        """_summary_

        Args:
            widget (None): This argument is required by Toga, but it is not used.
        """        
        global saved
        # Open a save file dialog
        save_path = await self.main_window.dialog(
            toga.SaveFileDialog(title='Save Output', suggested_filename='output.txt')
        )

        if save_path:
            # Write the content of the document window to the selected file
            try:
                with open(save_path, 'w', encoding='utf-8') as file:
                    file.write(self.output_window.value)
                await self.main_window.dialog(
                    toga.InfoDialog('Save Successful', f'Output saved to {save_path}')
                )
                saved = True
            except Exception as e:
                await self.main_window.dialog(
                    toga.ErrorDialog('Save Failed', f'An error occurred: {e}')
                )
                saved = False
        else:
            await self.main_window.dialog(
                toga.InfoDialog('Save Cancelled', 'File save operation was cancelled.')
            )
            saved = False

    async def on_exit(self):
        global saved
        # Check if there's content in the document window
        if self.output_window.value.strip() and not saved:  # If there's unsaved output
            save_choice = await self.main_window.dialog(
                toga.QuestionDialog(
                    title='Unsaved Changes',
                    message='You have unsaved output. Do you want to save it before exiting?',
                )
            )

            if save_choice:  # If the user chooses to save
                await self.save_to_file()
                if not saved:  # If save was unsuccessful or canceled
                    return False  # Cancel app exit
        return True  # Allow app exit

def main():
    return MplusT()
