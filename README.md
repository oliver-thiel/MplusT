# MplusT
MplusT is an app for researchers who use Mplus by Muthén and Muthén for data analysis. Mplus is great for many cases but does not calculate reliabilities and t-tests. This app adds these features.

## Purpose and features
The app will be used in addition to Mplus. For more information about Mplus, see https://www.statmodel.com/index.shtml
However, you do not need Mplus to run this application. It can be used independently to calculate reliabilities and t-tests from a dataset in CSV format.
The application is built using the Beeware Toga library, which allows for cross-platform GUI development.
It uses Pandas to handle data, regular expressions to parse Mplus input files, and SciPy to perform statistical tests.
The application has a simple GUI with options to load data, calculate reliabilities, perform t-tests, and save the output. The output is displayed in a text box within the application.

## How to use it
Prepare the data and the input file before you run the app.

### The Mplus data file
The app uses the Pandas function pd.read_csv(data_file_path, sep=',', header=None, names=variables) to read the Mplus data. So, please make sure that the data file has the correct format. It should be a CSV file that uses a comma as the separator. The file should not contain a header. The variable names are provided in the Mplus input file.
The path and filename of the data file are specified in the Mplus input file.

### The Mplus input file
The Mplus input is a text file. The app reads the whole file but uses only a few parts. The other parts are ignored. You can choose if you want to use the whole Mplus input file or only the parts that the app needs. The app extracts information from the following sections:

```
DATA:
  FILE IS "C:\Users\oth\Documents\Mplus\data.dat";
```
This section provides the path and name of the Mplus data file. Version 1.1 reads the file name correctly even if DATA: FILE IS is a single line. FILE = is not supported.

```
VARIABLE:
  NAMES ARE var1 var2 var3
            var4 var5 var6 var7;
```
This section provides the variable names. It may have several lines. The app reads all names until it reaches the semicolon.

```
  USEVARIABLES ARE var1-var3 var4 var5 var6;
```
This section defines which variables are used.

```
  MISSING ARE var1-var3(9) var4(99) var5(99);
```
This section defines the codes of the missing values.

```
MODEL:  factor1 BY var1-var4; ! A latent variable
        factor2 BY var3 var5 var6; ! A second factor
```
This section provides information about the latent variables. The app MplusT will calculate Cronbach's alpha for each latent variable specified in this section.

```
T-TESTS: factor1 WITH factor2;
         var1 WITH var2;
```
This section does not belong to the original Mplus input file. It has to be added. It provides information about the t-tests that shall be conducted. The app can compare both latent and observed variables. For each pair of variables factor1 WITH factor2, it conducts two t-tests:
1. It compares the mean of factor1 for all observed cases with the mean of factor1 for only those cases that have data in factor2.
2. It compares the mean of factor1 with the mean of factor2.

For each t-test, it calculates the means *M*, standard deviations *SD*, sample size *N*, mean difference *&#916;M*, *t*-value, degrees of freedom *df*, significance level (*p*-value), and effect size (Cohen's *d*).

### Installing and running the app on Windows
You can install the app under Windows by downloading the installer.

After the app is installed, you start it like any Windows application.

A window opens:

![MplusTscreenshot](https://github.com/user-attachments/assets/b410e751-53f7-4a4c-b50b-954c1e74552d)

When you click on [Open an Mplus input file], you can choose the input file you will use.
The app calculates the statistics specified in the file immediately after you open the input file.
Click on [Save output to file] to save the output.

### Running the app using Python
If you do not use the Windows app but run MplusT.py with Python, ensure all dependencies are installed. The app requires
* Pandas https://pandas.pydata.org/
* SciPy https://scipy.org/

Furthermore, you must specify the path of the input file. The default filename is `input.inp`. If your file has another name or a different path, you must edit line 328 of the file MplusT.py.

```
    # Define the filename of the Mplus input file (if needed including the path, e.g. 'C:\\Users\\oth\\Documents\\Mplus\\input_file.inp')
    inp_file = 'input_file.inp'
```

## The output
The app prints the output to the terminal. It looks like this:

```
Cronbach's alpha for factor1: 0.893
Cronbach's alpha for item var1 removed: 0.878
Cronbach's alpha for item var2 removed: 0.876
Cronbach's alpha for item var3 removed: 0.870
Cronbach's alpha for item var4 removed: 0.856
Removing any item would not increase Cronbach's alpha.

Cronbach's alpha for factor2: 0.904
Cronbach's alpha for item var3 removed: 0.889
Cronbach's alpha for item var5 removed: 0.911
Cronbach's alpha for item var6 removed: 0.880
Removing item var5 would increase Cronbach's alpha to 0.911


-------------------------------------
t-test for factor1 vs. factor2:

Mean of factor1: 4.805, SD: 1.141; N = 339
Mean of factor1* (for factor2 participants): 4.856, SD: 1.120; N = 195
Mean of factor2: 4.560, SD: 1.106; N = 199

Mean difference (factor1 - factor1*): -0.052
t-test: t(532) = -0.511, p = 0.609
Effect size: 0.05, very small effect

Mean difference (factor2 - factor1): -0.290
t-test: t(194) = -3.680, p = 0.000
Effect size: 0.26, small effect
-------------------------------------
t-test for var1 vs. var2:

Mean of var1: 4.652, SD: 1.161; N = 341
Mean of var1* (for var2 participants): 4.727, SD: 1.137; N = 203
Mean of var2: 4.492, SD: 1.201; N = 207

Mean difference (var1 - var1*): -0.074
t-test: t(542) = -0.729, p = 0.466
Effect size: 0.06, very small effect

Mean difference (var2 - var1): -0.241
t-test: t(202) = -2.740, p = 0.007
Effect size: 0.19, small effect
```
