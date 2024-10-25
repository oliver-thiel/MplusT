# MplusT
I use Mplus by Muthén and Muthén for data analysis. Mplus is great for many cases, but it does not calculate reliabilties and t-test. This Python program adds these features.

## Purpose and features
The program will be used in addition to Mplus. For more information about Mplus, see https://www.statmodel.com/index.shtml
It reads a Mplus input file and extracts information about the data and the tasks. Then, it reads the Mplus data file and calculates reliabilities (Cronbach's alpha) and t-tests.

## How to use it
Prepare the data and the input file before you run the program.

### The Mplus data file
The program uses the Pandas function pd.read_csv(data_file_path, sep=',', header=None, names=variables) to read the Mplus data. Therefore, you should ensure that the data file has the correct format. It should be a CSV file that uses a comma as the separator. The file should not contain a header. The variable names are provided in the Mplus input file.
The path and filename of the data file are specified in the input file.

### The Mplus input file
The Mplus input is a text file. The program reads the whole file but uses only a few parts. The other parts are ignored. You can choose if you want to use the whole Mplus input file og only the parts that the program needs. The program extracts information from the following sections:

```
DATA:
  FILE IS "C:\Users\oth\Documents\Mplus\data.dat";
```
This section provides the path and name of the Mplus data file.

```
VARIABLE:
  NAMES ARE var1 var2 var3
            var4 var5;
```
This section provides the variable names. It may have several lines. The program reads all names until it reaches the semicolon.

```
  MISSING ARE var1-var3(9) var4(99) var5(99);
```
This section defines the codes of the missing values.

```
MODEL:  
        factor1 BY var1 var2 var5; ! A latent variable
        factor2 BY var3 var4; ! A second factor
```
This section provides information about the latent variables. The program MplusT will calculate Cronbach's alpha for each latent variable specified in this section.

```
T-TESTS:
        factor1 WITH factor2;
        var1 WITH var2;
```
This section does not belong to the original Mplus input file. It has to be added. It provides information about the t-tests that shall be conducted. The program can compare both latent and observed variables. For each pair of variables factor1 WITH factor2, it conducts two t-tests:
1. It compares the mean of factor1 for all observed cases with the mean of factor1 for only those cases that have data in factor2.
2. It compares the mean of factor1 with the mean of factor2.

For each t-test, it calculates the means *M*, standard deviations *SD*, sample size *N*, mean difference *&#916;M*, *t*-value, degrees of freedom *df*, significance level (*p*-value), and effect size (Cohen's *d*).

The default filename of the input file is `input.inp`. If your file has another name or a different path, you must edit line 328 of the file MplusT.py.

```
    # Define the filename of the Mplus input file (if needed including the path, e.g. 'C:\\Users\\oth\\Documents\\Mplus\\input_file.inp')
    inp_file = 'input_file.inp'
```
### Running the program
Before you run the program, ensure that all dependencies are installed. The program requires
* Pandas https://pandas.pydata.org/
* SciPy https://scipy.org/

## The output
The program prints the output to the terminal. It looks like this:

```
Cronbach's alpha for factor1: 0.651
Cronbach's alpha for factor2: 0.904
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
