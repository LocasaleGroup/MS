## C13 abundance correction
1. download the following files or folder:
  - aply_pynac.py
  - apply_config.txt
  - pynac
2. change the settings in *apply_config.txt* file
3. run: python apply_pynac.py apply_config.txt

**Note**: The default C13 abundance is 0.01109, users can change this value at line 30 of this file: pynac/core/isotopelabels.py 

## calculate percentage of metabolites
1. download the *calculate_percentage.py* and *calculate_percentage_conf.txt*
2. change the settings in *calculate_percentage_conf.txt*
3. run: python *calculate_percentage.py*

## Data

The *inputFile_test.xlsx* can be used by above programs