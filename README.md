# Upstroke

Code used to determine upstroke frame number in Garg et al. manuscript. 

The following dependencies must first be installed:
numpy
scipy
peakutils   

Input file must be a csv formatted file with cell names in the first row and values for each cell in the columns. An example data file has been uploaded. 

To run, enter the following at the command line:

1. Process all cells in input file
python upstroke_pipeline2.py <input_file>

2. Process specific cell in input file
python upstroke_pipeline2.py <input file> -c <cell name>

example: python upstroke_pipeline2.py data1.csv -c Cell1

3. Save image after processing file
python upstroke_pipeline2.py <input file> -s



