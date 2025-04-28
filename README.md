# EDRAD-CNV
EDRAD-CNV: A new approach for Copy Number Variation (CNV) detection from next-generation sequencing data, combining global density and local entropy strategies, with an intuitive Tkinter GUI interface for simplified interaction.

# Usage
### How to Run EDRAD-CNV
(1) There are 2 options:
#### Option 1: Run from Python Script
Open the `EDRAD-CNV-Main.py` file and execute it directly using your Python environment.

#### Option 2: Run as an Executable Program (Linux)
Alternatively, you can run the pre-built executable version by following these steps:
- This software is available on request.
```bash
cd /mnt/c/EDRAD-CNV
./dist/EDRAD-CNV-Main
```
(2) After the GUI appears, as you see in Screenshot.png upload the required input files:

- Upload the input .bam file

- Upload the reference genome .fa or .fasta file

- Upload the ground truth file

(3) Set important parameters inside the GUI with default settings:

- k (number of neighbors) = 20

- bandwidth = 1

- bin size = 1000

- The column size used for segmentation = 50

(4) Choose the segmentation method (Python-based or R-based);

(5) for the output files must do this before clich on Run Analysis button:
- CNV Output File = output_file.txt

- P-value Output File = P-value.txt

- Result File ""this is the same CNV Output File for calculating the scores" = output_file.txt

- Score Result File = Score_Result.txt

(6) Click the Run Analysis button to start the CNV detection analysis!

(7) After the analysis completes, check the output files:

- Results file containing detected CNVs

- Evaluation metrics file (Precision, Sensitivity, F1-score)

(8) Visualize the results directly inside the GUI:

- Metrics such as Precision, Sensitivity, and F1-score

- Graphical display of the detected CNVs

Then interact with the GUI window to upload your files and run detection easily!

# Features
- Easy upload of BAM files through file dialog

- User-configurable parameters for flexible analysis

- Supports both Python and R-based segmentation

- Visualizes results with scores and CNV plots

- Beginner-friendly, no command-line skills needed

# GUI Screenshot
  - See the Screenshot.png

# Required Dependencies

1. Python 3.8            
   - biopython     
   - combo         
   - numpy         
   - pandas        
   - pysam        
   - pyod         
   - rpy2          
   - scikit-learn  
   - scipy         
   - tkinter (comes preinstalled with Python)
2. R 3.4.4
   - DNAcopy
     
# Real Datasets

The real datasets can be obtained in the following 2 ways.

- clink this link：https://pan.baidu.com/s/1Ja4XH2wZupeAcwc9qhZn8A extraction code：29to
- [1000 Genomes Project](https://www.internationalgenome.org/)

# Citing EDRAD-CNV

# For more details, please contact us via email:
- mohammedfisal18@gmail.com
