# EDRAD-CNV
EDRAD-CNV: A new approach for Copy Number Variation (CNV) detection from next-generation sequencing data, combining global density and local entropy strategies, with an intuitive Tkinter GUI interface for simplified interaction.

# Usage
Run the GUI Application

(1) Open the file main_EDRAD.py (or your main script), and run it;

(2) After the GUI appears, upload the required input files:

- Upload the input .bam file

- Upload the reference genome .fa or .fasta file

- Upload the ground truth file

(3) Set important parameters inside the GUI:

- k (number of neighbors)

- bandwidth

- bin size 
- The column size used for segmentation

(4) Choose the segmentation method (Python-based or R-based);

(5) Click the Run Analysis button to start the CNV detection analysis!

(6) After the analysis completes, check the output files:

Results file containing detected CNVs

Evaluation metrics file (Precision, Sensitivity, F1-score)

(7) Visualize the results directly inside the GUI:

Metrics such as Precision, Sensitivity, and F1-score

Graphical display of the detected CNVs

Quick Example

python
Copy
Edit
# Open Anaconda prompt or terminal
```python
python main_EDRAD.py
```
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

Your own BAM files prepared for analysis

# Citing EDRAD-CNV


