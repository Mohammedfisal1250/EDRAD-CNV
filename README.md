# EDRAD-CNV
EDRAD-CNV: A new approach for Copy Number Variation (CNV) detection from next-generation sequencing data, combining global density and local entropy strategies, with an intuitive Tkinter GUI interface for simplified interaction.

# Usage
Run the GUI Application

(1) Open the main_EDRAD.py file (or the main script you are using).

(2) Configure important parameters inside the GUI after launching:

- Upload input data (.bam files)

- Set parameters: k (number of neighbors), bandwidth, bin size, reference genome path, etc.

- Choose segmentation method (Python or R-based)

- Click the "Run" button to start the CNV detection analysis!

(3) Visualize the results directly inside the GUI:

- Metrics like Precision, Sensitivity, and F1-score

- Graphical display of detected CNVs

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


