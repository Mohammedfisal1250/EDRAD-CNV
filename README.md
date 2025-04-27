EDRAD-CNV
EDRAD-CNV: A new approach for Copy Number Variation (CNV) detection from next-generation sequencing data, combining global density and local entropy strategies, with an intuitive Tkinter GUI interface for simplified interaction.

Usage
Run the GUI Application

(1) Open the main_EDRAD.py file (or the main script you are using).

(2) Configure important parameters inside the GUI after launching:

Upload input data (.bam files)

Set parameters: k (number of neighbors), bandwidth, bin size, reference genome path, etc.

Choose segmentation method (Python or R-based)

Click the "Run" button to start the CNV detection analysis!

(3) Visualize the results directly inside the GUI:

Metrics like Precision, Sensitivity, and F1-score

Graphical display of detected CNVs

Quick Example

python
Copy
Edit
# Open Anaconda prompt or terminal
python main_EDRAD.py
Then interact with the GUI window to upload your files and run detection easily!

Features
Easy upload of BAM files through file dialog

User-configurable parameters for flexible analysis

Supports both Python and R-based segmentation

Visualizes results with scores and CNV plots

Beginner-friendly, no command-line skills needed

GUI Screenshot

Required Dependencies
Python 3.8 (recommended)

numpy

pandas

scipy

biopython

pysam

pyod

sklearn

rpy2

tkinter (comes preinstalled with Python)

R 3.4.4 (for CBS segmentation if selected)

DNAcopy

Real Datasets
You can test EDRAD-CNV using public datasets like:

1000 Genomes Project

Your own BAM files prepared for analysis

Citing EDRAD-CNV
If you use EDRAD-CNV in your research or publication, please cite this project appropriately. (Add your future paper info here if you publish.)

