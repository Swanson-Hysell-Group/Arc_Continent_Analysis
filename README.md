# Arc-continent collisions in the tropics set Earth’s climate state

This repository contains data and code associated with the following paper:

Macdonald, F.M.,  Swanson-Hysell, N.L., Park, Y., Lisiecki, L., and Jagoutz, O. (2019) **Arc-continent collisions in the tropics set Earth’s climate state** *Science* 10.1126/science.aav5300

## Description of code

### Paleogeographic analysis and calculations code

Within the main repository folder is a Jupyter notebook entitled **suture_analysis.ipynb** that contains code that utilizes the suture compilation and the paleogeographic models to conduct the analysis and develop the visualizations associated with the study. This notebook relies on the **recon_tools.py** function library developed in conjunction with this research that is also within the main repository folder. This code, and that of the notebook, relies on the pyGPlates module (https://www.gplates.org/docs/pygplates/) that enables the functionality of the GPlates software package to be programmatically accessed using Python. With the exception of the pyGPlates module, which needs to be installed locally and added to the Python path, the computational environment is specified within the **suture_analysis.yml** file. The **LIP_analysis.ipynb** conducts the LIP area analysis also using pyGPlates and the recon_tools.py functions.

A static html version of the suture_analysis notebook can be viewed here:
[suture_analysis.html](http://htmlpreview.github.io/?https://github.com/Swanson-Hysell-Group/Arc_Continent_Analysis/blob/master/suture_analysis.html?raw=true)


### Statistical tests code

Within the main repository folder is a Matlab script entitled **ice_suture_stats.m** that implements the statistical tests described in the Report and its supplementary materials and develops associated figures.

## data folder

### continental_arcs folder

This folder contains a .csv files with the compilation of continental arc length from Cao et al. (2017).

### ice folder

This folder contains the compilation of latitudinal extent of ice sheets through time.

### sutures folder

This folder contains the compilation of ophiolite-bearing sutures that was constructed for this study and used for the analysis. The main shapefile is Suture_Lines.shp and this shapefile was modified for reconstruction with the paleogeographic models in subfolders. This compilation is visualized within the **suture_analysis.ipynb** Jupyter notebook.

## paleogeo_models folder

This folder contains rotation files and polygon outlines for the paleogeographic models used for the reconstructions. The individual datafiles are described in more detail within **suture_analysis.ipynb**.

## code_output folder

This folder contains figures and tables that are output from the code. The associated visualizations were used to make the figures within the manuscript and supplementary materials.
