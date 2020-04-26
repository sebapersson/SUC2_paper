# About 

This is the GitHub repository for the manuscript *Work in progress* [1]. This repository contains all the code required for reproducing the modelling results presented in the manuscript. 

The aim of the modelling was to elucidate the mechanism behind recently collected single-cell time-lapse data of the *SUC2*-gene expression upon long term glucose starvation (Fig. 1) in *Saccharomyces cerevisiae*. To model the single-cell time-lapse data in an ODE-framework, two methods were used and compare. The standard two stage (STS) approach [2], and non-linear mixed effects (NLME) modelling [3]. What causes the decrease in *SUC2*-expression (Fig. 1) is currently unknown. In the manuscript we propose, by combining dynamic modelling and experimental data, a mechanism for this decrease.  


![fig:1](https://i.imgur.com/GfE4uA5.jpg)
<font size="1.5"> **Figure 1:** *Single cell SUC2 gene expression measured via YFP in *S. cerevisiae* when the external glucose is reduced from 4 % to 0.1 % at time zero. The black lines corresponds to the observed 0.05, 0.5 and 0.95 quantiles, while the blue line are four randomly chosen cells.* </font>

To fully replicate the result presented in the manuscript the requirements in the **Requirements for replication of result** section should be fulfilled. Given this, the result should be reproducible by running the *Run_all.sh*-script, which can be found in the Code directory. 

## Repository structure

Each directory contains a README.md describing the role of that directory. These roles can be summarised as:

* **Data**: Contains the original SUC2 data in CSV-format. Upon cloning the project the data is not in the folder, but it can be retreived from [ADD].
* **Intermediate**: Contains intermediate which files are produced by the analysis. Files here are not counted as results, but are required for producing the figures and tables in the result directory. 
* **Code**: Contains the code for processing the data, running STS and NLME, and plotting the end result. The Monolix files, e.g mlxtran-files, are also located here. 
* **Result**: Contains, after running *Run_all.sh*, the figures presented in the manuscript. It also contains tables with estimated parameter values for the NLME and STS approaches. 

## Requirements for reproducing the result

This section contains information of operating system, programming languages (with libraries) and software required to reproduce the result. 

### Operating system 


This repository was created on a computer running on Ubuntu 18.04.4 (Linux). In the code, all file-paths were encoded as relative paths using Unix-paths (e.g `dir_data = ./../Data/`). The run-all script is a shell-script (bash). Hence, the code should be able to run on any Unix-based system (e.g Linux, Mac). However, problems can arise with the file-paths and shell-script on Windows. One way to resolve this for a Windows user might be (I have not tested this) to download the [Ubuntu-terminal](https://ubuntu.com/tutorials/tutorial-ubuntu-on-windows#1-overview) for Windows. If this fails, a Windows user can change the paths in the scripts and runs the files in the order described in the README-file in the Code-directory. 

### Programming languaes (and libraries)

Three programming languages are required for reproducing the results, Julia (version 1.3.1), Python (version 3.7.4) and R (version 3.6.3)

**[Julia](https://julialang.org/)** (version 1.3.1) [4] was used for STS approach, and for simulating cells given estimated population parameters. As Julia is language under heavy development, reproducing the results might fail if not the **exact** same Julia and library versions are used as those listed here. The libraries used are: 

* **CSV** (version 0.5.24)
* **Distributions** (version 0.22.4)
* **NLopt** (version 0.5.1)
* **ProgressMeter** (version 1.2.0)
* **DifferentialEquations** (version 6.11.0)
* **DataFrames** (version 0.20.2)

**[R](https://www.r-project.org/)** (version 3.6.3) was used for processing the data and plotting. Mostly the tidyverse was used for this [5], and as only standard functions were used, the result should be reproducible if newer (but also older to a limit) versions of R and the listed libraries are used. The libraries used are: 

* **tidyverse** (version 1.3.0)
* **stringr** (version 1.4.0)
* **latex2exp** (version 0.4.0)
* **ggthemes** (version 4.2.0)
* **readxl** (version 1.3.1)

**[Python](https://www.python.org/)** (version 3.7.4) was used to process parts of the Monolix result. As only standard functions were used, the result should be reproducible if newer (but also older to a limit) versions of Python and the listed libraries are used. The libraries used are: 
* **numpy** (version 1.17.2)
* **pandas** (version 0.25.1)

### Software

**[Monolix](http://lixoft.com/products/monolix/)** (version 2019R2) [6] was used for the NLME-parameter estimation. I do not know how much changes between versions, hence I recommend using the **exact** same version if attempting to reproduce the results. When running the parameter estimation, the default seed in Monolix was used. 


## References 
1. Work in progress 
2. Karlsson M, Janzén DL, Durrieu L, Colman-Lerner A, Kjellsson MC, Cedersund G. Nonlinear mixed-effects modelling for single cell estimation: when, why, and how to use it. BMC systemsbiology. 2015;9(1):52.
3. Davidian M, Giltinan DM. Nonlinear models for repeated measurement data: an overview and update. Journal of agricultural, biological, and environmental statistics. 2003;8(4):387.
4. Bezanson J, Edelman A, Karpinski S, Shah VB. Julia: A Fresh Approach to Numerical Computing. SIAM. 2017;59(1):65–98.
5. Grolemund G, Wickham H. R for data science. 2018.
6. Monolix version 2019R2. Antony, France: Lixoft SAS; 2019. http://lixoft.com/products/monolix/.
