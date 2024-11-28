# Union-Free Generic Depth for Non-Standard Data - Application

## Introduction
This repository contains R-code and data sets corresponding to the "Union-Free Generic Depth for Non-Standard Data (Section 5.)" article. We apply the ufg-depth to two application examples considering spatial-categorical-numerical data and hierarchically ordered data.

The structure of the repository is as follows:
- File _setup_session.R installs all needed R-packages.
- File main_categorical_numerical_spatial.R is the main file to compute the spatial-categorical-numerical data example of Section 5.1.
- File main_hierarchical.R is the main file to compute the hierarchically-ordered data example of Section 5.2.

The code was tested with
- R version 4.2.2
- R version 4.3.2
on
- Linux Ubuntu 20.04.5
- Windows 11

## Setup
1. Please install all necessary R-packages (can be found in _setup_session.R).
2. Download the files main_categorical_numerical_spatial.R and main_hierarchical.R
3. Please download the data set corresponding to hierarchically-ordered data example and store in the same folder as main_hierarchical.R file. This data is freely accessible, but only after registration at the following [online portal](https://search.gesis.org/research_data/ZA5240) (accessed: 08.11.2024). Please download the file ZA5280_v2-0-1.sav (4.40MB) there.

Now, main_categorical_numerical_spatial.R reproduces the results of Section 5.1. and main_hierarchical.R reproduces the results of Section 5.2.

## References to the data sets
- spatial-categorical-numerical data: Baddeley, Adrian, and Rolf Turner. "Spatstat: an R package for analyzing spatial point patterns." Journal of statistical software 12 (2005): 1-42.
- hierarchically ordered data: GESIS. Allgemeine Bevölkerungsumfrage der Sozialwissenschaften Allbus 2014. GESIS Datenarchiv, Köln. ZA5240 Datenfile Version 2.2.0, https://doi.org/10.4232/1.13141, 2018 (last accessed: 28.11.2024)
