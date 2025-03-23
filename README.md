# R Template for Meta-Analysis

## Overview
This repository contains an R template designed to enable scientists to conduct meta-analyses on their data. The template is built to be versatile and user-friendly, allowing for application to continuous outcome variables to estimate effects across user-defined blocking variables.

## Files
- `r-template-meta-analysis.Rmd`: The main R Markdown template
- `functions.R`: Helper functions for data analysis and visualization
- `study.csv`: Example dataset for demonstration

## Requirements
- R (>= 4.0.0 recommended)
- RStudio
- Required R packages are installed automatically through the `renv` package

## How to Use
1. Clone or download this repository
2. Open the project in RStudio
3. Run `renv::restore()` in the console to set up the correct package environment
4. Modify the parameters in the YAML header of the R Markdown file:
   - `file_name`: Name of your data file
   - `block`: Column name representing your block of interest
   - `ref_trt`: String representing the reference group
   - `trt_interest`: Strings representing treatments of interest
   - `measure`: Column representing the outcome variable
   - `y_axis_string`: Label for the graph representing the outcome variable
5. Click the "Knit" button to generate the report

## Data Format Requirements
- Data should be in long format with a single header
- Required columns: outcome variable, block, treatment
- Extra columns will be ignored
- Missing values should be specified as empty cells or NA without quotes
- The outcome variable should be a continuous quantitative variable
- The block variable should have a number of defined levels (e.g., donor1, donor2, donor3)

## Features
- Boxplot visualizations of data by blocking variable
- Linear regression models and treatment contrasts
- Linear mixed regression models
- Diagnostic checks for model validity
- Random effects meta-analysis
- Forest plots for visualization
- Heterogeneity assessment


## Contact
[e-mail](gabbyteku@gmail.com)
