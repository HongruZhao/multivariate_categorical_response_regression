# Subspace Decompositions for Association Structure Learning in Multivariate Categorical Response Regression

This repository contains the code for simulation study and real data analysis accompanying the paper **"Subspace decompositions for association structure learning in multivariate categorical response regression"**. <!--The material herein has been prepared for submission to *JRSS-B*.-->

---

## Repository Structure

The repository is organized into two main sections: simulation studies and real data analyses. Both sections employ SLURM for parallel computing; please note that the SLURM scripts are designed for high-performance computing clusters and may execute slowly depending on the computational resources available.

### Simulation Study

- **SLURM Parallel Computing Script**  
  A SLURM job script is provided to run the simulation study in parallel. This script is configured for a typical high-performance cluster environment.  
  **Note:** The code within this file might run slowly due to the heavy computational load.

- **Summary Files**  
  - `Cross_entropy_summary_100rep.rds`: Summarizes the Hellinger distances and misclassification rates of the following estimators:
  - `overlap_summary_100rep.rds`: Summarizes the Hellinger distances and misclassification rates of the following estimators:  
    - **O-Mult**
    - **O-Pois**
    - **L-Mult**
    - **L-Pois**
    - **G-Mult**
    - **G-Pois**
    - **Sep-Mult**
    - **G-Mult-θ**

- **Plotting Script**  
  - `summary_plot_overlap copy.R`: Generates plots based on the output from the `.rds` files mentioned above.

### Real Data Analysis

- **Download Real Data**  
  - `GDS3268_data.R`: Run this script to download the real dataset.

- **SLURM Parallel Computing Script**  
  Similar to the simulation study, a SLURM job script is provided to run the real data analysis in parallel.  
  **Note:** Execution might be slow due to the nature of the computations.

- **Summary Files**  
  - `sep_summary_1000rep.rds`: Summarizes the empirical cross entropy and misclassification rates.
  - `Cross_entropy_summary_1000rep.rds`: Summarizes the empirical cross entropy and misclassification rates.
  - `overlap_summary_1000rep.rds`:Summarizes the empirical cross entropy and misclassification rates.

  The estimators considered in the real data analysis are:
    - **O-Mult**
    - **O-Pois**
    - **L-Mult**
    - **L-Pois**
    - **G-Mult**
    - **G-Pois**
    - **Sep-Mult**
    - **G-Mult-θ**

- **Table Generation Script**  
  - `final_output.R`: Processes the `.rds` files and generates tables for the real data analysis.

---

## Quick Start

To quickly generate the main outputs, follow these steps:

### Simulation Study
Run the plotting script to generate plots from the simulation study results:
```bash
Rscript "summary_plot_overlap copy.R"
```

### Real Data Analysis
Run the table generation script to create tables from the real data analysis results:
```bash
Rscript "final_output.R"
```


