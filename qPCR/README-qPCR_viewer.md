# Interactive qPCR Data Visualization App

This is a Shiny application written in R for the interactive visualization of qPCR data from Quantstudio Excel exports.

## Features

- **Four Interactive Plots**:
    1.  Melt Curve: Fluorescence vs. Temperature
    2.  Melt Curve: Derivative vs. Temperature
    3.  Amplification: Rn vs. Cycle
    4.  Amplification: Delta Rn vs. Cycle
---

## Setup and Installation

This application is designed to run within a self-contained `conda` environment to ensure all dependencies are met.

### Prerequisites

- You must have [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Miniforge](https://www.anaconda.com/products/distribution](https://conda-forge.org/download/) installed on your system.

### Installation Steps

1.  **Open your terminal** (e.g., Command Prompt, PowerShell, or Terminal on macOS/Linux).

2.  **Create a new conda environment** for this project. This isolates the required packages from your other projects. We will name it `qpcr_viz`.

    ```bash
    conda create --name lab_tools
    ```

3.  **Activate the new environment**. You must do this every time you want to run the application.

    ```bash
    conda activate lab_tools
    ```

4.  **Install R and all required R packages** into the environment from the `conda-forge` channel. This single command installs everything you need.

    ```bash
    conda install -c conda-forge r-base r-shiny r-plotly r-dplyr r-readr r-bslib r-scales r-tidyr r-readxl
    ```

Your environment is now set up and ready to use.

---

## How to Run the Application

1.  **Activate the conda environment** (if you haven't already):
    ```bash
    conda activate lab_tools
    ```

2.  **Navigate to the directory** where you have saved the `app.R` script.

3.  **Run the script** using the `Rscript` command, followed by the path to your qPCR Excel file. **It is important to enclose the file path in quotes**, especially if it contains spaces.

    **Example on Windows:**
    ```bash
    Rscript app.R "C:\Users\YourUser\Documents\qPCR Data\my_experiment.xls"
    ```

    **Example on macOS/Linux:**
    ```bash
    Rscript app.R "/Users/youruser/Documents/qPCR Data/my_experiment.xls"
    ```

The script will launch, and the interactive application will automatically open in your default web browser.

---

## Using the Application

The interface is designed to be intuitive:

- **Filtering**: Use the checkboxes in the sidebar on the left to control which samples are displayed. You can filter by **Target Name**, **Sample ID**, and **Well ID**. Use the "Select All" and "Deselect All" buttons for quick selections.
- **Plot Interaction**: Hover your mouse over any line on any of the four plots to see a detailed tooltip with all relevant information for that specific sample.
- **Zooming and Panning**:
    - Click and drag on a plot to zoom into a specific region.
    - Your zoomed view will be remembered even if you change the filters.
    - **To reset a plot's view, simply double-click anywhere on that plot.**
