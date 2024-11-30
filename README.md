# plateQC

High-throughput screening plate quality control with NRFE (Normalized Residual Fit Error) analysis and visualization tools.

## Description

`plateQC` is a comprehensive R package designed for analyzing high-throughput screening (HTS) plate data. It provides robust quality control metrics and visualization tools, with a particular focus on the Normalized Residual Fit Error (NRFE) metric for evaluating dose-response curve quality.

### Key Features:
- Advanced quality control metrics calculation (NRFE, Z-factor, SSMD)
- Automated dose-response curve fitting
- Interactive plate visualizations including:
  - Inhibition heatmaps
  - Residual error heatmaps
  - Row-wise distribution plots
- Parallel processing support for large datasets
- Robust outlier detection and handling

## Installation

You can install the development version from GitHub:

```r
# Install devtools if not already installed
if (!require("devtools")) install.packages("devtools")

# Install plateQC
devtools::install_github("IanevskiAleksandr/plateQC")
```

## Usage

### Basic Usage

```r
library(plateQC)

# Process plate data with default settings
results <- process_plate_data(plate_data)

# Access quality metrics
plate_stats <- results$plate_statistics
print(plate_stats)
```

### Advanced Usage

```r
# Full analysis with visualizations and parallel processing
results <- process_plate_data(
  plate_data,
  plot_dose_response = TRUE,    # Generate dose-response curves
  plot_plate_summary = TRUE,    # Create plate visualization plots
  directory_curves = "./fits",  # Directory for saving plots
  remove_empty_wells = TRUE,    # Remove wells without concentrations
  verbose = TRUE,              # Show processing messages
  cores_n = 4                  # Use 4 cores for parallel processing
)

# Access detailed results
plate_metrics <- results$plate_statistics
processed_data <- results$detailed_results
meta_info <- results$metadata
```

### Data Format

The input data frame should contain the following columns:
- `BARCODE`: Unique identifier for each plate
- `DRUG_NAME`: Name of the drug or control (use "POS_CTRL" and "NEG_CTRL" for controls)
- `CONC`: Drug concentration in nM
- `INTENSITY`: Measured response intensity
- `WELL`: Well position identifier (e.g., "A1", "B2")

### Quality Metrics

The package calculates several quality metrics:
- **NRFE**: Normalized Residual Fit Error for evaluating dose-response curve quality
- **Z-factor**: Classical plate quality metric based on controls
- **SSMD**: Strictly Standardized Mean Difference
- **Robust Z-prime**: Robust version of Z-factor using median and MAD
- **Signal vs Background**: Ratio between positive and negative controls

## Visualization Examples

The package generates three types of visualizations:
1. Dose-response curve plots for each compound
2. Plate heatmaps showing inhibition patterns
3. Residual error heatmaps highlighting problematic areas

## Contributing

Contributions are welcome! Please feel free to submit pull requests or create issues for bugs and feature requests.

## Citation

If you use this package in your research, please cite:
```
Ianevski, A. (2024). plateQC: High-throughput screening plate quality control with NRFE analysis.
GitHub repository: https://github.com/IanevskiAleksandr/plateQC
```

## Author

Aleksandr Ianevski (aleksandr.ianevski@helsinki.fi)

## License

MIT License
