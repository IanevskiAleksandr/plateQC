## Description
`plateQC` is a comprehensive R package designed to calculate plate quality control metrics for high-throughput experiments. It includes a novel metric called NRFE (Normalized Residual Fit Error) that is based on dose-response data, rather than relying solely on controls (such as DMSO). The package also computes other common metrics like Z-factor, SSMD, and signal-to-background ratio, which are typically calculated using control data.

### Key Features:
- Example data included in the package
- Advanced quality control metrics calculation (NRFE, Z-factor, SSMD, signal-to-background ratio)
- Interactive plate visualizations including:
  - Inhibition heatmaps
  - Residual error heatmaps
  - Row-wise distribution plots
- Parallel processing support for large datasets
- Robust outlier detection and handling

## Installation

You can install the package from GitHub:

```r
# Install devtools if not already installed
if (!require("devtools")) install.packages("devtools")

# Install plateQC
devtools::install_github("IanevskiAleksandr/plateQC")
```

### Data Format

The input data frame should contain the following columns:
- `BARCODE`: Unique identifier for each plate
- `DRUG_NAME`: Name of the drug or control (use "<b>POS_CTRL</b>" and "<b>NEG_CTRL</b>" for positive and negative controls, respectively). A negative control would mean that it will not have any impact on the cells while positive control is the treatment which will have maximum response. DMSO is a classical example of a negative control (that does not affect cell viability). A toxic compound like Benzethonium chloride(BzCl) which is a potent proteosome inhibitor and kills all the cells in the well is used as a positive control. 
- `CONC`: Drug concentration in nM
- `INTENSITY`: Measured response intensity
- `WELL`: Well position identifier (e.g., "A1", "B2"). Optional column for the core analysis; Required only for generating plate visualizations. If not provided, make sure to run without visualization `process_plate_data(plate_data, plot_dose_response = FALSE, plot_plate_summary = FALSE)`

## Usage

### Basic Usage

```r
library(plateQC)

# View example data
head(plate_data) 

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

### Quality Metrics

The package calculates several quality metrics:
- **NRFE**: Normalized Residual Fit Error calculated based on normalized dose-response curve fitting residuals
- **Z-factor**: Classical plate quality metric based on controls
- **SSMD**: Strictly Standardized Mean Difference
- **Robust Z-prime**: Robust version of Z-factor using median and MAD
- **Signal vs Background**: Ratio between positive and negative controls

## Visualization Examples

The package generates four types of visualizations:
1. Plate heatmaps showing inhibition patterns
2. Residual error heatmaps highlighting problematic areas
3. Scatter plots showing inhibition patterns
4. Dose-response curve plots for each compound

## Contributing

Contributions are welcome! Please feel free to submit pull requests or create issues for bugs and feature requests.

## Citation

If you use this package in your research, please cite:
```
Will be added later
GitHub repository: https://github.com/IanevskiAleksandr/plateQC
```

## License

MIT License
