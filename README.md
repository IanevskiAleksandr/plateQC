# plateQC

[![R-CMD-check](https://img.shields.io/badge/R%20CMD%20check-passing-brightgreen)](https://github.com/IanevskiAleksandr/plateQC)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

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
- **Comprehensive quality control workflow** with decision frameworks

## Quality Metrics

The package calculates several quality metrics:
- **NRFE**: Normalized Residual Fit Error calculated based on normalized dose-response curve fitting residuals
- **Z-factor**: Classical plate quality metric based on controls
- **SSMD**: Strictly Standardized Mean Difference
- **Robust Z-prime**: Robust version of Z-factor using median and MAD
- **Signal vs Background**: Ratio between positive and negative controls

| Metric | Calculation | Interpretation |
|--------|-------------|----------------|
| **Z-factor** | `1 - (3σ_pos + 3σ_neg)/|μ_pos - μ_neg|` | >0.5: excellent, 0.3-0.5: acceptable |
| **SSMD** | `(μ_neg - μ_pos)/√(σ²_neg + σ²_pos)` | >2: good separation |
| **S/B** | `μ_neg/μ_pos` | >5: adequate dynamic range |
| **NRFE** | Mean normalized residuals from dose-response fits | <10: good spatial quality |

## Installation

You can install the package from GitHub:

```r
# Install devtools if not already installed
if (!require("devtools")) install.packages("devtools")

# Install plateQC
devtools::install_github("IanevskiAleksandr/plateQC")
```

## Data Format

The input data frame should contain the following columns:
- `BARCODE`: Unique identifier for each plate
- `DRUG_NAME`: Name of the drug or control (use "**POS_CTRL**" and "**NEG_CTRL**" for positive and negative controls, respectively). A negative control would mean that it will not have any impact on the cells while positive control is the treatment which will have maximum response. DMSO is a classical example of a negative control (that does not affect cell viability). A toxic compound like Benzethonium chloride(BzCl) which is a potent proteosome inhibitor and kills all the cells in the well is used as a positive control. 
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

<br><hr><br>

## Quality Control Workflow

### Step 1: Calculate Quality Metrics

```r
library(plateQC)

# Process your plate data
results <- process_plate_data(
  plate_data,
  plot_dose_response = FALSE,    # Set TRUE to generate dose-response curves
  plot_plate_summary = FALSE,    # Set TRUE for plate visualizations
  verbose = TRUE                 # Show processing messages
)

# Extract quality metrics
quality_metrics <- results$plate_statistics
print(quality_metrics)
```

### Step 2: Interpret Quality Metrics

#### Traditional Control-Based Metrics:
- **Z-factor**: >0.5 (excellent), 0.3-0.5 (acceptable), <0.3 (poor)
- **SSMD**: >3 (excellent), 2-3 (good), 1-2 (acceptable), <1 (poor)  
- **S/B ratio**: >10 (excellent), 5-10 (good), 2-5 (acceptable), <2 (poor)

#### NRFE (Normalized Residual Fit Error):
- **<10**: Acceptable spatial quality
- **10-15**: Borderline - requires review
- **>15**: Poor spatial quality - likely systematic artifacts

### Step 3: Quality Decision Framework

```r
# Classify plate quality
library(dplyr)  # Required for case_when function

classify_plate_quality <- function(quality_metrics) {
  quality_metrics$quality_category <- with(quality_metrics, {
    case_when(
      # Excellent quality
      zfactor > 0.5 & NRFE < 10 ~ "Excellent",
      
      # Good quality  
      zfactor > 0.3 & NRFE < 10 ~ "Good",
      
      # Spatial artifacts detected (like GDSC1 plate 101416)
      zfactor > 0.3 & NRFE > 15 ~ "Spatial_Artifacts",
      
      # Control issues but spatial quality OK
      zfactor < 0.3 & NRFE < 10 ~ "Control_Issues", 
      
      # Borderline cases
      zfactor > 0.3 & NRFE >= 10 & NRFE <= 15 ~ "Borderline",
      
      # Poor overall quality
      TRUE ~ "Poor"
    )
  })
  
  quality_metrics$recommended_action <- with(quality_metrics, {
    case_when(
      quality_category == "Excellent" ~ "Accept",
      quality_category == "Good" ~ "Accept", 
      quality_category == "Spatial_Artifacts" ~ "Exclude - systematic errors",
      quality_category == "Control_Issues" ~ "Flag - manual review of controls",
      quality_category == "Borderline" ~ "Flag - visual inspection required",
      quality_category == "Poor" ~ "Exclude"
    )
  })
  
  return(quality_metrics)
}

# Apply classification to example data
classified_results <- classify_plate_quality(results$plate_statistics)
print(classified_results[, c("barcode", "zfactor", "NRFE", "quality_category", "recommended_action")])

#      barcode zfactor      NRFE    quality_category           recommended_action
# 1 PSCREEN_M2   -1.28 23.989170                Poor                      Exclude
# 2     115858    0.88  4.999876           Excellent                       Accept
# 3     101416    0.58 26.530984   Spatial_Artifacts Exclude - systematic errors
```

**Analysis of the three example plates:**

**Plate 115858** - Excellent Quality:
- Z-factor: 0.88 (EXCELLENT - above 0.5)
- SSMD: 25 (EXCELLENT - above 2.0) 
- NRFE: 5.0 (EXCELLENT - below 10)
- **Verdict: ACCEPT** - High confidence data

**Plate 101416** - Spatial Artifacts:
- Z-factor: 0.58 (PASS - above 0.5)
- SSMD: 7 (PASS - above 2.0)  
- NRFE: 26.5 (FAIL - systematic spatial artifacts)
- **Verdict: EXCLUDE** - Traditional metrics pass but NRFE detects spatial problems

**Plate PSCREEN_M2** - Poor Controls:
- Z-factor: -1.28 (FAIL - poor control separation)
- SSMD: 1 (FAIL - below 2.0)
- NRFE: 24.0 (FAIL - spatial artifacts + control issues)
- **Verdict: EXCLUDE** - Multiple quality issues

**Key Insight:** Plate 101416 demonstrates NRFE's unique value - traditional metrics suggest acceptable quality, but NRFE correctly identifies systematic spatial artifacts that would compromise dose-response measurements.

### Step 4: Quality Assessment Summary

```r
# Create quality summary for example data
quality_summary <- classified_results %>%
  group_by(quality_category) %>%
  summarise(
    n_plates = n(),
    mean_zfactor = round(mean(zfactor, na.rm = TRUE), 2),
    mean_NRFE = round(mean(NRFE, na.rm = TRUE), 1),
    .groups = 'drop'
  )

print("Quality Distribution:")
print(quality_summary)

# A tibble: 3 × 4
#   quality_category  n_plates mean_zfactor mean_NRFE
#   <chr>                <int>        <dbl>     <dbl>
# 1 Excellent                1         0.88       5.0
# 2 Poor                     1        -1.28      24.0
# 3 Spatial_Artifacts        1         0.58      26.5

# Summary for example data
acceptable_plates <- classified_results %>%
  filter(quality_category %in% c("Excellent", "Good"))

exclude_plates <- classified_results %>%
  filter(quality_category %in% c("Poor", "Spatial_Artifacts"))

cat(sprintf("Quality Assessment Summary for Example Data:
- Excellent/Good quality: %d plates (%.1f%%)
- Recommend exclusion: %d plates (%.1f%%)
", 
nrow(acceptable_plates), 100*nrow(acceptable_plates)/nrow(classified_results),
nrow(exclude_plates), 100*nrow(exclude_plates)/nrow(classified_results)))

# Quality Assessment Summary for Example Data:
# - Excellent/Good quality: 1 plates (33.3%)
# - Recommend exclusion: 2 plates (66.7%)
```

## Visualization Examples

The package generates four types of visualizations:
1. Plate heatmaps showing inhibition patterns
2. Residual error heatmaps highlighting problematic areas
3. Scatter plots showing inhibition patterns
4. Dose-response curve plots for each compound

## Key Advantages of NRFE

1. **Detects spatial artifacts** missed by control-based metrics
2. **Uses actual compound data** rather than just control wells  
3. **Identifies systematic errors** affecting dose-response relationships
4. **Complements traditional metrics** for comprehensive quality assessment
5. **Improves cross-study reproducibility** by filtering problematic plates

## Best Practices

1. **Always combine metrics**: Use NRFE alongside traditional control-based metrics
2. **Visual inspection**: Review borderline plates manually before final decisions
3. **Document decisions**: Keep records of plates excluded and reasons
4. **Assay-specific thresholds**: Adjust quality thresholds based on your experimental setup

## Contributing

Contributions are welcome! Please feel free to submit pull requests or create issues for bugs and feature requests.

## Citation

If you use this package in your research, please cite:
```
Ianevski A, et al. "Spatial artifact detection improves reproducibility of drug screening experiments." (2024)
GitHub repository: https://github.com/IanevskiAleksandr/plateQC
```

## License

MIT License
