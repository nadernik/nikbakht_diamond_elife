# Rat Vision Analysis - Python Implementation

This directory contains a Python implementation of the rat vision analysis from Nikbakht, N., & Diamond, M. E. (2021). Conserved visual capacity of rats under red light. eLife, 10, e66429.

## Overview

The code analyzes behavioral and spectral data from rat vision experiments under different LED light conditions. It generates publication-ready figures showing:

- Psychometric curves for orientation discrimination
- Daily performance across training sessions
- Performance distributions
- Spectral power density and photon flux calculations
- Retinal sensitivity functions

## Files

- `beautify.py` - Professional plotting utilities
- `utils.py` - Helper functions for data loading and analysis
- `plot_behavior.py` - Behavioral data analysis and visualization
- `plot_spect.py` - Spectral analysis and photon calculations
- `main.py` - Command-line interface to run analyses
- `requirements.txt` - Python package dependencies

## Installation

1. Install Python dependencies:
```bash
pip install -r requirements.txt
```

2. Ensure the MATLAB data files are in the parent directory:
   - `behavior_data.mat`
   - `spectral_data.mat`
   - `plos_s_labda.mat`

## Usage

### Command Line Interface

Run analyses using the main script:

```bash
# Run behavioral analysis only
python main.py --behavior

# Run spectral analysis only
python main.py --spectral

# Run both analyses
python main.py --all

# Save figures to files instead of displaying
python main.py --behavior --save-figs

# Specify output directory for saved figures
python main.py --all --save-figs --output-dir my_figures
```

### Individual Scripts

You can also run the analysis scripts directly:

```bash
# Behavioral analysis
python plot_behavior.py

# Spectral analysis
python plot_spect.py
```

## Output

The scripts generate several figures:

### Behavioral Analysis (`plot_behavior.py`)
- **Figure 1**: Psychometric curves for each LED condition (2x3 subplot)
- **Figure 2**: Daily performance across sessions
- **Figure 3**: Average performance distributions
- **Figure 4**: Cumulative performance comparison

### Spectral Analysis (`plot_spect.py`)
- **Figure 1**: Spectral power density plots
- **Figure 2**: Spectral photon count plots
- **Figure 3**: Normalized spectra with opsin sensitivity curves
- **Figure 4**: Retinal spectral sensitivity function

## Data Requirements

The scripts expect MATLAB `.mat` files with specific structures:

### `behavior_data.mat`
Contains behavioral data with fields:
- `DATA.results.percent_v` - Proportion of "vertical" responses
- `DATA.results.fit` - Fitted psychometric parameters
- `DATA.perfinday` - Daily performance data
- `DATA.stats` - Statistical summaries

### `spectral_data.mat`
Contains spectral measurements:
- `spectrum_raw` - Raw spectral data [wavelength, intensity]
- `led` - LED condition labels

### `plos_s_labda.mat`
Contains retinal sensitivity data:
- `s_lambda.wl` - Wavelength values
- `s_lambda.measured` - Measured sensitivity
- `s_lambda.fit_fft` - Fitted sensitivity curve

## Dependencies

- **numpy** - Numerical computing
- **scipy** - Scientific computing and optimization
- **matplotlib** - Plotting and visualization
- **seaborn** - Statistical visualization
- **pandas** - Data manipulation
- **h5py** - HDF5 file support (for newer .mat files)

## Notes

- The code uses scipy.constants for physical constants (speed of light, Planck's constant)
- Figures are formatted using the `beautify()` function for publication-ready appearance
- Error handling is included for robust data loading and processing
- The implementation maintains the same mathematical formulations as the original MATLAB code

## Citation

If you use this code, please cite:
Nikbakht, N., & Diamond, M. E. (2021). Conserved visual capacity of rats under red light. eLife, 10, e66429. 