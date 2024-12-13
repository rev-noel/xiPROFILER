# xiPROFILER

**xiPROFILER** is a web-based framework designed to streamline the visualization and analysis of crosslinking complexome profiling (XL-CP) data. By combining crosslinking mass spectrometry (MS) with SDS-PAGE separation, xiPROFILER provides researchers with an intuitive tool to explore protein-protein interactions and complex formation within biological systems.

## Features

- **Interactive Visualization**: Dynamically visualize protein intensity profiles across gel slices.
- **Molecular Weight Calibration**: Convert slice numbers to approximate molecular weights using a user-guided calibration.
- **Data Normalization**: Toggle between raw and normalized data views for detailed analysis.
- **Combined Molecular Weight Display**: Display vertical guide lines for combined molecular weight peaks, highlighting potential crosslinking events.
- **Protein Filtering**: Filter and select proteins using an integrated search and table view.
- **Comparison Across Conditions**: Compare intensity profiles across different experimental conditions.
- **Interpretation Mode**: Highlight known interactors using UniProt information.

## Hosting

xiPROFILER is hosted at: 

## Installation

To use xiPROFILER locally, follow these steps:

### Prerequisites

- Python (>= 3.9)
- Shiny for Python
- Plotly
- Pandas

### Steps

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/xiprofiler.git
   cd xiprofiler
   ```

2. Install the required dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Run the application:
   ```bash
   shiny run app:app
   ```

4. Open your browser and navigate to:
   ```
   http://127.0.0.1:8000
   ```

## Usage

1. **Upload Data**: Load your XL-CP data file in CSV format.
2. **Visualization**:
   - Choose between raw or normalized data.
   - Adjust sliders to calibrate molecular weights.
3. **Analysis**:
   - Use interpretation mode to highlight known protein interactors.
   - Compare data across conditions in the "Comparison" tab.

## Data Requirements

Input data should be a tab-separated values (TSV) file with the following structure:

| Protein.Group | Protein.Ids | Protein.Names | Genes | Slice_1| Slice_2  | ... |
|---------------|-------------|---------------|-------|--------|----------|-----|
| C1P619        | C1P619      | ILVX_ECOLI    | ilvX  | 12345  | 236125   | ... |
| P00350        | P00350      | 6PGD_ECOLI    | gnd   | 345256 | 45734734 | ... |

- **Metadata Columns**:
  - `Protein.Group`, `Protein.Ids`, `Protein.Names`, `Genes`: Describe the proteins in the dataset.
  - `First.Protein.Description` (optional): Provides additional information about the protein.

- **Intensity Columns**:
  - Columns labeled `Slice_1`, `Slice_2`, ..., represent the intensity values for each gel slice or experimental condition.


---

### Data Format Notes:
1. The data must include at least the metadata columns (`Protein.Group`, `Protein.Ids`, `Protein.Names`, `Genes`) and intensity columns.
2. The intensity values should correspond to protein signals across gel slices or experimental conditions.
3. Use tab-separated values (`.tsv`) format for input files.

## Known Limitations

- Organism-specific proteome mass ranges may limit molecular weight calculations.
- High crosslinker concentrations may introduce inaccuracies.

## Future Features

- Spike-in experiment support for enhanced molecular weight accuracy.
- Automated interaction prediction using machine learning.

## Contributing

Contributions are welcome! Please follow these steps:

1. Fork the repository.
2. Create a new branch for your feature/bugfix.
3. Submit a pull request.

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.





