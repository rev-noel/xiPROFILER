from shiny import App, Inputs, Outputs, Session, render, ui, reactive, req
import pandas as pd
import matplotlib.pyplot as plt
from shiny.types import FileInfo, ImgData
from scipy.ndimage import gaussian_filter1d
import numpy as np
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import io
from Bio import SeqIO
import shinyswatch
import math
import plotly.graph_objs as go
from pathlib import Path
import re
from typing import List
from shinywidgets import output_widget, render_widget

# Define the main user interface for the xiPROFILER application
app_ui = ui.page_fluid(
    ui.navset_card_tab(
        # Welcome Tab
        ui.nav_panel(
            "Welcome",
            ui.markdown('''
            # Welcome to xiPROFILER

            xiPROFILER is designed to assist in the interpretation of crosslinking data. Developed as part of a Master thesis, it integrates DIA-NN analysis of mass spectrometry (MS) data for enhanced visualization and interpretation.

            ## Key Features

            - **Protein Display Tab**: Visualize pg-matrix data, filter proteins, and activate interpretation mode to display interactors for selected proteins.
            - **Molecular Weight Calculation**: Map slice numbers to molecular weights for easier analysis.
            - **Multiple Condition Comparison (Beta)**: Compare intensity profiles across multiple experiments.

            ## How to Start

            1. Upload your pg-matrix file and select the appropriate organism.
            2. Visualize and filter proteins in the Protein Display tab.
            3. Use the Molecular Weight Calculation tab to parameterize molecular weight mappings.

            Enjoy exploring your data with xiPROFILER!
            '''),
            ui.output_image("image", height="150px"),
            ui.markdown('''
            <details>
              <summary>How to Download a Proteome File from UniProt</summary>
              <p>
                - Visit UniProt and navigate to the Proteomes section.
                - Search for your desired organism (e.g., Homo sapiens).
                - Apply filters to show reviewed entries and customize columns as needed.
                - Download the file in TSV format.
              </p>
            </details>
            ''')
        ),

        # Protein Display Tab
        ui.nav_panel(
            "Protein Display",
            ui.layout_sidebar(
                ui.sidebar(
                    ui.accordion(
                        ui.accordion_panel(
                            "File Upload",
                            ui.input_file("file_upload_pg_matrix", "Upload pg matrix", accept=[".tsv"], multiple=False),
                            ui.input_select(
                                "select_organism",
                                "Select Organism:",
                                {
                                    "e_coli": "Escherichia coli (reviewed)",
                                    "human": "Homo Sapiens (reviewed)",
                                    "b_subtilis": "Bacillus subtilis (reviewed)",
                                    "mouse": "Mouse (reviewed)",
                                    "other": "Other"
                                }
                            ),
                            ui.input_file("file_upload_mol_mass", "Upload Proteome Mass File", accept=[".tsv"], multiple=False),
                        ),
                        ui.accordion_panel(
                            "Plot Mode",
                            ui.input_select(
                                "plot_mode",
                                "Data Processing/ Y-Axis",
                                {
                                    "show_raw": "Raw Data",
                                    "show_raw_log": "Raw Data [Log10]",
                                    "show_raw_bl": "Raw Data + baseline",
                                    "show_median_normalized": "Median Normalized",
                                    "show_median_normalized_bl": "Median Normalized + baseline",
                                    "show_min_max_normalized": "Min-Max Normalized",
                                }
                            ),
                            ui.input_switch("show_interactors", "Interpretation Mode", False),
                            ui.input_radio_buttons(
                                "x_axis",
                                "X-Axis",
                                {
                                    "slice_axis": "Slice Number",
                                    "mol_weight_axis": "Apparent MW"
                                }
                            ),
                            ui.input_switch("combined_mw", "Mark combined MW", False),
                            ui.input_slider("peak_intensity_threshold", "Peak Intensity Threshold", 0, 150, 0, step=1, ticks=True),
                            ui.input_slider("gausian_smoothing", "Gaussian Smoothing", 0, 2, 0, step=0.1, ticks=True),
                        ),
                        ui.accordion_panel(
                            "Protein Search List",
                            ui.input_radio_buttons(
                                "protein_search_input_type", "Select Input Type",
                                {
                                    "uniprot_id": "Uniprot ID",
                                    "protein_name": "Protein Name",
                                    "gene_name": "Gene Name",
                                    "description": "Protein Description",
                                }
                            ),
                            ui.input_text_area("search_protein_list", "Search Proteins", placeholder="Enter Uniprot IDs, protein names or genes"),
                        ),
                        ui.accordion_panel(
                            "Sequence Coverage",
                            ui.input_file("file_upload_pr_matrix", "Upload pr matrix", accept=[".tsv"], multiple=False),
                            ui.input_file("file_upload_fasta", "Upload Proteome FASTA", accept=[".fasta"], multiple=False),
                            ui.input_slider("sequence_coverage_slice", "Slice Selection", 1, 192, 1, step=1, ticks=True),
                        ),
                        multiple=True, open=['Plot Mode', 'File Upload']
                    ),
                    open='always', width=300
                ),
                output_widget('plot_selected_proteins', width='100%', fill=False),
                ui.output_data_frame('prot_selection'),
                ui.output_plot('sequence_coverage', width='100%', height=200, fill=False),
            ),
        ),

        # Molecular Weight Calculation Tab
        ui.nav_panel(
            "Mol. Weight Calculation",
            ui.layout_sidebar(
                ui.sidebar(
                    ui.accordion(
                        ui.accordion_panel(
                            "Curve Fitting",
                            ui.input_slider("filter_intens_mol_weight", "Intensity Threshold", 50, 100, 98, step=1, ticks=True),
                            ui.input_slider("skip_slice", "Exclude first x Slices", 0, 50, 28, step=1, ticks=True),
                            ui.download_button("download_molweight", "Download Mol. Weight Data"),
                        ),
                    ),
                open='always', width=300
                ),
                ui.output_plot('plot_molweight', width=700, height=400, fill=False),
            ),
        ),

        # Multiple Condition Comparison Tab
        ui.nav_panel(
            "Multiple condition comparison (beta)",
            ui.layout_sidebar(
                ui.sidebar(
                    ui.accordion(
                        ui.accordion_panel(
                            "Upload",
                            ui.input_file("multiple_condition_upload", "Upload Files", multiple=True),
                            ui.input_file("file_upload_mol_mass_multiple", "Upload Proteome Mass File", accept=[".tsv"], multiple=False),
                        ),
                    ),
                open='always', width=300
                ),
                output_widget('plot_selected_proteins_multiple', width='100%', fill=False),
                ui.output_data_frame('display_dataframe_multiple'),
            ),
        ),

        selected='Welcome'
    ),
    theme=shinyswatch.theme.yeti
)

# Global storage for selected proteins
selected_proteins_global = set()



def server(input, output, session):
    """
    Define the server logic for the xiPROFILER application. This includes handling user inputs, 
    processing data, and rendering outputs such as plots and data frames.
    """

    @render.image
    def image():
        """
        Render the xiPROFILER logo as an image in the app.

        Returns:
            ImgData: The image data dictionary including the source path and dimensions.
        """
        from pathlib import Path

        dir = Path(__file__).resolve().parent
        img: ImgData = {"src": str(dir / "images/RapLabTextLogo.png"), "width": "200px"}
        return img

    @reactive.Calc
    def uploaded_pg_file():
        """
        Handle the uploaded pg-matrix file and read it as a pandas DataFrame.

        Returns:
            pd.DataFrame: The uploaded pg-matrix data or an empty DataFrame if no file is uploaded.
        """
        file: list[FileInfo] | None = input.file_upload_pg_matrix()
        if file is None:
            return pd.DataFrame()
        return pd.read_csv(file[0]["datapath"], sep='\t')

    @reactive.Calc
    def uploaded_file_mass():
        """
        Handle the uploaded proteome mass file or provide a default file based on organism selection.

        Returns:
            pd.DataFrame: The proteome mass data or an empty DataFrame if no file is uploaded.
        """
        file: list[FileInfo] | None = input.file_upload_mol_mass()
        if file is None:
            organism_files = {
                'e_coli': "ecoli_proteome.tsv",
                'human': "uniprotkb_taxonomy_id_9606_AND_reviewed_2024_10_21.tsv",
                'b_subtilis': "uniprotkb_taxonomy_id_224308_AND_review_2024_10_24.tsv",
                'mouse': "uniprotkb_taxonomy_id_10090_AND_reviewe_2024_11_04.tsv",
            }
            infile = Path(__file__).parent / organism_files.get(input.select_organism(), "")
            return pd.read_csv(infile, sep='\t') if infile.exists() else pd.DataFrame()

        return pd.read_csv(file[0]["datapath"], sep='\t')

    @reactive.Calc
    def uploaded_file_pr_matrix():
        """
        Handle the uploaded pr-matrix file and read it as a pandas DataFrame.

        Returns:
            pd.DataFrame: The uploaded pr-matrix data or an empty DataFrame if no file is uploaded.
        """
        file: list[FileInfo] | None = input.file_upload_pr_matrix()
        if file is None:
            return pd.DataFrame()
        return pd.read_csv(file[0]["datapath"], sep='\t')

    def read_protein_sequences(fasta_file: str) -> dict:
        """
        Parse a FASTA file to extract protein sequences.

        Args:
            fasta_file (str): Path to the FASTA file.

        Returns:
            dict: A dictionary mapping protein IDs to their sequences.
        """
        protein_sequences = {}
        for record in SeqIO.parse(fasta_file, "fasta"):
            protein_sequences[record.id.split('|')[2]] = str(record.seq)
        return protein_sequences

    def visualize_sequence_coverage(protein_name: str, protein_seq: str, detected_peptides: list,
                                     highlighted_peptides: list, intensity_column: str) -> plt.Figure:
        """
        Visualize the sequence coverage of a protein.

        Args:
            protein_name (str): Name of the protein.
            protein_seq (str): Amino acid sequence of the protein.
            detected_peptides (list): List of peptides detected for the protein.
            highlighted_peptides (list): List of peptides with measured intensities.
            intensity_column (str): Column name representing the slice intensity.

        Returns:
            plt.Figure: Matplotlib figure of the sequence coverage plot.
        """
        protein_length = len(protein_seq)
        coverage = [0] * protein_length

        for peptide in detected_peptides:
            start_index = protein_seq.find(peptide)
            if start_index != -1:
                for i in range(start_index, start_index + len(peptide)):
                    coverage[i] += 1

        fig, ax = plt.subplots(figsize=(15, 2))
        ax.barh(0, protein_length, left=0, height=0.6, color='gainsboro', edgecolor='none')

        for peptide in detected_peptides:
            start_index = protein_seq.find(peptide)
            if start_index != -1:
                color = 'orangered' if peptide in highlighted_peptides else 'blue'
                alpha = 1 if peptide in highlighted_peptides else 0.3
                ax.barh(0, len(peptide), left=start_index, height=0.6, color=color, alpha=alpha, edgecolor='none')

        ax.set_xlim([0, protein_length])
        ax.set_yticks([])
        ax.set_xlabel('Amino Acid Position')
        ax.set_title(f'Sequence Coverage for Protein {protein_name} in Slice {intensity_column}')

        general_patch = plt.Line2D([0], [0], color='blue', alpha=0.3, lw=4, label='Total Sequence Coverage')
        in_slice_patch = plt.Line2D([0], [0], color='orangered', lw=4, label=f'Sequence Coverage in Slice {intensity_column}')
        ax.legend(handles=[general_patch, in_slice_patch], loc='upper right', bbox_to_anchor=(1, 1))

        return fig

    def process_protein_coverage(fasta_file: str, pr_matrix_df: pd.DataFrame, protein_name: str,
                                  intensity_column: str) -> plt.Figure:
        """
        Main function to process protein coverage visualization.

        Args:
            fasta_file (str): Path to the FASTA file.
            pr_matrix_df (pd.DataFrame): pr-matrix data as a DataFrame.
            protein_name (str): Name of the protein.
            intensity_column (str): Column name representing the slice intensity.

        Returns:
            plt.Figure: Matplotlib figure of the sequence coverage plot.
        """
        protein_sequences = read_protein_sequences(fasta_file)

        detected_peptides = pr_matrix_df.loc[pr_matrix_df['Protein.Names'] == protein_name, 'Stripped.Sequence'].tolist()
        highlighted_peptides = pr_matrix_df.loc[pr_matrix_df['Protein.Names'] == protein_name].dropna(
            subset=[intensity_column])['Stripped.Sequence'].tolist()

        if protein_name in protein_sequences:
            protein_seq = protein_sequences[protein_name]
            return visualize_sequence_coverage(protein_name, protein_seq, detected_peptides, highlighted_peptides,
                                               intensity_column)
        else:
            print(f"Protein {protein_name} not found in the provided FASTA file.")

    @render.plot
    def sequence_coverage():
        """
        Render the sequence coverage plot based on user-uploaded pr-matrix and FASTA files.

        Returns:
            plt.Figure: The sequence coverage plot or a placeholder plot if files are missing.
        """
        pr_matrix_df = uploaded_file_pr_matrix()
        if pr_matrix_df is not None and not pr_matrix_df.empty:
            fasta_file = input.file_upload_fasta()[0]["datapath"]
            protein_name = get_protein_name()
            intensity_column = str(input.sequence_coverage_slice())

            return process_protein_coverage(fasta_file, pr_matrix_df, protein_name, intensity_column)
        else:
            fig, ax = plt.subplots()
            ax.set_title('Sequence Coverage (pr matrix and fasta file not uploaded)')
            return fig

    
    def get_protein_name():
        """
        Retrieve the name of the selected protein based on the input plot mode.

        Returns:
            str: The name of the selected protein or an error message if multiple proteins are selected.
        """
        df = pd.DataFrame()
        if input.plot_mode() in ["show_median_normalized", "show_median_normalized_bl"]:
            df = filtered_normalized_dataframe()
            indices = selected_proteins()

        if input.plot_mode() in ["show_raw", "show_raw_bl", "show_raw_log"]:
            if input.prot_selection_selected_rows() is not None:
                df = add_masses()
                df_filter = filtered_normalized_dataframe()
                selected_indices = selected_proteins()
                selected_protein_names = df_filter.iloc[selected_indices]['Protein.Names'].tolist()
                indices = df.index[df['Protein.Names'].isin(selected_protein_names)].tolist()

        if len(indices) == 1:
            protein_name = df['Protein.Names'].iloc[indices[0]]
        else:
            protein_name = "ERROR: Multiple Proteins selected"
        return protein_name


    def multiple_condition_uploaded_files():
        """
        Handle multiple uploaded files for condition comparison.

        Returns:
            dict: A dictionary where keys are file names and values are DataFrames of the file content.
        """
        files: list[FileInfo] | None = input.multiple_condition_upload()
        if files is None:
            return {}

        file_data = {}
        for file_info in files:
            df = pd.read_csv(file_info["datapath"], sep='\t')
            file_data[file_info["name"]] = df

        return file_data

    @reactive.Calc
    def uploaded_file_mass_multiple():
        """
        Handle the uploaded proteome mass file for multiple conditions.

        Returns:
            pd.DataFrame: The uploaded proteome mass data or an empty DataFrame if no file is uploaded.
        """
        file: list[FileInfo] | None = input.file_upload_mol_mass_multiple()
        if file is None:
            return pd.DataFrame()
        return pd.read_csv(file[0]["datapath"], sep='\t')

    def add_masses_multiple():
        """
        Add mass information to the uploaded files for multiple conditions.

        Returns:
            dict: A dictionary where keys are file names and values are DataFrames with added mass information.
        """
        files_data = multiple_condition_uploaded_files()
        mass_df = uploaded_file_mass_multiple()

        if files_data:
            if mass_df is not None and not mass_df.empty:
                if 'Entry Name' in mass_df.columns:
                    df_mass_renamed = mass_df.rename(columns={'Entry Name': 'Protein.Names'})
                elif 'Entry_Name' in mass_df.columns:  # Early testing fallback
                    df_mass_renamed = mass_df.rename(columns={'Entry_Name': 'Protein.Names'})

                for file_name, df in files_data.items():
                    merged_df = df.merge(df_mass_renamed[['Protein.Names', 'Mass']], on='Protein.Names', how='left')
                    files_data[file_name] = merged_df
        else:
            files_data = {}

        return files_data

    @render.data_frame
    def display_dataframe_multiple():
        """
        Display the DataFrame for multiple uploaded files with selected columns.

        Returns:
            DataGrid: A rendered DataGrid of the processed data.
        """
        file_data = add_masses_multiple()
        if not file_data:
            data = pd.DataFrame({
                'No ': [' ', ' ', '   '],
                'Files': [' ', ' ', '   '],
                'Uploaded': [' ', ' ', '   '],
            })
            return render.DataGrid(data)

        first_file_name = list(file_data.keys())[0]
        df = file_data[first_file_name]

        columns_to_include = ['Protein.Names', 'Genes', 'Mass']
        filtered_df = df[[column for column in columns_to_include if column in df.columns]]

        return render.DataGrid(filtered_df, width="100%", selection_mode="rows", filters=True)

    @reactive.Calc
    def selected_proteins_multiple():
        """
        Track the indices of selected proteins across multiple conditions.

        Returns:
            list: List of selected protein indices.
        """
        if input.display_dataframe_multiple_selected_rows() is not None:
            indices = list(input.display_dataframe_multiple_selected_rows())
            indices.sort()
        else:
            indices = []
        return indices

    @render_widget
    def plot_selected_proteins_multiple():
        """
        Generate a plot for selected proteins across multiple conditions.

        Returns:
            go.Figure: A Plotly figure showing intensity of selected proteins across slices.
        """
        files_data = add_masses_multiple()
        if not files_data:
            fig = go.Figure()
            fig.update_layout(title='No File Uploaded')
            return fig

        descriptive_headers = ['Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'First.Protein.Description',
                            'Mass', 'Proteotypic', 'Stripped.Sequence', 'Modified.Sequence', 'Precursor.Charge',
                            'Precursor.Id', 'Entry Name', 'Protein names', 'Gene Names', 'Subcellular location [CC]']

        fig = go.Figure()

        for file_name, df in files_data.items():
            plot_data = df.drop(columns=descriptive_headers, errors='ignore')
            plot_data.columns = range(1, len(plot_data.columns) + 1)
            indices = selected_proteins_multiple()

            for index in indices:
                if index not in plot_data.index:
                    continue
                row_data = plot_data.iloc[index]
                protein_name = df['Protein.Names'].iloc[index]
                text_list = [protein_name] * len(plot_data.columns)

                fig.add_trace(go.Scatter(
                    x=list(plot_data.columns),
                    y=row_data.values,
                    mode='lines+markers',
                    name=f"{protein_name} ({file_name})",
                    text=text_list,
                    hovertemplate='<b>%{text}</b><br>Slice: %{x}<br>Intensity: %{y:.2f}<extra></extra>'
                ))

        fig.update_layout(
            xaxis_title='Slices',
            yaxis_title='Intensity',
            title='Intensity of Selected Proteins across Slices',
            legend=dict(x=1, y=1, traceorder="normal"),
            template='plotly_white',
            height=700,
            hovermode='x unified',
            xaxis=dict(range=[-1, len(plot_data.columns)])
        )

        return fig




    @reactive.Calc
    def selected_proteins():
        """
        Reactive expression to store the indices of selected rows in the table.

        Returns:
            list: A sorted list of selected row indices.
        """
        if input.prot_selection_selected_rows() is not None:
            indices = list(input.prot_selection_selected_rows())
            indices.sort()
        else:
            indices = []
        return indices

    def logarithmic_decay(x: np.ndarray, a: float, b: float, c: float) -> np.ndarray:
        """
        Define a logarithmic decay function.

        Args:
            x (np.ndarray): Independent variable values.
            a (float): Parameter representing the offset.
            b (float): Parameter scaling the logarithmic term.
            c (float): Parameter shifting the logarithmic term.

        Returns:
            np.ndarray: Computed values of the logarithmic decay function.
        """
        return a - b * np.log(x + c)

    def apply_logarithmic_decay_model(x: np.ndarray, y: np.ndarray, refitting: bool = False) -> tuple:
        """
        Fit a logarithmic decay model to the data and optionally refit after removing outliers.

        Args:
            x (np.ndarray): Independent variable data.
            y (np.ndarray): Dependent variable data.
            refitting (bool): Whether to refit the model after removing outliers.

        Returns:
            tuple: Fitted x and y values, model parameters, and R-squared value.
        """
        initial_guess = [max(y), 1, 1]
        params, _ = curve_fit(logarithmic_decay, x, y, p0=initial_guess)

        x_fit = np.linspace(min(x), max(x), 100)
        y_fit = logarithmic_decay(x_fit, *params)

        y_pred = logarithmic_decay(x, *params)
        r2_initial = r2_score(y, y_pred)

        if refitting:
            residuals = y - y_pred
            std_dev = np.std(residuals)
            threshold = 2 * std_dev

            outliers = np.abs(residuals) > threshold

            x_no_outliers = x[~outliers]
            y_no_outliers = y[~outliers]

            params_no_outliers, _ = curve_fit(logarithmic_decay, x_no_outliers, y_no_outliers, p0=params)
            y_pred_no_outliers = logarithmic_decay(x_no_outliers, *params_no_outliers)
            r2_no_outliers = r2_score(y_no_outliers, y_pred_no_outliers)

            params = params_no_outliers
            r2_initial = r2_no_outliers
            x_fit = np.linspace(min(x_no_outliers), max(x_no_outliers), 100)
            y_fit = logarithmic_decay(x_fit, *params)

        return x_fit, y_fit, params, r2_initial

    def add_masses() -> pd.DataFrame:
        """
        Add protein mass information to the main DataFrame by merging with the uploaded mass file.

        Returns:
            pd.DataFrame: The merged DataFrame containing protein mass and interaction information.
        """
        df = uploaded_pg_file()
        mass_df = uploaded_file_mass()
        if df is not None and not df.empty:
            if mass_df is not None and not mass_df.empty:
                if 'Entry Name' in mass_df.columns:
                    df_mass_renamed = mass_df.rename(columns={'Entry Name': 'Protein.Names'})
                elif 'Entry_Name' in mass_df.columns:
                    df_mass_renamed = mass_df.rename(columns={'Entry_Name': 'Protein.Names'})

                merged_df = df.merge(df_mass_renamed[['Protein.Names', 'Mass', 'Interacts with']], on='Protein.Names', how='left')
                return merged_df

    def median_normalize_df() -> pd.DataFrame:
        """
        Perform median normalization on the intensity data within the DataFrame.

        Returns:
            pd.DataFrame: The median-normalized DataFrame.
        """
        df = add_masses()
        if df is not None and not df.empty:
            descriptive_headers = ['Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'First.Protein.Description',
                                'Mass', 'Proteotypic', 'Stripped.Sequence', 'Modified.Sequence', 'Precursor.Charge',
                                'Precursor.Id', 'Interacts with']
            data_to_normalize = df.drop(columns=[col for col in descriptive_headers if col in df.columns], errors='ignore')
            data_to_normalize = data_to_normalize.apply(pd.to_numeric, errors='coerce')

            normalized_df = data_to_normalize.div(data_to_normalize.median(axis=1), axis=0)

            for col in descriptive_headers:
                if col in df.columns:
                    normalized_df[col] = df[col]

            median_normalized_df = normalized_df[df.columns]
            return median_normalized_df

    def min_max_normalized_df() -> pd.DataFrame:
        """
        Perform min-max normalization on the intensity data within the DataFrame.

        Returns:
            pd.DataFrame: The min-max normalized DataFrame.
        """
        df = add_masses()
        if df is not None and not df.empty:
            descriptive_headers = ['Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'First.Protein.Description',
                                'Mass', 'Proteotypic', 'Stripped.Sequence', 'Modified.Sequence', 'Precursor.Charge',
                                'Precursor.Id', 'Interacts with']

            data_to_normalize = df.drop(columns=[col for col in descriptive_headers if col in df.columns], errors='ignore')
            data_to_normalize = data_to_normalize.apply(pd.to_numeric, errors='coerce')

            normalized_data = data_to_normalize.sub(data_to_normalize.min(axis=1), axis=0).div(
                (data_to_normalize.max(axis=1) - data_to_normalize.min(axis=1)).replace(0, 1), axis=0)

            for col in descriptive_headers:
                if col in df.columns:
                    normalized_data[col] = df[col]

            normalized_data = normalized_data[df.columns]
            return normalized_data

    def filtered_normalized_dataframe() -> pd.DataFrame:
        """
        Filter the normalized DataFrame based on a peak intensity threshold.

        Returns:
            pd.DataFrame: The filtered DataFrame with rows meeting the intensity criteria.
        """
        df = median_normalize_df()
        if df is not None and not df.empty:
            descriptive_headers = ['Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'First.Protein.Description',
                                'Mass', 'Proteotypic', 'Stripped.Sequence', 'Modified.Sequence', 'Precursor.Charge',
                                'Precursor.Id', 'Interacts with']

            temp_df = df.drop(columns=descriptive_headers, errors='ignore')
            mask = temp_df.apply(lambda x: x > input.peak_intensity_threshold()).any(axis=1)
            filtered_df = df[mask]
            filtered_df.reset_index(drop=True, inplace=True)
            return filtered_df

    def calc_molweight():
        """
        Calculate molecular weights for protein slices using median-normalized data.

        Returns:
            tuple: Contains fitted values, parameters, R-squared values, and masked data.
        """
        df = median_normalize_df()
        descriptive_headers = ['Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'First.Protein.Description',
                            'Mass', 'Proteotypic', 'Stripped.Sequence', 'Modified.Sequence', 'Precursor.Charge',
                            'Precursor.Id', 'Interacts with']

        filter_amount = input.filter_intens_mol_weight()

        filtered_df = df.drop(columns=descriptive_headers, errors='ignore')
        filtered_df.columns = range(1, len(filtered_df.columns) + 1)
        max_value = filtered_df.max().max()
        filtered_df = filtered_df.map(lambda x: x if x >= max_value * ((100 - filter_amount) / 100) else None)

        long_df = filtered_df.reset_index().melt(id_vars=['index'], value_vars=filtered_df.columns)
        long_df.rename(columns={'index': 'OriginalIndex', 'variable': 'Slice', 'value': 'FilteredValue'}, inplace=True)
        long_df = long_df.merge(df['Mass'], left_on='OriginalIndex', right_index=True)
        long_df.dropna(subset=['FilteredValue'], inplace=True)

        median_mass_per_slice = long_df.groupby('Slice')['Mass'].median()

        x_data, y_data = median_mass_per_slice.index.to_numpy(), median_mass_per_slice.to_numpy()
        valid_mask = np.isfinite(x_data) & np.isfinite(y_data)
        x_data_masked, y_data_masked = x_data[valid_mask], y_data[valid_mask]

        filtered_pairs = [(x, y) for x, y in zip(x_data_masked, y_data_masked) if x > input.skip_slice()]
        x_data_masked, y_data_masked = map(list, zip(*filtered_pairs)) if filtered_pairs else ([], [])

        refitting = False
        x_fit, y_fit, params_filtered, r2_filtered = apply_logarithmic_decay_model(x_data_masked, y_data_masked, refitting)
        return x_fit, y_fit, params_filtered, r2_filtered, x_data, y_data, x_data_masked, y_data_masked

    @render.plot
    def plot_molweight():
        """
        Render a plot of molecular weights across slices.

        Returns:
            plt.Figure: The rendered plot or a placeholder plot if no data is available.
        """
        df = median_normalize_df()
        if df is None or df.empty:
            fig, ax = plt.subplots()
            ax.set_title('No File Uploaded')
            return fig

        x_fit, y_fit, params_filtered, r2_filtered, x_data, y_data, x_data_masked, y_data_masked = calc_molweight()
        logarithmic_decay_function = f"y = {params_filtered[0]:.4f} - {params_filtered[1]:.4f} * log(x + {params_filtered[2]:.4f})"
        r2_filtered = f"R² = {r2_filtered:.4f}"

        fig, ax = plt.subplots()
        ax.scatter(x_data, y_data, label='Original Data', color='blue', alpha=0.5)
        ax.scatter(x_data_masked, y_data_masked, label='Filtered Data', color='red', alpha=0.5)

        ax.plot(x_fit, y_fit, color='red', label=f'Fitted Logarithmic Decay Curve \n{logarithmic_decay_function}\n {r2_filtered}')

        ax.set_xlabel('Slices')
        ax.set_ylabel('Mass [Da]')
        ax.set_title('Mass distribution of highest Intensity Proteins, shown per Slice')
        plt.legend()
        plt.tight_layout()
        return fig
    
    def update_selected_proteins(selected_indices):
        """
        Update the global storage with the newly selected protein indices.

        Args:
            selected_indices (list): List of selected protein indices to be added to the global storage.
        """
        global selected_proteins_global
        selected_proteins_global.update(selected_indices)

    @render.data_frame
    def prot_selection():
        """
        Render the filtered DataFrame for protein selection and include previously selected proteins.

        Returns:
            DataGrid: A DataGrid of the filtered proteins with selectable rows.
        """
        df = filtered_normalized_dataframe()

        if df is None or df.empty:
            data = pd.DataFrame({
                'No ': [' ', ' ', '   '],
                'File': [' ', ' ', '   '],
                'Uploaded': [' ', ' ', '   '],
            })
            return render.DataGrid(data)

        columns_to_include = ['Mass', 'Protein.Names', 'Protein.Ids', 'Interacts with']
        filtered_df = df[[column for column in columns_to_include if column in df.columns]]

        selected_rows = filtered_df.index[filtered_df['Protein.Names'].isin(selected_proteins_global)].tolist()

        return render.DataGrid(filtered_df, selection_mode="rows", filters=True, width='100%')

    @reactive.Effect
    def track_selected_proteins():
        """
        Track currently selected proteins in the DataGrid and update the global storage.
        """
        current_selection = input.prot_selection_selected_rows()
        if current_selection is not None:
            update_selected_proteins(current_selection)

    def split_proteins_from_inputlist(protein_string: str) -> List[str]:
        """
        Split a string of protein identifiers into a list based on common delimiters.

        Args:
            protein_string (str): String of protein identifiers.

        Returns:
            List[str]: A list of protein identifiers.
        """
        proteins = re.split(r'[,\t\n\s;:|]+', protein_string.strip())
        return [protein for protein in proteins if protein]

    def prot_names_from_ids(protein_ids: list, df_mass: pd.DataFrame) -> list:
        """
        Map a list of protein IDs to their corresponding names using a mass DataFrame.

        Args:
            protein_ids (list): List of protein IDs.
            df_mass (pd.DataFrame): DataFrame containing protein information.

        Returns:
            list: List of protein names corresponding to the provided IDs.
        """
        protein_names = []
        for protein_id in protein_ids:
            protein_name = df_mass[df_mass['Entry'] == protein_id]['Entry Name'].values
            if len(protein_name) > 0:
                protein_names.append(protein_name[0])
        return protein_names

    def prot_names_from_genes(protein_genes: list, df_mass: pd.DataFrame) -> list:
        """
        Map a list of gene names to their corresponding protein names using a mass DataFrame.

        Args:
            protein_genes (list): List of gene names.
            df_mass (pd.DataFrame): DataFrame containing protein information.

        Returns:
            list: List of protein names corresponding to the provided gene names.
        """
        protein_names = []
        for protein_gene in protein_genes:
            protein_name = df_mass[df_mass['Gene Names'] == protein_gene]['Entry Name'].values
            if len(protein_name) > 0:
                protein_names.append(protein_name[0])
        return protein_names

    def index_from_prot_name(protein_names: list, df: pd.DataFrame) -> int:
        """
        Retrieve indices of proteins in the DataFrame based on their names.

        Args:
            protein_names (list): List of protein names.
            df (pd.DataFrame): DataFrame containing protein information.

        Returns:
            int: Indices of proteins matching the provided names.
        """
        protein_indexes = []
        for protein_name in protein_names:
            prot_index = df.index[df['Protein.Names'] == protein_name].tolist()
            if len(prot_index) > 0:
                protein_indexes.append(prot_index[0])
        return protein_indexes

    def prot_names_from_description(protein_descriptions: list, df: pd.DataFrame) -> list:
        """
        Map a list of protein descriptions to their corresponding names in the DataFrame.

        Args:
            protein_descriptions (list): List of protein descriptions.
            df (pd.DataFrame): DataFrame containing protein information.

        Returns:
            list: List of protein names corresponding to the provided descriptions.
        """
        protein_names = []
        for protein_description in protein_descriptions:
            protein_name = df[df['First.Protein.Description'] == protein_description]['Protein.Names'].values
            if len(protein_name) > 0:
                protein_names.append(protein_name[0])
        return protein_names

    def split_proteins_from_inputlist_no_space(protein_string: str) -> List[str]:
        """
        Split a string of protein identifiers into a list, excluding spaces as delimiters.

        Args:
            protein_string (str): String of protein identifiers.

        Returns:
            List[str]: A list of protein identifiers.
        """
        proteins = re.split(r'[,\t\n;:|]+', protein_string.strip())
        return [protein for protein in proteins if protein]

    @render.download(filename="download_csv")
    def download_molweight():
        """
        Generate a downloadable CSV file of molecular weights calculated for each slice.

        Returns:
            str: Path to the generated CSV file.
        """
        input_file = uploaded_pg_file()
        _, _, params_filtered, _, _, _, _, _ = calc_molweight()

        slice_nr = len(input_file.columns) - 6
        slices = list(range(1, slice_nr + 1))
        masses = []
        for slice_value in slices:
            try:
                mass_value = params_filtered[0] - params_filtered[1] * np.log(slice_value + params_filtered[2])
                masses.append(mass_value.round(1))
            except Exception:
                masses.append(None)

        df = pd.DataFrame({'Slice': slices, 'Apparent MW [Da]': masses})

        filename = "molweight.csv"
        df.to_csv(filename, index=False)
        return filename


    def secondary_to_primary(secondary_value):
        """
        Convert a secondary axis value to the primary axis value based on calculated molecular weight parameters.

        Args:
            secondary_value (float): The value from the secondary axis (e.g., apparent molecular weight).

        Returns:
            float or None: The corresponding primary axis value (e.g., slice number) or None if an error occurs.
        """
        _, _, params_filtered, _, _, _, _, _ = calc_molweight()
        try:
            result = np.exp((secondary_value * 1000 - params_filtered[0]) / -params_filtered[1]) + params_filtered[2]
            return result
        except Exception:
            return None

    @render_widget
    def plot_selected_proteins():
        """
        Render an interactive Plotly graph of the selected proteins' intensity profiles across slices.

        Returns:
            go.Figure: Plotly figure displaying intensity profiles of selected proteins.
        """
        df = pd.DataFrame()
        mass_df = uploaded_file_mass()
        
        if input.plot_mode() in ["show_median_normalized", "show_median_normalized_bl"]:
            df = filtered_normalized_dataframe()
            indices = selected_proteins()

        if input.plot_mode() in ["show_raw", "show_raw_bl", "show_raw_log"]:
            if input.prot_selection_selected_rows() is not None:
                df = add_masses()
                df_filter = filtered_normalized_dataframe()
                selected_indices = selected_proteins()
                selected_protein_names = df_filter.iloc[selected_indices]['Protein.Names'].tolist()
                indices = df.index[df['Protein.Names'].isin(selected_protein_names)].tolist()

        if input.plot_mode() in ["show_min_max_normalized"]:
            if input.prot_selection_selected_rows() is not None:
                df = min_max_normalized_df()
                df_filter = filtered_normalized_dataframe()
                selected_indices = selected_proteins()
                selected_protein_names = df_filter.iloc[selected_indices]['Protein.Names'].tolist()
                indices = df.index[df['Protein.Names'].isin(selected_protein_names)].tolist()

        if df is None or df.empty:
            fig = go.Figure()
            fig.update_layout(title="No File Uploaded")
            return fig

        descriptive_headers = ['Protein.Group', 'Protein.Ids', 'Protein.Names', 
                               'Genes', 'First.Protein.Description', 
                               'Mass', 'Proteotypic', 'Stripped.Sequence', 
                               'Modified.Sequence', 'Precursor.Charge', 'Precursor.Id',
                               'Interacts with']

        sigma = input.gausian_smoothing()
        window_size = max(3, int(sigma * 7) | 1)

        plot_data = df.drop(columns=descriptive_headers, errors='ignore')
        plot_data = plot_data.apply(pd.to_numeric, errors='coerce')
        plot_data.columns = range(1, len(plot_data.columns) + 1)

        if input.plot_mode() in ["show_median_normalized_bl", "show_raw_bl"]:
            plot_data.fillna(0, inplace=True)

        fig = go.Figure()

        max_intensity = 0
        min_intensity = float('inf')

        interacting_proteins_indices = []
        if input.show_interactors():
            if len(indices) == 1:
                for index in indices:
                    protein_name = df['Protein.Names'].iloc[index]
                    entry_name = mass_df[mass_df['Entry Name'] == protein_name]

                    if not entry_name.empty:
                        interaction_string = entry_name['Interacts with'].values[0]
                        if isinstance(interaction_string, str):
                            interacting_proteins_entry_names_list = [protein.strip() for protein in interaction_string.split(';')]

                            for interacting_protein in interacting_proteins_entry_names_list:
                                interacting_protein_name = mass_df[mass_df['Entry'] == interacting_protein]['Entry Name'].values

                                if len(interacting_protein_name) > 0:
                                    interacting_protein_name = interacting_protein_name[0]
                                    prot_index = index_from_prot_name([interacting_protein_name], df)

                                    if prot_index:
                                        interacting_proteins_indices += prot_index

        search_protein_indices = []
        if input.search_protein_list() != '':
            protein_string = input.search_protein_list()
            if input.protein_search_input_type() == 'uniprot_id':
                prot_index = index_from_prot_name(prot_names_from_ids(split_proteins_from_inputlist(protein_string), mass_df), df)
            elif input.protein_search_input_type() == 'protein_name':
                prot_index = index_from_prot_name(split_proteins_from_inputlist(protein_string), df)
            elif input.protein_search_input_type() == 'gene_name':
                prot_index = index_from_prot_name(prot_names_from_genes(split_proteins_from_inputlist(protein_string), mass_df), df)
            elif input.protein_search_input_type() == 'description':
                prot_index = index_from_prot_name(prot_names_from_description(split_proteins_from_inputlist_no_space(protein_string), df), df)

            if prot_index:
                search_protein_indices += prot_index

        def plotting(indices, marker_symbol, marker_size, max_intensity, min_intensity):
            for index in indices:
                if index not in plot_data.index:
                    continue

                row_data = plot_data.iloc[index]
                original_indices = np.where(~np.isnan(row_data))[0]

                max_intensity = max(max_intensity, row_data.max())
                min_intensity = min(min_intensity, row_data[row_data > 0].min())

                if sigma > 0:
                    row_data_interpolated = row_data.interpolate(method='linear', limit_direction='both')
                    smoothed_data = gaussian_filter1d(row_data_interpolated, sigma=sigma)
                    fig.add_trace(go.Scatter(x=original_indices, y=smoothed_data[original_indices], mode='lines+markers', name=df['Protein.Names'].iloc[index]))
                else:
                    protein_name = df['Protein.Names'].iloc[index]
                    text_list = [protein_name] * len(plot_data.columns)

                    fig.add_trace(
                        go.Scatter(
                            x=plot_data.columns,
                            y=row_data.values,
                            mode='lines+markers',
                            marker=dict(size=marker_size, symbol=marker_symbol),
                            name=protein_name,
                            text=text_list,
                            hovertemplate='<b>%{text}</b><br>Slice: %{x}<br>Intensity: %{y:.2f}<extra></extra>'
                        )
                    )

            return max_intensity, min_intensity

        max_intensity, min_intensity = plotting(indices, 'circle', 10, max_intensity, min_intensity)
        max_intensity, min_intensity = plotting(search_protein_indices, 'triangle-up', 8, max_intensity, min_intensity)
        max_intensity, min_intensity = plotting(interacting_proteins_indices, 'x', 6, max_intensity, min_intensity)

        if input.plot_mode() in ["show_median_normalized", "show_median_normalized_bl"]:
            fig_title = "Median Normalized Intensity of Selected Proteins across Slices"
            y_axis_title = "Median Normalized Intensity"
        elif input.plot_mode() in ["show_min_max_normalized"]:
            fig_title = "Min-Max Normalized Intensity of Selected Proteins across Slices"
            y_axis_title = "Min-Max Normalized Intensity"
        else:
            fig_title = "Intensity of Selected Proteins across Slices"
            y_axis_title = "Intensity"

        fig.update_layout(
            xaxis_title="Slices",
            yaxis=dict(
                title=y_axis_title,
            ),
            title=fig_title,
            legend=dict(x=1, y=1, traceorder="normal"),
            template='plotly_white',
            height=500,
            hovermode='x unified',
            hoverdistance=1
        )

        if input.plot_mode() in ["show_raw_log"]:
            fig.update_yaxes(type="log", title="Log10(Intensity)")

        try:
            _, _, params_filtered, _, _, _, _, _ = calc_molweight()
        except:
            params_filtered = [0, 0, 0]

        if input.x_axis() == "mol_weight_axis":
            new_ticks = range(abs(int(params_filtered[2])), len(plot_data.columns), int((len(plot_data.columns) - abs(int(params_filtered[2]))) / 15))

            new_tick_labels = [round((int(params_filtered[0]) - params_filtered[1] * np.log(x + params_filtered[2])) / 1000, 1)
                            for x in new_ticks if x > params_filtered[2] and (params_filtered[0] - params_filtered[1] * np.log(x + params_filtered[2])) > 0]
            valid_ticks = [x for x in new_ticks if x > params_filtered[2] and (params_filtered[0] - params_filtered[1] * np.log(x + params_filtered[2])) > 0]

            fig.update_xaxes(
                tickvals=valid_ticks,
                ticktext=new_tick_labels,
                overlaying='x',
                side='top',
                title_text='Apparent Molecular Weight (kDa)'
            )
            fig.update_layout(margin=dict(t=150))

        if input.combined_mw():
            if 'Mass' in df.columns:
                combined_mass = df.loc[indices, 'Mass'].sum()
                combined_mass += df.loc[interacting_proteins_indices, 'Mass'].sum()

                x_val = np.exp((combined_mass - params_filtered[0]) / -params_filtered[1]) - params_filtered[2]

                if x_val > math.ceil(abs(params_filtered[2])) and (params_filtered[0] - params_filtered[1] * np.log(x_val + params_filtered[2])) > 0:
                    if input.plot_mode() in ["show_raw_log"]:
                        fig.add_shape(
                        type="line",
                        x0=x_val, x1=x_val,
                        y0=min_intensity, y1=max_intensity,
                        line=dict(
                            color="Red",
                            width=3,
                            dash="dash"
                        )
                        )
                    else:
                        fig.add_shape(
                        type="line",
                        x0=x_val, x1=x_val,
                        y0=0, y1=max_intensity,
                        line=dict(
                            color="Red",
                            width=3,
                            dash="dash"
                        )
                        )
        return fig


app = App(app_ui, server)

