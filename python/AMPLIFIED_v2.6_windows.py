#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Refactored Script for Processing CLC Workbench OTU Tables

The final project directory is computed as:
    project_dir = permanent_dir + dynamic_dir
All other paths (data_dir, output_dir, config) are always treated as relative paths
to that project directory.
"""

import os
import sys
import argparse
import configparser
import pandas as pd

# =============================================================================
# Configuration Section – Update these variables as needed before running.
# =============================================================================

# Permanent (fixed) base directory – mapped network drive path.
PERMANENT_DIR = r"M:\Microbiology\2- Operation\1- Management Materials, Reports, Projects\Validation\Validation 2023\Targeted Metagenomic Sequencing\Assay Development"

# Dynamic subdirectory (update to the current Windows subfolder using backslashes)
DYNAMIC_DIR = r"BFTMS\Hesam Test"

# Batch name for the current run
BATCH_NAME = "batch1_20230720"

# Default subdirectory names (relative to project_dir)
DEFAULT_DATA_SUBDIR   = "data"
DEFAULT_OUTPUT_SUBDIR = r"output_tables_53"
DEFAULT_CONFIG_SUBDIR = os.path.join("Python", "config.parameters")

# =============================================================================
# Argument Parsing – All user-specified directory paths are treated as relative.
# =============================================================================
def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Process CLC Workbench OTU tables and generate output CSV files."
    )
    # These two determine the project directory
    parser.add_argument('--permanent_dir', default=PERMANENT_DIR,
                        help='Permanent base directory (unchanging part).')
    parser.add_argument('--dynamic_dir', default=DYNAMIC_DIR,
                        help='Dynamic subdirectory (changeable part).')
    # Other parameters (all user-specified paths will be interpreted relative to project_dir)
    parser.add_argument('--batch_name', default=BATCH_NAME,
                        help='Name of the batch (e.g., batch1_20230720).')
    parser.add_argument('--data_dir', default=None,
                        help='Path to the data directory (relative to the project directory).')
    parser.add_argument('--output_dir', default=None,
                        help='Directory to save output CSV files (relative to the project directory).')
    parser.add_argument('--config', default=None,
                        help='Path to the configuration file (relative to the project directory).')
    args = parser.parse_args()

    # Compute project directory
    project_dir = os.path.join(args.permanent_dir, args.dynamic_dir)

    # Always treat user-specified paths as relative to project_dir.
    if args.data_dir is None:
        args.data_dir = os.path.join(project_dir, DEFAULT_DATA_SUBDIR)
    else:
        args.data_dir = os.path.join(project_dir, args.data_dir)

    if args.output_dir is None:
        base_output_dir = os.path.join(project_dir, DEFAULT_OUTPUT_SUBDIR)
    else:
        base_output_dir = os.path.join(project_dir, args.output_dir)
    
    # Create batch-specific output directory
    args.output_dir = os.path.join(base_output_dir, args.batch_name)
    os.makedirs(args.output_dir, exist_ok=True)

    if args.config is None:
        args.config = os.path.join(project_dir, DEFAULT_CONFIG_SUBDIR)
    else:
        args.config = os.path.join(project_dir, args.config)

    return args

# =============================================================================
# Helper Functions
# =============================================================================
def load_config(config_path):
    config = configparser.ConfigParser()
    if not os.path.exists(config_path):
        sys.exit(f"Configuration file not found: {config_path}")
    config.read(config_path)
    return config

def mark_species(df, id_column, false_pos_tag, species_list):
    for species in species_list:
        df.loc[df[id_column].str.startswith(species), false_pos_tag] = 'Yes'
    return df

def normalize_and_correct(abund_tb, trim_tbl, samples, control, resid_val):
    for samp in samples:
        trimmed_reads = trim_tbl[trim_tbl["Name"] == samp]["Trimmed sequences"].values[0]
        abund_tb[f"{samp}_Normalized"] = abund_tb[samp] * trimmed_reads / 1_000_000
    for samp in samples:
        abund_tb[f"{samp}_Normalized_corrected"] = abund_tb[f"{samp}_Normalized"] - abund_tb[f"{control}_Normalized"]
        abund_tb[f"{samp}_S:N"] = (abund_tb[f"{samp}_Normalized_corrected"] + resid_val) / (abund_tb[f"{control}_Normalized"] + resid_val)
    return abund_tb

def apply_filtering(abund_tb, samples, control, criteria, standard_col_names, mode):
    select_IDs = pd.DataFrame(columns=standard_col_names)
    for sample in samples:
        select_columns = ['ID', 'Name', 'Taxonomy', 'reference_lib', 'batch',
                          'false_pos_prone', 'Combined', 'Min', 'Max', 'Mean',
                          'Median', 'Std', sample, control,
                          f"{sample}_Normalized", f"{control}_Normalized",
                          f"{sample}_Normalized_corrected", f"{sample}_S:N"]
        sample_abund = abund_tb[select_columns].copy()
        sample_abund['sample'] = sample
        sample_abund['control'] = control
        sample_abund.columns = standard_col_names
        sample_abund = sample_abund.sort_values(by='normalized_corrected', ascending=False).head(15)
        condition = (
            ((sample_abund['false_pos_prone'] == 'Yes') &
             (sample_abund['sample_abund'] >= float(criteria['yes_abund'])) &
             ((sample_abund['control_abund'] == 0) | (sample_abund['S:N'] >= float(criteria['sn_threshold']))) &
             (sample_abund['normalized_corrected'] >= float(criteria['norm_corrected'])))
            |
            ((sample_abund['false_pos_prone'] == 'No') &
             (sample_abund['sample_abund'] >= float(criteria['no_abund'])) &
             ((sample_abund['control_abund'] == 0) | (sample_abund['S:N'] >= float(criteria['sn_threshold']))) &
             (sample_abund['normalized_corrected'] >= float(criteria['norm_corrected'])))
        )
        filtered = sample_abund[condition]
        select_IDs = pd.concat([select_IDs, filtered], axis=0)
    select_IDs['mode'] = mode
    return select_IDs

def process_library(library, batch, data_dir, output_dir, config, integrated_data, input_summary_list):
    abund_tb_path = os.path.join(data_dir, batch, library, "OTU.xlsx")
    trim_tbl_path = os.path.join(data_dir, batch, "Trim_summary.xlsx")
    control_sample_tbl_path = os.path.join(data_dir, "control_sample_all_batches.xlsx")
    
    for path in [abund_tb_path, trim_tbl_path, control_sample_tbl_path]:
        if not os.path.exists(path):
            print(f"Required file not found: {path}. Skipping {library}.")
            return
    
    abund_tb = pd.read_excel(abund_tb_path)
    trim_tbl = pd.read_excel(trim_tbl_path)
    control_sample_tbl = pd.read_excel(control_sample_tbl_path)
    
    # Process columns for consistency
    abund_tb.columns = [col.replace("Abundance", "") for col in abund_tb.columns]
    abund_tb.columns = [col.replace("Species (Aggregated)", "Name") for col in abund_tb.columns]
    abund_tb.columns = [col.replace(" ", "") for col in abund_tb.columns]
    
    trim_tbl['Name'] = trim_tbl['Name'].str.replace(' (single)', '', regex=False)
    samples = trim_tbl["Name"].astype(str).tolist()
    
    try:
        control = control_sample_tbl.loc[control_sample_tbl['batch'] == batch, 'control_sample'].iloc[0]
    except IndexError:
        print(f"Control sample for batch '{batch}' not found. Skipping {library}.")
        return
    
    abund_tb['reference_lib'] = library
    abund_tb['batch'] = batch
    
    species_list = config[library]['false_pos_prone_species'].split(',')
    abund_tb['false_pos_prone'] = 'No'
    abund_tb = mark_species(abund_tb, 'ID', 'false_pos_prone', species_list)
    
    resid_val = float(config['common'].get('resid_val', '0.00'))
    abund_tb = normalize_and_correct(abund_tb, trim_tbl, samples, control, resid_val)
    
    standard_col_names = ['ID', 'Name', 'Taxonomy', 'reference_lib', 'batch',
                          'false_pos_prone', 'Combined', 'Min', 'Max', 'Mean',
                          'Median', 'Std', 'sample_abund', 'control_abund',
                          'sample_normalized', 'control_normalized',
                          'normalized_corrected', 'S:N', 'sample', 'control']
    
    filter_modes = config[library]['output_modes'].split(',')
    library_output_dir = os.path.join(output_dir, library)
    os.makedirs(library_output_dir, exist_ok=True)
    
    # Capture input statistics (per library)
    input_sample_count = len(samples)
    input_species_count = len(abund_tb)
    
    for mode in filter_modes:
        # Record input statistics for this library and mode.
        input_stats_record = {
            'library': library,
            'mode': mode,
            'num_input_samples': input_sample_count,
            'num_input_species': input_species_count
        }
        input_summary_list.append(input_stats_record)
        
        mode_section = f"{library}_{mode}"
        if mode_section not in config.sections():
            print(f"Warning: No configuration found for {mode_section}. Skipping.")
            continue
        criteria = {
            'yes_abund': config[mode_section]['yes_abund'],
            'no_abund': config[mode_section]['no_abund'],
            'norm_corrected': config[mode_section]['norm_corrected'],
            'sn_threshold': config[mode_section]['sn_threshold']
        }
        output_filename = f"{batch}_{library}-{mode}.csv"
        output_path = os.path.join(library_output_dir, output_filename)
        filtered_df = apply_filtering(abund_tb, samples, control, criteria, standard_col_names, mode)
        filtered_df.to_csv(output_path, index=False)
        integrated_data.append(filtered_df)
    
    print(f"Processing for {library} completed. Outputs are saved in {library_output_dir}.")


def main():
    args = parse_arguments()
    config = load_config(args.config)
    
    # Identify library types from sections that do not include an underscore.
    library_types = [section for section in config.sections() if '_' not in section]
    if not library_types:
        sys.exit("No library types found in the configuration file.")
    
    integrated_data = []
    input_summary_list = []  # New list to hold input OTU stats per library & mode.
    
    for library in library_types:
        print(f"Starting processing for library type: {library}")
        process_library(library, args.batch_name, args.data_dir, args.output_dir, config,
                        integrated_data, input_summary_list)
        print(f"Finished processing for library type: {library}\n")
    
    if integrated_data:
        integrated_df = pd.concat(integrated_data, ignore_index=True)
        # Save integrated filtered outputs in the batch output directory.
        integrated_csv_path = os.path.join(args.output_dir, "integrated_results.csv")
        integrated_df.to_csv(integrated_csv_path, index=False)
        print(f"Integrated CSV file created at: {integrated_csv_path}")
        
        # Call the updated summary_stat script to generate a summary table that combines input and output metrics.
        from summary_stat import generate_summary_stat
        summary_csv_path = os.path.join(args.output_dir, "summary_statistics.csv")
        generate_summary_stat(integrated_df, input_summary_list, summary_csv_path)
    else:
        print("No data to integrate.")

if __name__ == '__main__':
    main()