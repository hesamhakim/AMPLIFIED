#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script to generate a summary statistics table combining input OTU information
with filtered output results and then generate visualizations.

For each library and filtering mode, this script reports:
  - Number of samples in the input OTU file.
  - Number of species (rows) in the input OTU file.
  - Number of species that passed the filter (unique species among the filtered output).
  - Total number of filtered entries.
  - Count of species flagged as false positives.

Two visualizations are generated:
  1. A bar chart comparing input species vs. filtered species.
  2. A bar chart comparing input sample count vs. total filtered entries.

The summary CSV and the PNG files are saved in the specified output directory.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def generate_summary_stat(integrated_df, input_stats, output_path):
    """
    Generates a summary statistics table and visualizations using both filtered
    (output) data and input OTU file metrics.
    
    Parameters:
      integrated_df (pd.DataFrame): The concatenated filtered results DataFrame.
                                    Expected columns include: 'reference_lib', 'mode',
                                    'ID', and 'false_pos_prone'.
      input_stats (list of dict): A list of dictionaries. Each dictionary has keys:
                                    'library', 'mode', 'num_input_samples', 'num_input_species'
      output_path (str): File path where the summary CSV will be saved.
      
    Returns:
      summary_df (pd.DataFrame): The final summary statistics DataFrame.
    """
    
    # Summarize the filtered (output) data by library and mode.
    output_summary = (
        integrated_df.groupby(['reference_lib', 'mode'])
        .agg(
            filtered_species_count = ('ID', pd.Series.nunique),
            total_entries_output = ('ID', 'size'),
            false_positive_species_count = ('false_pos_prone', lambda x: (x == 'Yes').sum())
        )
        .reset_index()
    )
    
    # Convert input stats (list of dicts) to DataFrame.
    input_summary_df = pd.DataFrame(input_stats)
    
    # Merge input and output summaries. (If a given library/mode did not produce any
    # filtered data, the merge will generate NaN; these are replaced by 0.)
    summary_df = pd.merge(
        input_summary_df,
        output_summary,
        left_on=['library', 'mode'],
        right_on=['reference_lib', 'mode'],
        how='left'
    )
    
    # Clean up the merged table.
    summary_df.drop(columns=['reference_lib'], inplace=True)
    summary_df.rename(
        columns={
            'library': 'reference_lib',
            'num_input_samples': 'input_samples',
            'num_input_species': 'input_species'
        },
        inplace=True
    )
    summary_df.fillna(0, inplace=True)
    
    # Save the summary table.
    summary_df.to_csv(output_path, index=False)
    print(f"Summary statistics CSV written to: {output_path}")
    
    # Create a helper label for plotting.
    summary_df["label"] = summary_df["reference_lib"] + " (" + summary_df["mode"].astype(str) + ")"
    
    # -------------------------------------------------------------------------
    # Visualization 1: Bar chart for Input Species vs. Filtered Species Count
    # -------------------------------------------------------------------------
    plt.figure(figsize=(10,6))
    width = 0.35  # width of the bars
    
    x = np.arange(len(summary_df))
    plt.bar(x - width/2, summary_df["input_species"], width, label="Input Species", color='#4c72b0')
    plt.bar(x + width/2, summary_df["filtered_species_count"], width, label="Filtered Species", color='#55a868')
    
    plt.xticks(x, summary_df["label"], rotation=45, ha='right')
    plt.ylabel("Species Count")
    plt.title("Input vs. Filtered Species Count by Library and Mode")
    plt.legend()
    plt.tight_layout()
    
    viz1_path = os.path.join(os.path.dirname(output_path), "species_comparison.png")
    plt.savefig(viz1_path)
    plt.close()
    print(f"Visualization 1 saved to: {viz1_path}")
    
    # -------------------------------------------------------------------------
    # Visualization 2: Bar chart for Input Sample Count vs. Total Filtered Entries
    # -------------------------------------------------------------------------
    plt.figure(figsize=(10,6))
    plt.bar(x - width/2, summary_df["input_samples"], width, label="Input Samples", color='#c44e52')
    plt.bar(x + width/2, summary_df["total_entries_output"], width, label="Filtered Entries", color='#8172b3')
    
    plt.xticks(x, summary_df["label"], rotation=45, ha='right')
    plt.ylabel("Count")
    plt.title("Input Samples vs. Total Filtered Entries by Library and Mode")
    plt.legend()
    plt.tight_layout()
    
    viz2_path = os.path.join(os.path.dirname(output_path), "samples_vs_entries.png")
    plt.savefig(viz2_path)
    plt.close()
    print(f"Visualization 2 saved to: {viz2_path}")
    
    return summary_df

if __name__ == '__main__':
    # Allow running this script independently for testing.
    import sys
    if len(sys.argv) != 3:
        print("Usage: python summary_stat.py <integrated_csv_path> <output_summary_csv_path>")
        sys.exit(1)
    integrated_csv_path = sys.argv[1]
    output_summary_csv_path = sys.argv[2]
    integrated_df = pd.read_csv(integrated_csv_path)
    
    # In standalone mode, input metrics are not available so we pass an empty list.
    generate_summary_stat(integrated_df, [], output_summary_csv_path)