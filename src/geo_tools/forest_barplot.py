### LIBRARIES ###

# System
from pathlib import Path
import argparse 
import pdb

# Dataframe
import pandas as pd

# Geospatial
import xarray as xr
import geopandas as gpd
from geo_tools.utils import define_class, clip_array, import_forest

# Plotting
import matplotlib.pyplot as plt

### DIRECTORY ### 

# Specify the desired output folder path for figures
output_figure_dir = Path.cwd() / ".." / "output_figures"
# output_figure_dir = "/mnt/PROVIDE/firzar/output_figures"

# Check if the output folder exists, and create if not
Path.mkdir(output_figure_dir, exist_ok=True, parents=True)

# Specify the desired output folder path for processed dataset
output_data_dir = Path.cwd() / ".." / "output_data"

# Check if the output folder exists, and create if not
Path.mkdir(output_data_dir, exist_ok=True, parents=True)

### PIPELINE ###

# Calculate area from forest cover map - I think this can be moved to utils so I can use it 
# in plot_iso_forest.py as well
def forest_cover_area(iso_code):

    # Clipping forest cover data to country border
    iso_forest_cover, iso_border = clip_array(
        iso_code,
        xarray_id="Forest_cover"
    )

    # Calculate area 
    forest_clean = iso_forest_cover.where((iso_forest_cover <= 1) & (iso_forest_cover > 0))

    convert_to_mha = 1e+6
    iso_forest_cover_area = (
        forest_clean
        .count(["latitude", "longitude"])
        .values
        .flatten()[0]
    ) / convert_to_mha

    return iso_forest_cover_area

# Get file description
def get_file_description(xarray_id):
    database = import_forest(xarray_id)
    description = database.long_name

    return description

# Add data labels
def add_value_labels(ax, spacing=5):
    """Add labels to the end of each bar in a bar chart.

    Arguments:
        ax (matplotlib.axes.Axes): The matplotlib object containing the axes
            of the plot to annotate.
        spacing (int): The distance between the labels and the bars.
    """

    # For each bar: Place a label
    for rect in ax.patches:
        # Get X and Y placement of label from rect.
        y_value = rect.get_height()
        x_value = rect.get_x() + rect.get_width() / 2

        # Number of points between bar and label. Change to your liking.
        space = spacing
        # Vertical alignment for positive values
        va = 'bottom'

        # If value of bar is negative: Place label below bar
        if y_value < 0:
            # Invert space to place label below
            space *= -1
            # Vertically align label at top
            va = 'top'

        # Use Y value as label and format number with one decimal place
        label = "{:.1f}".format(y_value)

        # Create annotation
        ax.annotate(
            label,                      # Use `label` as label
            (x_value, y_value),         # Place label at end of the bar
            xytext=(0, space),          # Vertically shift label by `space`
            textcoords="offset points", # Interpret `xytext` as offset in points
            ha='center',                # Horizontally center label
            va=va)                      # Vertically align label differently for
                                        # positive and negative values.

# Plot forest fraction classes
def plot_fraction_class(fraction_class_df, iso_code, xarray_id):
    
    fig, ax = plt.subplots(figsize = (8,5))

    ax.bar(fraction_class_df["Variable"], fraction_class_df["Value"])

    # Set xticks size
    ax.tick_params(axis='x', labelsize=8, rotation=345)

    # Add data label on each bar
    add_value_labels(ax)

    # Set y label
    x_label = get_file_description(xarray_id) # Need to strip Mg ha-1 in Forest_agb
    ax.set(
        ylabel = "Percentage of area over forest cover (%)",
        xlabel = f"Classes of {x_label}",
        title = f"{iso_code}\n{x_label}"
    )

    plt.tight_layout()

    return print("Data description:", x_label)

# Settings 
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
                    prog='ProgramName',
                    description='What the program does',
                    epilog='Text at the bottom of help')
    
    parser.add_argument('--iso_code', help='ISO3 code for the selected country', default="DEU")
    parser.add_argument('--xarray_id', help='NetCDF identification', default="Forest_cover_fraction")
    parser.add_argument('--start_step', help='Input start step', default=0)
    parser.add_argument('--end_step', help='Input end step', default=0.1)
    parser.add_argument('--delta_diff', help='Timestep between two classes', default=0.1)
    parser.add_argument('--max_val', help='Maximum value in the dataframe', default=1.0)
    parser.add_argument('--var_unit', help='Unit for the variable being mapped', default="adimensional")
    
    args = parser.parse_args()

    iso_code = args.iso_code
    xarray_id = args.xarray_id
    start_step = float(args.start_step)
    end_step = float(args.end_step)
    delta_diff = float(args.delta_diff)
    max_val = float(args.max_val)
    unit = args.var_unit

    # Forest cover area
    # pdb.set_trace()
    deu_cover_area = forest_cover_area(iso_code)

    # Classes of selected forest variable
    deu_classes = define_class(
        iso_code,
        xarray_id,
        deu_cover_area,
        start_step,
        end_step,
        delta_diff,
        max_val        
    ) # Should return a CSV

    # Plotting
    plot_fraction_class(deu_classes, iso_code, xarray_id, unit)
