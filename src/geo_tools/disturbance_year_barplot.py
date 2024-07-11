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
from geo_tools.utils import clip_array
from geo_tools.forest_barplot import get_file_description

# Plotting
import matplotlib.pyplot as plt

### PIPELINE ###

# Calculate total forest area affected by disturbances in every year
def yearly_disturbance(
        iso_code, 
        xarray_id,
        delta_diff=1,
        start_step=1986,
        max_val=2020
        ): 

    # Clipping forest cover fraction
    country_forest, iso_border = clip_array(
        iso_code,
        xarray_id
    )
    
    start = start_step
    area_mha = []
    delta_step = delta_diff

    while start <= max_val:
        # print(start, end)

        disturbance_year = country_forest.where(country_forest == start)

        calc_area = (
            disturbance_year
            .count(["latitude", "longitude"]) 
            .values
            .flatten()[0]
        ) / 1e+6 # in Mha

        # print(start, calc_area)

        column_name = f"{start}"
        area_series = pd.Series(calc_area, name=column_name)

        area_mha.append(area_series)

        start += delta_step
    
    # Convert to dataframe 
    dist_year_df = pd.DataFrame(area_mha)

    # Rename column
    dist_year_df.rename(columns={0 : "Value"}, inplace=True)

    # Reset index and rename it into variable
    dist_year_df = dist_year_df.reset_index()
    dist_year_df.rename(columns={"index" : "Year"}, inplace=True)

    # Add a column for iso code
    dist_year_df.insert(
        loc=0,
        column="Region",
        value=iso_code
    )

    # Add a column for variable
    dist_year_df.insert(
        loc=1,
        column="Variable",
        value="Disturbance Year"
    )

    # Add a column for var_uid 
    var_uid = xarray_id.lower()
    dist_year_df.insert(
        loc=2,
        column="Var_uid",
        value=var_uid
    )

    # Add a column for unit
    dist_year_df.insert(
        loc=4,
        column="Unit",
        value="Million hectares"
    )

    return dist_year_df

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
def plot_fraction_class(fraction_class_df, iso_code, xarray_id, unit):
    fig, ax = plt.subplots(figsize = (15,5))

    ax.bar(fraction_class_df.index, fraction_class_df["Perc. Forest Cover (%)"])

    # Set xticks size
    ax.tick_params(axis='x', labelsize=8, rotation=345)

    # Add data label on each bar
    add_value_labels(ax)

    # Set y label
    x_label = get_file_description(xarray_id)
    ax.set(
        ylabel = "Percentage of area over forest cover (%)",
        xlabel = f"Classes of {x_label} ({unit})",
        title = f"{iso_code}\n{x_label}"
    )

    plt.tight_layout()
    
    return print("Data description:", x_label)
