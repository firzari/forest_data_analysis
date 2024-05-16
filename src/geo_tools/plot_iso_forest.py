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

# Plotting
import matplotlib.pyplot as plt

### PIPELINE ###

# Plotting and export figure
def plot_forest(iso_code:str, xarray_id: xr.DataArray) -> plt:

    """Creating a distribution map of a relevant forest data

    Args:
        iso_code (str): 
            The ISO3 code for the selected country i.e. "DEU"
        xarray_id (xr.DataArray): 
            File name of the forest data that can be identified in the words following
            the naming convention "Forest4model_v1_" i.e. "Canopy_height" is the words following
            "Forest4model_v1_Canopy_height.nc"

    Returns:
        plt: A distribution map of a relevant forest data
    """

    # Clip array here
    country_forest, border = clip_array(
        iso_code,
        xarray_id
    )

    # Figure template
    figure, ax = plt.subplots()

    # Plot the distribution of the relevant forest data
    non_zero = country_forest.where(country_forest > 0)

    # Specially for aboveground biomass, we want to customise the variable label
    if xarray_id == "Forest_agb":
        ax_p = non_zero.plot.imshow(cmap='viridis', ax=ax, add_colorbar=False)
        cbar = figure.colorbar(ax_p)
        cbar.set_label('Forest aboveground biomass\n(AGB) 2020 in t ha-1\n[adimensional]')

    # For the rest of forest data, we keep the original variable name
    else:
        non_zero.plot.imshow(cmap='viridis', ax=ax)

    # Plot administration border
    border.boundary.plot(ax=ax, edgecolor="black", linewidth=0.5)

    # Additional information
    ax.axis("off")
    ax.set(title=f"{xarray_id}_{iso_code} in 2020")


# Calculate total area
def agg_total_area(
        iso_code: str, 
        xarray_id: xr.DataArray, 
        forest_layer_year: int
    ) -> pd.DataFrame:

    """Aggregating raster values to a selected spatial extent

    Args:
        iso_code (str): 
            The ISO3 code for the selected country i.e. "DEU"
        xarray_id (xr.DataArray): 
            File name of the forest data that can be identified in the words following
            the naming convention "Forest4model_v1_" i.e. "Canopy_height" is the words following
            "Forest4model_v1_Canopy_height.nc"
        forest_layer_year (int):
            The data acquisition year as shown in the attribute table of xarray_id.

    Returns:
        pd.DataFrame: A dataframe containing aggregated raster values to the spatial extent of a selected country 
    """
    # Subset EU data to country level
    country_forest, border = clip_array(
        iso_code,
        xarray_id
    )

    # Don't forget to select area where values > 0
    country_forest_clean = country_forest.where(country_forest > 0)

    convert_to_mha = 1e+6

    # Using the count function for canopy height because it contains non-binary values
    if xarray_id == "Canopy_height":
        total_area = (
            country_forest_clean
            .count(["latitude", "longitude"])
            .values
            .flatten()[0]
        ) / convert_to_mha
    
    # Using the sum function for the rest of forest data with binary values
    else:
        total_area = (
            country_forest_clean
            .sum(["latitude", "longitude"])
            .values
            .flatten()[0]
        ) / convert_to_mha

    # Remove underscore in xarray_id
    var_name = xarray_id.replace("_", " ")
    var_unit = "Million hectares"

    total_area_df = pd.DataFrame(
        {
            "Region" : [iso_code],
            "Variable" : [var_name],
            "Year" : [forest_layer_year],
            "Unit" : [var_unit],
            "Value" : [total_area]
        }
    )

    return total_area_df


# # Calculate total area and statistics
# # Need to think the data arch to include forest cover area in the stats description file
# def calc_stats(iso_code, xarray_id, var_unit):
    
#     # Clip array here
#     country_forest, border = clip_array(
#         iso_code,
#         xarray_id
#     )

#     # Don't forget to select area where values > 0
#     country_forest_clean = country_forest.where(country_forest > 0)

#     convert_to_mha = 1e+6

#     # Calculate total area 
#     if xarray_id == "Canopy_height":
#         total_area = (
#             country_forest_clean
#             .count(["latitude", "longitude"])
#             .values
#             .flatten()[0]
#         ) / convert_to_mha
    
#     else:
#         total_area = (
#             country_forest_clean
#             .sum(["latitude", "longitude"])
#             .values
#             .flatten()[0]
#         ) / convert_to_mha

#     # Calculate average 
#     average = (
#         country_forest_clean
#         .mean(["latitude", "longitude"])
#         .values
#         .flatten()[0]
#     )

#     # Standard deviation
#     std_dev = (
#         country_forest_clean
#         .std(["latitude", "longitude"])
#         .values
#         .flatten()[0]
#     )

#     # Minimum
#     minimum = (
#         country_forest_clean
#         .min(["latitude", "longitude"])
#         .values
#         .flatten()[0]
#     )

#     # Maximum
#     maximum = (
#         country_forest_clean
#         .max(["latitude", "longitude"])
#         .values
#         .flatten()[0]
#     )

#     # Range
#     range_vals = maximum - minimum

#     stats_df = pd.DataFrame(
#         {
#             "ISO3" : [iso_code],
#             "Total_area (Mha)" : [total_area],
#             f"Average ({var_unit})" : [average],
#             f"Standard deviation ({var_unit})" : [std_dev],
#             f"Maximum ({var_unit})" : [maximum],
#             f"Minimum ({var_unit})" : [minimum],
#             f"Value range ({var_unit})" : [range_vals]
#         }
#     )

#     return stats_df

