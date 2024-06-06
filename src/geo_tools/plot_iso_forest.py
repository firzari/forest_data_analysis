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

# Plotting forest cover
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


# Plotting forest_type
def plot_forest_type(iso_code:str, xarray_id: xr.DataArray, forest_type_cat: int) -> plt:

    # Clip array here
    country_forest, border = clip_array(
        iso_code,
        xarray_id
    )

    # Figure template
    figure, ax = plt.subplots()

    # Subset to the selected forest type
    sel_forest_type = country_forest.where(country_forest == forest_type_cat)

    # Plotting forest type
    sel_forest_type.plot.imshow(ax=ax)

    # Plot administration border
    border.boundary.plot(ax=ax, edgecolor="black", linewidth=0.5)

    # Additional information
    ax.axis("off")
    ax.set(title=f"{xarray_id}_{iso_code} in 2020")


# Calculate total area for forest cover
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

    # Removing fill_value in Forest_fragment
    if xarray_id == "Forest_fragment":
        country_forest_clean = country_forest.where(country_forest < 28)
    else:
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

    # TODO: add a unique identifier for variable
    # uid = xarray_id; characters are in lower cases, spaces are replaced by "_"

    return total_area_df


# Calculate area for forest type
def agg_area_forest_type(iso_code: str, xarray_id: str, forest_type_cat: int, forest_layer_year: int):

    # Clip array here
    country_forest, border = clip_array(
        iso_code,
        xarray_id
    )

    # Subset to the selected forest type
    sel_forest_type = country_forest.where(country_forest == forest_type_cat)

    convert_to_mha = 1e+6

    total_area = (
        sel_forest_type
        .sum(["latitude", "longitude"])
        .values
        .flatten()[0]
    ) / convert_to_mha

    # Identify variable names
    if forest_type_cat == 1:
        var_name = xarray_id.replace("_", " ") + ": Undefined"
    elif forest_type_cat == 2:
        var_name = xarray_id.replace("_", " ") + ": Broadleaved"
    elif forest_type_cat == 3:
        var_name = xarray_id.replace("_", " ") + ": Coniferous"
    elif forest_type_cat == 4:
        var_name = xarray_id.replace("_", " ") + ": Mixed"

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


# Calculate average
def average_values(iso_code: str, xarray_id: str, forest_layer_year: int):

    # Clip array here
    country_forest, border = clip_array(
        iso_code,
        xarray_id
    )

    if xarray_id == "Canopy_height":
        avg_vals = (
            country_forest
            .mean(["latitude", "longitude"])
            .values
            .flatten()[0]
        )
        var_unit = "Meters"

    elif xarray_id == "Forest_agb":
        agb_above_zero = country_forest.where(country_forest > 0)
        avg_vals = (
            agb_above_zero
            .mean(["latitude", "longitude"])
            .values
            .flatten()[0]
        )
        var_unit = "t/ha"

    # Identify the variable name and the unit
    var_name = xarray_id.replace("_", " ")

    avg_vals_df = pd.DataFrame(
        {
            "Region" : [iso_code],
            "Variable" : [var_name],
            "Year" : [forest_layer_year],
            "Unit" : [var_unit],
            "Value" : [avg_vals]
        }
    )

    return avg_vals_df


def agg_area_natural_forest(iso_code: str, xarray_id: str, forest_layer_year: int, natural_forest:int = 2):

    # Clip array here
    country_forest, border = clip_array(
        iso_code,
        xarray_id
    )

    # Subset to the selected forest type
    sel_natural = country_forest.where(country_forest == natural_forest)

    convert_to_mha = 1e+6

    total_area = (
        sel_natural
        .count(["latitude", "longitude"])
        .values
        .flatten()[0]
    ) / convert_to_mha

    # Identify the variable name and the unit
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