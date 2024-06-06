### TODO: STATS CALCULATION FROM PLOT_ISO_FOREST.PY TO BE MOVED HERE

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

# Calculate area for fragmentation change
def agg_fragment_change_area(
        iso_code: str, 
        xarray_id: xr.DataArray, 
        forest_layer_year: int,
        calculate_for="decrease"
    ) -> pd.DataFrame:

    # Subset EU data to country level
    country_forest, border = clip_array(
        iso_code,
        xarray_id
    )

    # Removing fill_value in Forest_fragment_change
    country_forest = country_forest.where(country_forest < 28)

    # Subset data: increase or decrease in fragmentation
    if calculate_for == "decrease":
        fragment_change = country_forest.where(country_forest < 0)
        var_name = "Forest fragment decrease"
        
    elif calculate_for == "increase":
        fragment_change = country_forest.where(country_forest > 0)
        var_name = "Forest fragment increase"

    convert_to_mha = 1e+6

    total_area = (
        fragment_change
        .sum(["latitude", "longitude"])
        .values
        .flatten()[0]
    ) / convert_to_mha

    # Remove underscore in xarray_id
    var_unit = "Million hectares"
    var_uid = var_name.replace(" ", "_").lower()

    total_area_df = pd.DataFrame(
        {
            "Region" : [iso_code],
            "Variable" : [var_name],
            "Var_uid" : [var_uid],
            "Year" : [forest_layer_year],
            "Unit" : [var_unit],
            "Value" : [total_area]
        }
    )

    # TODO: add a unique identifier for variable
    # uid = xarray_id; characters are in lower cases, spaces are replaced by "_"

    return total_area_df