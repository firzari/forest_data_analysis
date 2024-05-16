from pathlib import Path
import xarray as xr
import geopandas as gpd
import pandas as pd

### DIRECTORY ### 

# # Specify the desired output folder path for processed dataset
# output_data_dir = Path.cwd() / ".." / "output_data"

# # Check if the output folder exists, and create if not
# Path.mkdir(output_data_dir, exist_ok=True, parents=True)

### FUNCTIONS ###

# Importing NetCDF file 
def import_forest(file_id: str) -> xr.DataArray:

    """Importing forest dataset whose naming convention starts with "Forest4model_v1_".
    Data is downloadable in https://gitlab.iiasa.ac.at/forestnavigator/wp2/public

    Args:
        file_id (str): File name of the forest data that can be identified in the words following
        the naming convention "Forest4model_v1_" i.e. "Canopy_height" is the words following
        "Forest4model_v1_Canopy_height.nc"

    Returns:
        xr.DataArray: Forest data file is imported as an xr.DataArray
    """

    input_dir = Path.cwd() / ".." / "input_data"

    # Directory to the NetCDF file
    file_name = Path(f"{input_dir}/Forest4model_v1_{file_id}.nc")

    # Read the data
    ds = xr.open_dataset(file_name)[file_id].isel(time=0)

    return ds

# Importing subregion borders
def import_boundary() -> gpd.GeoDataFrame:

    """Importing the default shapefile containing global country borders

    Returns:
        gpd.DataFrame: Global country borders shapefile is imported as a geopandas.DataFrame
    """
    input_dir = Path.cwd() / ".." / "input_data"

    subregion_borders = gpd.read_file(
        f"{input_dir}/WORLD_BORDERS/TM_WORLD_BORDERS-0.3.shp"
    ).set_index("ISO3")

    return subregion_borders

# Clipping forest data to country boundaries
def clip_array(
        iso_code: str, 
        xarray_id: str
    ):

    """Clipping xarray data of larger bounds to a smaller bounds, specifically this function
    clips xarray data of larger bounds/spatial extent to the spatial extent of a selected country

    Args:
        iso_code (str): The ISO3 code for the selected country i.e. "DEU"
        xarray_id (str): File name of the forest data that can be identified in the words following
        the naming convention "Forest4model_v1_" i.e. "Canopy_height" is the words following
        "Forest4model_v1_Canopy_height.nc"

    Returns:
        country_forest (xr.DataArray): An xarray data clipped into the spatial extent of a selected country
        border (gpd.DataFrame): A shapefile containing the polygon and the spatial extent of a selected country
    """
    # Import NetCDF
    xarray_input = import_forest(file_id=xarray_id)

    # Import country border
    shp_border = import_boundary()

    # Reprojection to the CRS in subregion borders
    ds_reproject = xarray_input.rio.write_crs(shp_border.crs)

    border = shp_border.loc[[iso_code]]
    country_forest = ds_reproject.rio.clip_box(
        *border.total_bounds
    ).rio.clip(border.geometry, crs=ds_reproject.rio.crs)

    return country_forest, border

# Calculate area for each histogram bin
def define_class(
        iso_code: str, 
        xarray_id: str,
        iso_forest_cover_area: xr.DataArray,
        start_step: float,
        end_step: float,
        delta_diff: float,
        max_val: float,
        ) -> pd.DataFrame: 
    
    """Dividing data values to specific and customisable bin classes. Can be used to generate a bar plot
    visualising the distribution of data values into bin classes and the frequency of each bin class

    Args:
        iso_code (str): 
            The ISO3 code for the selected country i.e. "DEU"
        xarray_id (str): 
            File name of the forest data that can be identified in the words following
            the naming convention "Forest4model_v1_" i.e. "Canopy_height" is the words following
            "Forest4model_v1_Canopy_height.nc"
        iso_forest_cover_area (xr.DataArray): 
            Total forest area calculated from "Forest4model_v1_Forest_cover.nc"
        start_step (float): 
            The starting value to initiate the creation of the first bin class, which is usually the minimum value
            of the dataset. 
            The value range of each bin is specified by delta_diff (description below).
            The next bin classes will be created automatically using the initiated start_step.
            For example, the minimum value is 1.0, then we can pass it into the start_step argument for 
            creating the first bin class. Then we set the value range to be included in each bin class, for example
            4.0, to be passed into delta_diff. Then it will create a first bin with the following range: 1.0 - 5.0 
            (delta_diff = 4.0)
        end_step (float): 
            The ending value for the the first bin class. If you'd like to have a delta_diff of 4.0, then the
            end_value could be set to 5.0 # Need to be modified, I think we need only start and delta_diff # 
        delta_diff (float): 
            The value range for each bin class
        max_val (float): 
            Maximum value from the dataset.

    Returns:
        pd.DataFrame: 
            A data frame containing the range of each bin class and its frequency, as well as ISO3 identification,
            total forest cover for comparison, and percentage area for each bin class over total forest cover
    """

    # Clipping forest cover fraction
    country_forest, iso_border = clip_array(
        iso_code,
        xarray_id
    )
    
    start = start_step
    end = end_step
    area = []
    delta_step = delta_diff

    while end <= max_val:
        # print(start, end)

        forest_data = country_forest.where((country_forest > start) & (country_forest <= end))

        calc_area = (
            forest_data
            .count(["latitude", "longitude"]) 
            .values
            .flatten()[0]
        ) / 1e+6 # in Mha

        # print(calc_area)

        column_name = f">{round(start, 1)} to <= {round(end, 1)}"
        area_series = pd.Series(calc_area, name=column_name)

        area.append(area_series)

        start += delta_step
        end += delta_step
    
    # Convert to dataframe 
    forest_class_df = pd.DataFrame(area)

    # Rename column
    forest_class_df.rename(columns={0 : "Count pixels (Mha)"}, inplace=True)

    # Add a column for forest cover area 
    forest_class_df["Forest cover (Mha)"] = round(iso_forest_cover_area, 2)

    # Add a column for iso code
    forest_class_df["ISO3"] = iso_code

    # Calculate percentage area of each fraction class over forest cover
    forest_class_df["Perc. Forest Cover (%)"] = (
        forest_class_df["Count pixels (Mha)"] / forest_class_df["Forest cover (Mha)"]
        ) * 100

    # # Export dataframe to excel
    # forest_class_df.to_csv(output_data_dir/f"{iso_code}_{xarray_id}_class.csv")

    return forest_class_df
