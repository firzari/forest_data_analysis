from pathlib import Path
import xarray as xr
import geopandas as gpd
import pandas as pd

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

    xarray_input.close()

    return country_forest, border

# Calculate area for each histogram bin
def define_class(
        iso_code: str, 
        xarray_id: str,
        forest_layer_year: int,
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
        forest_layer_year (int): 
            The data acquisition year as shown in the attribute table of xarray_id
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

    if xarray_id == "Forest_fragment":
        country_forest_clean = country_forest.where(country_forest < 28)
        country_forest.close()

    else: 
        country_forest_clean = country_forest
        country_forest.close()

    while end <= max_val:
        # print(start, end)

        forest_data = country_forest_clean.where(
            (country_forest_clean > start) & (country_forest_clean <= end)
        )

        calc_area = (
            forest_data
            .count(["latitude", "longitude"]) 
            .values
            .flatten()[0]
        ) / 1e+6 # in Mha

        # print(calc_area)

        var_name = xarray_id.replace("_", " ")
        column_name = f"{var_name} >{round(start, 1)} to <= {round(end, 1)}"
        area_series = pd.Series(calc_area, name=column_name)

        # # Prepare var_uid for the classes
        # lower = xarray_id.lower()
        # var_uid = f"{lower}_{round(start, 1)}_{round(end, 1)}"

        area.append(area_series)

        start += delta_step
        end += delta_step
    
    # Convert to dataframe 
    forest_class_df = pd.DataFrame(area)

    # Rename column
    forest_class_df.rename(columns={0 : "Value"}, inplace=True)

    # Reset index and rename it into variable
    forest_class_df = forest_class_df.reset_index()
    forest_class_df.rename(columns={"index" : "Variable"}, inplace=True)

    # Add a column for iso code
    forest_class_df.insert(
        loc=0,
        column="Region",
        value=iso_code
    )

    # Add a column for var_uid
    var_uid = xarray_id.lower()
    forest_class_df.insert(
        loc=2,
        column="Var_uid",
        value=var_uid
    )

    # Add a column for year
    forest_class_df.insert(
        loc=3,
        column="Year",
        value=forest_layer_year
    )

    # Add a column for unit
    forest_class_df.insert(
        loc=4,
        column="Unit",
        value="Million hectares"
    )

    return forest_class_df
