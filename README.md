# forest_data_analysis
For Forest Navigator project

# To start
- Download and clone this repository
- Setup conda forge by running as below

```
conda config --add channels conda-forge
```

- Create project environment using the file `geo_environment.yml`

```
conda env create -f geo_environment.yml -n <give_new_name_if_needed>

conda activate <give_new_name_if_needed>
```

- Initialise the folder organisation by running `_init.ipynb`

# Clipping forest data layer to country/subregion boundaries
Requirements:
- Forest data layers in .netCDF format, covering a larger spatial extent than the area we want to crop.
- A shapefile containing multipolygon boundaries to use for cropping the forest data layers.

You need to import the module `utils` to your notebook, then use the function `clip_array`. Check `notebook/forest_cover.ipynb` for the example.

The function will return:
- `xarray.dataarray` for the forest data layer that is cropped into the specified boundaries/spatial extent
- `geopandas.GeoDataFrame` for the multipolygon boundaries that are used to crop the forest data layer

You can use these layers for plotting, analysis, and to be exported as a GeoTiff.

