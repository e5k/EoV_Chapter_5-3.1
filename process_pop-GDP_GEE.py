#%%
import geemap
import ee
ee.Initialize()
import pandas as pd 
import geopandas as gpd
import numpy as np
import seaborn as sns
sns.set_context('talk')
import matplotlib.pyplot as plt


#%% Read gpd and export them to shapefiles for ingestion in GEE
gdf_radius = gpd.read_file("Output/radius_footprint.gpkg", layer='radius').rename({'Significant':'sig'},axis=1)
gdf_countryA = gpd.read_file("Output/radius_footprint.gpkg", layer='countryA').rename({'Significant':'sig'},axis=1)
gdf_countryS = gpd.read_file("Output/radius_footprint.gpkg", layer='countryS').rename({'Significant':'sig'},axis=1)
gdf_countryA_diff = gpd.read_file("Output/radius_footprint.gpkg", layer='countryA_diff').rename({'Significant':'sig'},axis=1)
gdf_countryS_diff = gpd.read_file("Output/radius_footprint.gpkg", layer='countryS_diff').rename({'Significant':'sig'},axis=1)
gdf_country = gpd.read_file("Output/radius_footprint.gpkg", layer='country').rename({'Significant':'sig'},axis=1)

# Setup radii
rad = [10,30,100,300]

# Setup GDP
gdpT = np.arange(2000,2021,1)
gdp = ee.ImageCollection("users/sebastienbiasse/wang_gdp_h")
gdp_proj = gdp.first().projection().getInfo()
gpd_ppp_vis = {
    'min': 0,
    'max': 1e6,
    'palette': geemap.get_palette_colors(cmap_name='BuPu')
}

# Setup GHS
ghsT = np.arange(1975,2031,5)
ghs = ee.ImageCollection("JRC/GHSL/P2023A/GHS_POP")
ghs_proj = ghs.first().projection().getInfo()

ghs_vis = {
  'min': 0.0,
  'max': 20.0,
  'palette': ['060606', '337663', '337663', 'ffffff']
}

#%%
# Some cleanup required to remove empty geometries
def cleanGeom(dfTmp):
    toDrop = []
    # Remove parts with less than 5 points
    for i in range(len(dfTmp)):
        geometry = dfTmp.loc[i, "geometry"]
        if geometry.geom_type == 'Polygon':
            exterior = geometry.exterior
            if len(exterior.coords) < 5:
                # dfTmp.at[i, "geometry"] = None
                toDrop.append(i)
            else:
                interiors = []
                for interior in geometry.interiors:
                    if len(interior.coords) >= 5:
                        interiors.append(interior)
                dfTmp.at[i, "geometry"] = Polygon(exterior, interiors)
        elif geometry.geom_type == 'MultiPolygon':
            polygons = []
            for polygon in geometry.geoms:
                exterior = polygon.exterior
                if len(exterior.coords) >= 5:
                    interiors = []
                    for interior in polygon.interiors:
                        if len(interior.coords) >= 5:
                            interiors.append(interior)
                    polygons.append(Polygon(exterior, interiors))
                # else:
                #     toDrop.append(i)
            dfTmp.at[i, "geometry"] = MultiPolygon(polygons)
    
    return dfTmp.drop(toDrop).dropna()

def getStats(df, name, rungdp=True, runghs=True, country=False, rad=rad, gdpT=gdpT, ghsT=ghsT):
    # Setting up columns for batch export
    col = df.columns.tolist()
    col = col[:-1]+['sum']
    # Looping through radii is necessary to avoid memory errors on GEE
    for r in rad:
        print(r)
        # Added the dropna() to remove the empty geometries
        if country:
            tmp = df.copy().dropna().reset_index(drop=True)
        else:
            tmp = df.loc[df['Radius']==r].copy().dropna().reset_index(drop=True)
            tmp = cleanGeom(tmp).reset_index(drop=True)
        shp = geemap.gdf_to_ee(tmp)
        
        if rungdp:
            for y in gdpT:
                print(f'Processing {name} {r} {y}')
                img = gdp.filterDate(f'{y}-01-01', f'{y}-12-31').first().select('b1')
                # Remove negative values
                imgC = img.where(img.lte(0),0)
                # Sum
                df_sum = imgC.reduceRegions(collection=shp,reducer=ee.Reducer.sum(),
                                                crsTransform= gdp_proj['transform'])
                                                # crs= gdp_proj['crs'],
                t_sum = ee.batch.Export.table.toDrive(collection=df_sum,description=f'gdp_{name}_{r}_{y}_sum',folder='exposure', selectors=col, fileFormat='csv')
                t_sum.start()
        
        if runghs:
            for y in ghsT:
                print(f'Processing {name} {r} km year {y}')
                img = ghs.filterDate(f'{y}-01-01', f'{y}-12-31').first().select('population_count')
                # Remove negative values
                imgC = img.where(img.lte(0),0)
                # Sum
                df_sum = imgC.reduceRegions(collection=shp,reducer=ee.Reducer.sum(),
                                                crsTransform= ghs_proj['transform'])
                t_sum = ee.batch.Export.table.toDrive(collection=df_sum,description=f'ghs_{name}_{r}_{y}_sum',folder='exposure', selectors=col, fileFormat='csv')
                t_sum.start()
        
#%% Start the export process
# Get gpd
getStats(gdf_countryS, 'countryS', runghs=False)
getStats(gdf_countryA, 'countryA', runghs=False)
getStats(gdf_country, 'country', runghs=False, country=True, rad=[5])

# Get ghs
getStats(gdf_countryA, 'countryA', rungdp=False)
getStats(gdf_countryS, 'countryS', rungdp=False)
getStats(gdf_country, 'country', rungdp=False, country=True, rad=[5])
