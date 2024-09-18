#%%
import pandas as pd 
import geopandas as gpd
import numpy as np
import utm 
from pyproj import CRS

# import os 
# os.chdir('/Users/seb/Documents/WORK/Projects/Exposure/EoV-Chapter')
from antimeridian_splitter.split_polygon import split_polygon
from antimeridian_splitter.geopolygon_utils import OutputFormat

#%% Read and merge GVP and significant eruptions
volcPth = 'Data/GVP_Volcanoes.xlsx'
significantPth = 'Data/significant.tsv'
# Read GVP
volc = pd.read_excel(volcPth)
# Read significant eruptions
significant = pd.read_csv(significantPth, delimiter='\t')

# Filter volcanoes between 80S and 84N
volc = volc[(volc.Latitude>=-80) & (volc.Latitude<=84)]
# keep only relevant columns
volc = volc[['Volcano Number', 'Volcano Name', 'Country', 'Region', 'Subregion', 'Latitude', 'Longitude', 'Elevation (m)']]

# Add column to store if a given volcano has had a significant eruption
volc['Significant'] = volc['Volcano Name'].isin(significant.Name.unique())
volc = volc[volc.Region != 'Antarctica']

# Add UTM and EPSG info 
getUTM = lambda x: pd.Series(utm.from_latlon(x['Latitude'], x['Longitude']))
getEPSG = lambda x: pd.Series(CRS.from_dict({'proj': 'utm', 'zone': x['ZoneN'], 'south': True if x['Latitude'] < 0 else False}).to_authority()[1])
volc[['Easting', 'Northing', 'ZoneN', 'ZoneL']] = volc.apply(getUTM, axis=1)
volc[['EPSG']] = volc.apply(getEPSG, axis=1)
volc[['EPSG']] = volc[['EPSG']].astype(int)

# Convert to GDF to use with the exposure chapter 
gdf = gpd.GeoDataFrame(
    volc, geometry=gpd.points_from_xy(volc.Longitude, volc.Latitude), crs="EPSG:4326"
)
gdf = gdf.set_index('Volcano Number')


#%% ## Prepare the volcano db
# Add radius 
rad = [10,30,100,300]
gdf_radius = gpd.GeoDataFrame()
for v in gdf.index:
    tmp = gdf.loc[[v]]
    tmp = tmp.to_crs(epsg=gdf.loc[v,'EPSG'])
    for r in rad:
        tmp2 = tmp.copy()
        tmp2['geometry'] = tmp2.buffer(r*1000)
        tmp2['Radius'] = r 
        tmp2 = tmp2.to_crs(epsg=4326)
        gdf_radius = pd.concat([gdf_radius, tmp2])

#%% Correct polygons that cross the antimeridian --> https://github.com/guidorice/antimeridian_splitter
def splitPoly(poly):
    return MultiPolygon(split_polygon(poly, OutputFormat.GeometryCollection))

# Set multi-index to ensure unique values
gdf_radius = gdf_radius.reset_index().set_index(['Volcano Number', 'Radius'])
# List to store polygons that cross the antimeridian
rows_toDrop = []
# Empty gdf to store split polygons
gdf_radius_split = gpd.GeoDataFrame()

# Loop through all polygons
for row in range(0,gdf_radius.shape[0]):
    print(f'{row+1}/{gdf_radius.shape[0]}')
    # Isolate the current rows
    tmpRow = gdf_radius.iloc[row]
    # Get indices
    idx = gdf_radius.iloc[[row]].index.get_level_values(0).values[0]
    radius = gdf_radius.iloc[[row]].index.get_level_values(1).values[0]

    # Convert gdf polygon to json
    poly = json.loads(json.dumps(shapely.geometry.mapping(tmpRow.geometry)))
    # Split the polygon
    polySplit = splitPoly(poly)
    # Convert back to gpd
    polySplit = gpd.GeoDataFrame(index=[0], crs='epsg:4326', geometry=[polySplit])

    # If the polygon is split in 2
    if len(polySplit.explode(index_parts=True)) == 2:
        # Save the index of the row
        rows_toDrop.append(gdf_radius.iloc[[row]].index.values[0])
        # Append the new split polygon
        df_tmp = gpd.GeoDataFrame({'geometry': polySplit.geometry, 'Volcano Number': idx, 'Radius': radius}, index=[0]).set_index(['Volcano Number', 'Radius'])
        gdf_radius_split = pd.concat([gdf_radius_split,df_tmp])

#%% Replace each row containing a polygon crossing the antimeridian by its split counterpart
gdf_radius.loc[gdf_radius_split.index.tolist(), 'geometry'] = gdf_radius_split['geometry'].reset_index()
# Rename
gdf_radius = gdf_radius.rename({
    'Elevation (m)': 'Elevation',
    'Volcano Number': 'vID',
    'Volcano Name': 'vName',
    'Country': 'vCountry',
    'Region': 'vRegion',
    'Subregion': 'vSubregion'}, axis=1)
# Reset the index
gdf_radius = gdf_radius.reset_index().set_index('vName')

#%% Read Natural Earth countries
countries = gpd.read_file('Data/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp')
countries = countries[['NAME_LONG', 'WOE_ID', 'CONTINENT', 'REGION_UN', 'SUBREGION', 'geometry']]
# Adapt some column names
countries = countries.rename({'NAME_LONG': 'Country', 'WOE_ID': 'CountryID', 'CONTINENT': 'Continent', 'REGION_UN': 'Region', 'SUBREGION': 'Subregion'}, axis=1)
countries['Area'] = countries.area

# Dissolve polygons by country (i.e. one polygon per country)
gdf_country_A = gpd.GeoDataFrame()
gdf_country_S = gpd.GeoDataFrame()
# Symmetric difference
gdf_country_A_diff = gpd.GeoDataFrame()
gdf_country_S_diff = gpd.GeoDataFrame()

for r in rad:
    ## Start with the union of all radii
    maskA = gdf_radius[gdf_radius['Radius'] == r]
    maskS = gdf_radius[(gdf_radius['Radius'] == r) & (gdf_radius['Significant'])]

    # Take the overlay of radii by country
    maskA = gpd.overlay(maskA, countries).dissolve('Country').reset_index()
    maskS = gpd.overlay(maskS, countries).dissolve('Country').reset_index()
    
    # Get rid of all volcano information since we are now on country granularity
    maskA = maskA[['Country', 'Radius', 'CountryID', 'Continent', 'Region', 'Subregion', 'geometry']]
    maskS = maskS[['Country', 'Radius', 'CountryID', 'Continent', 'Region', 'Subregion', 'geometry']]
    
    # Add back to original df
    gdf_country_A = pd.concat([gdf_country_A, maskA])
    gdf_country_S = pd.concat([gdf_country_S, maskS])
    
    ## Now take the symmetric difference between countries and footprints
    def makeDiff(df):
    # Take the difference
        df = df.set_index('Country')
        df_diff = countries.overlay(df[['geometry']], how='symmetric_difference').dropna(subset='Country').set_index('Country')
        df_diff['Radius'] = r
        
        # Prepare datasets for the join
        df = df.reset_index().set_index(['Country','Radius', 'CountryID', 'Continent', 'Region', 'Subregion'])
        df_diff = df_diff.reset_index().set_index(['Country','Radius', 'CountryID', 'Continent', 'Region', 'Subregion'])
        
        # This adds back to _diff all the countries that disappear when taking the inverse of the circular buffer 
        # (i.e., those countries - mostly islands - that are entirely exposed by a given buffer)
        df_diff = df_diff.join(df, how='outer', rsuffix='_r').drop('geometry_r', axis=1).reset_index().set_index('Country')
        df_diff['Area'] = df_diff.area
        
        # Now keep only the countries that are affected by the radius 
        # (this is to attempt preventing memory errors on GEE)
        df_diff = df_diff.join(countries.set_index(['Country'])[['Area']], rsuffix='_r') 
        df_diff['dArea'] = df_diff['Area_r'] - df_diff['Area']
        df_diff = df_diff[df_diff['dArea']!=0].drop(['Area', 'Area_r', 'dArea'], axis=1)
        return df_diff.reset_index()
    
    maskA_diff = makeDiff(maskA)
    maskS_diff = makeDiff(maskS)
    
    # Add back to original df
    gdf_country_A_diff = pd.concat([gdf_country_A_diff, maskA_diff])
    gdf_country_S_diff = pd.concat([gdf_country_S_diff, maskS_diff])
    

#%% Save all
gdf_radius.to_file("Output/radius_footprint.gpkg", layer='radius', driver="GPKG")
gdf_country_A.to_file("Output/radius_footprint.gpkg", layer='countryA', driver="GPKG")
gdf_country_S.to_file("Output/radius_footprint.gpkg", layer='countryS', driver="GPKG")
gdf_country_A_diff.to_file("Output/radius_footprint.gpkg", layer='countryA_diff', driver="GPKG")
gdf_country_S_diff.to_file("Output/radius_footprint.gpkg", layer='countryS_diff', driver="GPKG")
countries.to_file("Output/radius_footprint.gpkg", layer='country', driver="GPKG")