# Emission-estimate

# Satellite data-driven algorithm for methane hotspot emission estimate
*Francesco Piroddu, researcher at Aarhus Universitet (Denmark)*

In this work a fast satellite data-driven algorithm for estimating hotspot methane emissions is presented and described step-by-step. At the scope of applying the developed methodology, a case study area located in the region of NW Canada, where an extensive fault zone, with both subparallel and criss-crossing faults, make up for several fields of faults of different cinematic behaviour and geometries has been delineated.
Briefly, it comprises the downloading of Sentinel-5p datastrips, the loading and creation of a faults database, the selection of source and background emission regions, the calling of the fundamental variables and the calculus of the regional emission estimate.
### 1. Fundamental libraries
The following list embraces the important and basic Python libraries utilized for developing the following algorithm. These libraries are ultimately imported and loaded and serve as the starting point.
```Python
import xarray as xr
import numpy as np
import pandas as pd 
import geopandas as gpd
import glob; import os
from sentinel5dl import search, download
from tqdm import tqdm
```
### 2. Data query and download from ESA Copernicus Hub
Satellite products query by area of interest vertices coordinates, date and time range, product type, processing level, processing mode and download of the products of interest.
```Python
result = search(
        polygon='POLYGON((-180.0 -90.0, 180.0 -90.0, 180.0 90.0, -180.0 90.0, -180.0 -90.0))',
        begin_ts='20XX-YY-01T00:00:00.000Z',
        end_ts='20XX-YY-01T23:59:59.999Z',
        product='L2__CH4___',
        processing_level='L2',
        processing_mode='Offline')

download(result.get('products'))
```
### 3. Faults database and emission regions
Loading of faults database (in this case the NW Canada fault zone) composed of hotspot emission regions selected by density, number and/or length of faults. Reading of the vertices coordinates of source regions of methane emissions and creation of the related background regions.
##### 3.1. Faults database import
```Python
hotsp = gpd.read_file('*/source_regions.geojson')
df = pd.DataFrame(hotsp)
geom = df['geometry']
δ_deg = 1.
```
##### 3.2. Source and background regions
```Python
geom_arr = []; sc_bound_lt = []; bg_bound_lt = []
for i in tqdm(range(len(geom))):
    sc_bound = np.array(geom[i].bounds, dtype=float)
    bg_bound = np.array([sc_bound[:2] - δ_deg, sc_bound[2:] + δ_deg], dtype=float)
    sc_bound_lt.append(sc_bound)
    bg_bound_lt.append(bg_bound)
```
### 4. Sentinel-5p database
Reading and loading of Sentinel-5p datastrips. Reading and calling of the following fundamental variables datasets: **Methane Mixing Ratio**, **Surface Pressure**, **Eastern Winds**, **Northern Winds** and **Surface Elevation**.
##### 4.1. Database reading
Initializing the variable list containers for collecting data variables from every datastrip file. Reading of the database of each NC file passed to reader as input, through file paths aiming at retrieving the fundamental variables and the related data.
```Python
ch4mm_lt = []; mixlon_lt = []; mixlat_lt = []
sfpres_lt = []; Ew_lt = []; Nw_lt = []; sfelv_lt = []
path = 'file_path'
for filename in tqdm(glob.glob(os.path.join(path, '*/*.nc'))):

    xr_data = xr.open_mfdataset(filename, group='PRODUCT',                                  
                                engine='netcdf4', data_vars='all', decode_coords=True)
    xr_data2 = xr.open_mfdataset(filename, group='PRODUCT/SUPPORT_DATA/INPUT_DATA',         
                                 engine='netcdf4', data_vars='all', decode_coords=True)                          
```
##### 4.2. Picking fundamental variables
Reading of Methane Mixing Ratio dataset.
```Python
    ch4mm = xr_data.methane_mixing_ratio
```
Reading the coordinates of Methane Mixing Ratio pixels to assign latitude and longitude of Methane Mixing Ratio as coordinate variables of the datasets of Surface Pressure, Eastern Winds, Northern Winds and Surface Elevation.
```Python
    latit = xr_data.methane_mixing_ratio.coords['latitude']
    long = xr_data.methane_mixing_ratio.coords['longitude']
```
**Methane mixing ratio dataset**: reading of coordinates and dataset calling.
```Python
    mixlon = np.array(ch4mm['longitude'], dtype=float)                         
    mixlat = np.array(ch4mm['latitude'], dtype=float)                          
    mixlon_lt.append(mixlon); mixlat_lt.append(mixlat)
    ch4mm = np.array(ch4mm, dtype=float)
    ch4mm_lt.append(ch4mm)
```
**Surface pressure dataset**: assignment of coordinates and dataset calling.
```Python
    sfpres = xr_data2.surface_pressure
    sfpres.coords == {}
    sfpres.coords['latitude'] = latit
    sfpres.coords['longitude'] = long
    sfpres = np.array(sfpres, dtype=float)
    sfpres_lt.append(sfpres)
```
**Eastern winds dataset**: assignment of coordinates and dataset calling.
```Python
    Ew = xr_data2.eastward_wind  
    Ew.coords == {}
    Ew.coords['latitude'] = latit
    Ew.coords['longitude'] = long
    Ew = np.array(Ew, dtype=float)
    Ew_lt.append(Ew)
```
**Northern winds dataset**: assignment of coordinates and dataset calling.
```Python
    Nw = xr_data2.northward_wind
    Nw.coords == {}
    Nw.coords['latitude'] = latit
    Nw.coords['longitude'] = long
    Nw = np.array(Nw, dtype=float)
    Nw_lt.append(Nw)
```
**Surface elevation dataset**: assignment of coordinates and dataset calling.
```Python
    sfelv = xr_data2.surface_altitude
    sfelv.coords == {}
    sfelv.coords['latitude'] = latit
    sfelv.coords['longitude'] = long
    sfelv = np.array(sfelv, dtype=float)
    sfelv_lt.append(sfelv)
```
### 5. Variables filtering
**Coordinate variables checking**: at this step we aim to establish which of the coordinate variables counts for the higher number of data, in order to check if the pixels of interest, belonging to the chosen emission source regions, have coupled values of coordinates, by selecting the coordinate variable with the lowest number of data in common.
In the next step, we select the source and background regions pixels by applying filtering conditions to the coordinate variables datasets, so as for picking the data contained in each source and background region. Then, out of this data only positive values and negative values in the latitude and longitude datasets are selected respectively.
Finally, the variables used to obtain the emission estimate each are just assigned their proper values belonging to the emission regions, by setting as a condition the filtered coordinate variable dataset chosen earlier after the checking step (latitude or longitude dataset).
##### 5.1. Checking coordinate dimensions
If statement which enables to choose between the two coordinate variables to work with for retrieving the data pixels of the emission regions and the fundamental variables.
```Python
if len(mixlon_lt) < len(mixlat_lt):
    print('input by longitude values')
elif len(mixlon_lt) > len(mixlat_lt):
    print('input by latitude values')
else:
    len(mixlon_lt) == len(mixlat_lt)
    print('input by either longitude or latitude is indifferent, all pixels of interest each have pixel centers with coordinates coupled!')
    print('working with longitude values...')
```
##### 5.2. Selection of source regions pixels
The following blocks serve as a filter of the dataset of the resulting coordinate dimension (latitude or longitude), filtered through the ranges set by the vertices coordinates of the incasing areas that contain the margins of the source regions surfaces, as can be seen in the attached figure 1.

- ***Latitude dataset filtering***
```Python
latlt_sc = []
el_sc = 0
while el_sc < len(sc_bound_lt):
    for i in range(len(mixlat_lt)):
        mixlat_lt[i] = np.resize(mixlat_lt[i], (1, 4172, 215))
        latvals_sc = np.where(mixlat_lt[i].any() > sc_bound_lt[el_sc][1].any() and \
                              mixlat_lt[i].any() < sc_bound_lt[el_sc][-1].any(), mixlat_lt[i], np.nan)
        lat_arr_sc = xr.DataArray(latvals_sc)
        lat_cond_sc = xr.where(lat_arr > 0, lat_arr_sc, 0)   
        latlt_sc.append(lat_cond_sc)
    el_sc += 1
```
- ***Longitude dataset filtering***
```Python
lonlt_sc = []
el_sc = 0; j = 0
while el_sc < len(sc_bound_lt):
    for j in range(len(mixlon_lt)):
        mixlon_lt[j] = np.resize(mixlon_lt[j], (1, 4172, 215))
        lonvals_sc = np.where(mixlon_lt[j].any() > sc_bound_lt[el_sc][0].any() and \
                              mixlon_lt[j].any() < sc_bound_lt[el_sc][2].any(), mixlon_lt[j], np.nan)
        lon_arr_sc = xr.DataArray(lonvals_sc)
        lon_cond = xr.where(lon_arr < 0, lon_arr_sc, 0)
        lonlt_sc.append(lon_cond_sc)
    el_sc += 1
```
- ***Backfiltering***

In this case, in order to target precisely the values belonging to the source surfaces ranges of methane emission, one of the two coordinate dimensions (latitude) is used to filter the other one (longitude, after par. 5.1.), thus selecting the very fundamental data pixels which account for every studied emission hotspot.
```Python
lgdat_sc_lt = []
el_sc = 0; k = 0
while el_sc < len(sc_bound_lt):
    for k in range(len(mixlon_lt)):
        for n in range(len(lonlt_sc)):
            lgdat_sc = np.where(mixlat_lt[k].any() > sc_bound_lt[el_sc][1].any() and \
                                mixlat_lt[k].any() < sc_bound_lt[el_sc][-1].any(), lonlt_sc[n], np.nan)
            lgdat_sc_lt.append(lgsc_dat)
    el_sc += 1
```
##### 5.3. Selection of background regions pixels

As stated in par. 5.2., similarly here the two dimension coordinates are filtered through the areas of background regions, setting the source regions as core areas and adding an incremental factor δ in degrees of coordinate values.

- ***Latitude dataset filtering***
```Python
latlt_bg = []
el_bg = 0
while el_bg < len(bg_bound_lt):
    for m in range(len(mixlat_lt)):
        mixlat_lt[m] = np.resize(mixlat_lt[m], (1, 4172, 215))
        latvals_bg = np.where(mixlat_lt[m].any() > bg_bound_lt[el_bg][0, 1].any() and \
                              mixlat_lt[m].any() < bg_bound_lt[el_bg][1, 1].any(), mixlat_lt[i], np.nan)
        lat_arr_bg = xr.DataArray(latvals_bg)
        lat_cond_bg = xr.where(lat_arr2 > 0, lat_arr_bg, 0)
        latlt_bg.append(lat_cond_bg)
    el_bg += 1
```
- ***Longitude dataset filtering***
```Python
lonlt_bg = []
el_bg = 0
while el_bg < len(bg_bound_lt):
    for n in range(len(mixlon_lt)):
        mixlon_lt[n] = np.resize(mixlon_lt[n], (1, 4172, 215))
        lonvals_bg = np.where(mixlon_lt[n].any() > bg_bound_lt[el_bg][0, 0].any() and \
                              mixlon_lt[n].any() < bg_bound_lt[el_bg][1, 0].any(), mixlon_lt[n], np.nan)
        lon_arr_bg = xr.DataArray(lonvals_bg)
        lon_cond_bg = xr.where(lon_arr2 < 0, lon_arr_bg, 0)
        lonlt_bg.append(lon_cond_bg)
    el_bg += 1
```
- ***Backfiltering***

As previously explained in par. 5.2., at this step the backfiltering of the longitude dimension via the latitude dataset is laid down, against the background surfaces ranges of values.
```Python
lgdat_bg_lt = []
el_bg = 0; k = 0
while el_sc < len(bg_bound_lt):
    for p in range(len(mixlon_lt)):
        for n in range(len(lonlt_bg)):
            lgdat_bg = np.where(mixlat_lt[p].any() > bg_bound_lt[el_bg][0, 0].any() and \
                                mixlat_lt[p].any() < bg_bound_lt[el_bg][1, 0].any(), lonlt_sc[n], np.nan)
            lgdat_bg_lt.append(lgbg_dat)
    el_sc += 1
```
##### 5.4. Variable values picking

The next two code blocks serve as assignment loops to assign the values belonging to the fundamental variables datasets that are indexed as in the coordinate dimension dataset (longitude or latitude). This dataset constitutes the baseline of the picked values of interest to be highlighted and drawn for laying down the layout of variables of any hotspot emission region, that will lead to evaluate its own emission estimate.

- ***Source regions values***

In the following statement the variables **Methane Mixing Ratio**, **Surface Pressure**, **Eastern Winds**, **Northern Winds** and **Surface Elevation** are assigned with their proper indexed values, by giving a source-region coordinate dimension as input, i.e. longitude values.
```Python
mixrt_sc_lt = []; supres_lt = []; Ewind_lt = []; Nwind_lt = []; suelev_lt = []
mixrt_sc_hsp = []; supres_hsp = []; Ewind_hsp = []; Nwind_hsp = []; suelev_hsp = []
el_sc = 0; j = 0
while el_sc < len(sc_bound_lt):
    for k in range(len(lgdat_sc_lt)):
        for sl in range(lgdat_sc_lt[k].shape[1]):
            for gp in range(lgdat_sc_lt[k].shape[2]):
                mixrt_sc = np.where(lgdat_sc_lt[k][0, sl, gp].values < 0, ch4mm[0, sl, gp], 0)
                supres = np.where(lgdat_sc_lt[k][0, sl, gp].values < 0, sfpres[0, sl, gp], 0)
                Ewind = np.where(lgdat_sc_lt[k][0, sl, gp].values < 0, Ew[0, sl, gp], 0)
                Nwind = np.where(lgdat_sc_lt[k][0, sl, gp].values < 0, Nw[0, sl, gp], 0)
                suelev = np.where(lgdat_sc_lt[k][0, sl, gp].values < 0, sfelv[0, sl, gp], 0)
                mixrt_sc_lt.append(mixrt_sc)
                supres_lt.append(supres)
                Ewind_lt.append(Ewind)
                Nwind_lt.append(Nwind)
                suelev_lt.append(suelev)
    mixrt_sc_hsp.append(mixrt_sc_lt)
    supres_hsp.append(supres_lt)
    Ewind_hsp.append(Ewind_lt)
    Nwind_hsp.append(Nwind_lt)
    suelev_hsp.append(suelev_lt)
    el_sc += 1
```
- ***Background regions values***

In the following statement the Methane Mixing Ratio is assigned with values from a background-region coordinate dimension dataset, i.e. longitude values.
```Python
mixrt_bg_lt = []; mixrt_bg_hsp = []
while el_bg < len(bg_bound_lt):
    for k in range(len(lgdat_bg_lt)):
        for sl in range(lgdat_bg_lt[k].shape[1]):
            for gp in range(lgdat_bg_lt[k].shape[2]):
                mixrt_bg = np.where(lgdat_bg_lt[k][0, sl, gp].values < 0, 0, ch4mm[0, sl, gp])
                mixrt_bg_lt.append(mixrt_bg)
    mixrt_bg_hsp.append(mixrt_bg_lt)
    el_bg += 1
```
### 6. Emission Estimate calculus

At this final step the process of calculating the **Emission Estimate** is defined. Starting from the definition and calculus of the **Methane Enhancement ΔXCH4**, then the parameters defining the **Correction Factor CF** are defined step-by-step. Finally, the Emission Estimate is calculated by the product of the Methane Enhancement and the Correction Factor.

##### 6.1. Methane Enhancement definition

Definition and calculus of Methane Enhancement by retrieving the average of the Methane Mixing Ratio values from the source regions dataset and the background regions dataset.

- ***Average for source regions Methane Mixing Ratio.***
```Python
mean_sc = []
for x in range(len(mixrt_sc_hsp)):
    avrg_sc = np.average(mixrt_sc_hsp[x])
    mean_sc.append(avrg_sc)
```
- ***Average for background regions Methane Mixing Ratio.***
```Python
mean_bg = []
for y in range(len(mixrt_bg_hsp)):
    avrg_bg = np.average(mixrt_bg_hsp[y])
    mean_bg.append(avrg_bg)
```
- ***Definition of Methane Enhancement.***
```Python
ΔXCH4 = []
for m in range(len(mean_sc)):
    for n in range(len(mean_bg)):
        met_enhanc = mean_sc[m] - mean_bg[n]                                
        ΔXCH4.append(met_enhanc)
```
##### 6.2. Correction Factor variables definition

Here the **Correction Factor CF** is defined step-by-step by the following parameters and variables. 
**Constant Conversion factor M**: this factor has a value of 5.345 x 10^-9 MtCH4 km^-12 ppb^-1 and is needed to convert a methane mole fraction variation into a methane mass variation per area for standard conditions, i.e. for surface pressure Psurf = 1013 hPa (Buchwitz et al., 2017).
```Python
M = 5.345e-09
```
Definition of **Surface Pressure** by averaging all grid cells data of the source regions, measured in hPa.
```Python
Pm = []
for lst2 in range(len(supres_hsp)):
    p_mean = np.average(np.array(supres_hsp[lst2]))/100                               
    Pm.append(p_mean)
```
The following code block defines the **Conversion Factor Mexp** by Surface Pressure as a dimensionless factor used to correct for the actual mass. It is obtained using the surface elevation map that is used for the determination of the *Elevation Correction* factor (Buchwitz et al., 2017). This correction factor is applied to reduce potential effects related to a location-dependent weighting of tropospheric and stratospheric contributions on XCH4 (Buchwitz et al., 2017).
```Python
Mexp_ps = []
for j in range(len(Pm)):
    Mexp = Pm[j]/1013.                                             
    Mexp_ps.append(Mexp)
```
In addition, here the **Conversion Factor Mexp** by Surface Elevation is defined as a dimensionless factor as well. H is the assumed scale height in km (8.5 km).
```Python
H = 8.5
Mexp_elv = []
for lst3 in range(len(suelev_hsp)):
    Mexp2 = np.average(np.exp(-(np.array(suelev_hsp[lst3]))/H))             
    Mexp_elv.append(Mexp2)
```
Definition of the transverse length of the determined source regions in kilometers through the conversion factor δ_km.
```Python
δ_km = 111.
el_sc = 0
dist = []
while el_sc < len(sc_bound_lt):
    L = -(sc_bound_lt[el_sc][0] - sc_bound_lt[el_sc][2])*δ_km                     
    dist.append(L)
    el_sc += 1
```
Definition of the effective wind velocity (**Northeastern Winds**) of air parcels travelling an effective **length L** over any source region along the NE direction, by summing up the two wind velocities in the N and E direction.
```Python
NE_wind = []
for lst1 in range(len(Ewind_hsp)):
    NEw = np.sqrt((np.array(Ewind_hsp[lst1]))**2 + (np.array(Nwind_hsp[lst1]))**2)
    NE_wind.append(NEw)
```
Statement for the dimensionless **factor C**.
```Python
C = 2.
```
##### 6.3. Conversion factor and Emission Estimate

The following blocks let us define and obtain values for the Correction Factor that, together with Methane Enhancement, leads the way for defining and calculating the Emission Estimate, both by Surface Pressure and Surface Elevation equivalently.

- ***Definition of Correction Factor by Surface Pressure***
```Python
CF_ps = []
for x in range(len(Mexp_ps)):
    for y in range(len(dist)):
        for z in range(len(NE_wind)):
            CF = M * Mexp_ps[x].any() * dist[y].any() * NE_wind[z].any() * C                    
            CF_ps.append(CF)
```
- ***Definition of Correction Factor by Surface Elevation***
```Python
CF_elv = []
for x in range(len(Mexp_elv)):
    for y in range(len(dist)):
        for z in range(len(Ne_wind)):
            CF = M * Mexp_elv[x].any() * dist[y].any() * NE_wind[z].any() * C                    
            CF_elv.append(CF)
```
- ***Definition of Emission Estimate by Surface Pressure***
```Python
Ee_ps = []
for a in range(len(ΔXCH4)):
    for b in range(len(CF_ps)):
        Ee = ΔXCH4[a].any() * CF_ps[b].any()
        Ee_ps.append(Ee)
```
- ***Definition of Emission Estimate by Surface Elevation***
```Python
Ee_elv = []
for a in range(len(ΔXCH4)):
    for b in range(len(CF_elv)):
        Ee = ΔXCH4[a].any() * CF_elv[b].any()
        Ee_elv.append(Ee)
```
##### 6.4. Emission estimate check

This final loop let us verify the equivalency of the two means for obtaining the Emission Estimate and check if the values obtained are equal or not.
```Python
for c in range(len(Ee_ps)):
    for d in range(len(Ee_elv)):
        if Ee_ps == Ee_elv:
            print('Emission estimate is the same using a mass correction factor calculated with either surface pressure or surface elevation')
        else:
            print('Emission estimate changes when using surface pressure or surface elevation for the mass correction factor')  
```
##### Acknowledgments

This product has been conceived and made in cooperation with Prof. Christoffer Karoff from Geoscience Department at Aarhus University (Denmark). In addition, funds had been enlarged by European Union via the Erasmus+ KA1 mobility programme for Traineeship, by the Italian Ministry of Education, University and Research and by Regione Autonoma della Sardegna.

##### References
Buchwitz, Michael, et al. "Satellite-derived methane hotspot emission estimates using a fast data-driven method." Atmospheric Chemistry and Physics 17.9 (2017): 5751-5774.
Government of Canada. CanadianGIS.com, Geographic Information and Geospatial Resources. © Copyright - Canadian GIS & Geomatics
https://canadiangis.com/data.php#British-Columbia

![faults_db2](https://user-images.githubusercontent.com/106487184/180883675-eb4af2ab-f44d-4345-ae6e-c78b634de0dd.png)
*Fig. 1. NW Canadian region and the related fault zone divided by 13 hotspot source emission areas as case studies for the methane emissions (courtesy of: Government of Canada).*
