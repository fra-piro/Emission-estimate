#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import xarray as xr
import numpy as np
import pandas as pd 
import geopandas as gpd
import glob; import os
from sentinel5dl import search, download
from tqdm import tqdm

result = search(
        polygon='POLYGON((-180.0 -90.0, 180.0 -90.0, 180.0 90.0, -180.0 90.0, -180.0 -90.0))',
        begin_ts='20XX-YY-01T00:00:00.000Z',
        end_ts='20XX-YY-01T23:59:59.999Z',
        product='L2__CH4___',
        processing_level='L2',
        processing_mode='Offline')

download(result.get('products'))

hotsp = gpd.read_file('*/source_regions.geojson')
df = pd.DataFrame(hotsp)
geom = df['geometry']
δ_deg = 1.

geom_arr = []; sc_bound_lt = []; bg_bound_lt = []
for i in tqdm(range(len(geom))):
    sc_bound = np.array(geom[i].bounds, dtype=float)
    bg_bound = np.array([sc_bound[:2] - δ_deg, sc_bound[2:] + δ_deg], dtype=float)
    sc_bound_lt.append(sc_bound)
    bg_bound_lt.append(bg_bound)
    

ch4mm_lt = []; mixlon_lt = []; mixlat_lt = []
sfpres_lt = []; Ew_lt = []; Nw_lt = []; sfelv_lt = []
path = 'file_path'
for filename in tqdm(glob.glob(os.path.join(path, '*/*.nc'))):

    xr_data = xr.open_mfdataset(filename, group='PRODUCT',                                  
                                engine='netcdf4', data_vars='all', decode_coords=True)
    xr_data2 = xr.open_mfdataset(filename, group='PRODUCT/SUPPORT_DATA/INPUT_DATA',         
                                 engine='netcdf4', data_vars='all', decode_coords=True)                          
    
    ch4mm = xr_data.methane_mixing_ratio
    
    latit = xr_data.methane_mixing_ratio.coords['latitude']
    long = xr_data.methane_mixing_ratio.coords['longitude']
    
    mixlon = np.array(ch4mm['longitude'], dtype=float)                         
    mixlat = np.array(ch4mm['latitude'], dtype=float)                          
    mixlon_lt.append(mixlon); mixlat_lt.append(mixlat)
    ch4mm = np.array(ch4mm, dtype=float)
    ch4mm_lt.append(ch4mm)
    
    sfpres = xr_data2.surface_pressure
    sfpres.coords == {}
    sfpres.coords['latitude'] = latit
    sfpres.coords['longitude'] = long
    sfpres = np.array(sfpres, dtype=float)
    sfpres_lt.append(sfpres)
    
    Ew = xr_data2.eastward_wind  
    Ew.coords == {}
    Ew.coords['latitude'] = latit
    Ew.coords['longitude'] = long
    Ew = np.array(Ew, dtype=float)
    Ew_lt.append(Ew)
    
    Nw = xr_data2.northward_wind
    Nw.coords == {}
    Nw.coords['latitude'] = latit
    Nw.coords['longitude'] = long
    Nw = np.array(Nw, dtype=float)
    Nw_lt.append(Nw)
    
    sfelv = xr_data2.surface_altitude
    sfelv.coords == {}
    sfelv.coords['latitude'] = latit
    sfelv.coords['longitude'] = long
    sfelv = np.array(sfelv, dtype=float)
    sfelv_lt.append(sfelv)
    

if len(mixlon_lt) < len(mixlat_lt):
    print('input by longitude values')
elif len(mixlon_lt) > len(mixlat_lt):
    print('input by latitude values')
else:
    len(mixlon_lt) == len(mixlat_lt)
    print('input by either longitude or latitude is indifferent, all pixels of interest each have pixel centers with coordinates coupled!')
    print('working with longitude values...')
    

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

lgdat_sc_lt = []
el_sc = 0; k = 0
while el_sc < len(sc_bound_lt):
    for k in range(len(mixlon_lt)):
        for n in range(len(lonlt_sc)):
            lgdat_sc = np.where(mixlat_lt[k].any() > sc_bound_lt[el_sc][1].any() and \
                                mixlat_lt[k].any() < sc_bound_lt[el_sc][-1].any(), lonlt_sc[n], np.nan)
            lgdat_sc_lt.append(lgsc_dat)
    el_sc += 1

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

lgdat_bg_lt = []
el_bg = 0; k = 0
while el_sc < len(bg_bound_lt):
    for p in range(len(mixlon_lt)):
        for n in range(len(lonlt_bg)):
            lgdat_bg = np.where(mixlat_lt[p].any() > bg_bound_lt[el_bg][0, 0].any() and \
                                mixlat_lt[p].any() < bg_bound_lt[el_bg][1, 0].any(), lonlt_sc[n], np.nan)
            lgdat_bg_lt.append(lgbg_dat)
    el_sc += 1

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

mixrt_bg_lt = []; mixrt_bg_hsp = []
while el_bg < len(bg_bound_lt):
    for k in range(len(lgdat_bg_lt)):
        for sl in range(lgdat_bg_lt[k].shape[1]):
            for gp in range(lgdat_bg_lt[k].shape[2]):
                mixrt_bg = np.where(lgdat_bg_lt[k][0, sl, gp].values < 0, 0, ch4mm[0, sl, gp])
                mixrt_bg_lt.append(mixrt_bg)
    mixrt_bg_hsp.append(mixrt_bg_lt)
    el_bg += 1

mean_sc = []
for x in range(len(mixrt_sc_hsp)):
    avrg_sc = np.average(mixrt_sc_hsp[x])
    mean_sc.append(avrg_sc)

mean_bg = []
for y in range(len(mixrt_bg_hsp)):
    avrg_bg = np.average(mixrt_bg_hsp[y])
    mean_bg.append(avrg_bg)

ΔXCH4 = []
for m in range(len(mean_sc)):
    for n in range(len(mean_bg)):
        met_enhanc = mean_sc[m] - mean_bg[n]                                
        ΔXCH4.append(met_enhanc)
        

M = 5.345e-09

Pm = []
for lst2 in range(len(supres_hsp)):
    p_mean = np.average(np.array(supres_hsp[lst2]))/100                               
    Pm.append(p_mean)

Mexp_ps = []
for j in range(len(Pm)):
    Mexp = Pm[j]/1013.                                             
    Mexp_ps.append(Mexp)

H = 8.5
Mexp_elv = []
for lst3 in range(len(suelev_hsp)):
    Mexp2 = np.average(np.exp(-(np.array(suelev_hsp[lst3]))/H))             
    Mexp_elv.append(Mexp2)

δ_km = 111.
el_sc = 0
dist = []
while el_sc < len(sc_bound_lt):
    L = -(sc_bound_lt[el_sc][0] - sc_bound_lt[el_sc][2])*δ_km                     
    dist.append(L)
    el_sc += 1

NE_wind = []
for lst1 in range(len(Ewind_hsp)):
    NEw = np.sqrt((np.array(Ewind_hsp[lst1]))**2 + (np.array(Nwind_hsp[lst1]))**2)
    NE_wind.append(NEw)

C = 2.

CF_ps = []
for x in range(len(Mexp_ps)):
    for y in range(len(dist)):
        for z in range(len(NE_wind)):
            CF = M * Mexp_ps[x].any() * dist[y].any() * NE_wind[z].any() * C                    
            CF_ps.append(CF)
            

CF_elv = []
for x in range(len(Mexp_elv)):
    for y in range(len(dist)):
        for z in range(len(Ne_wind)):
            CF = M * Mexp_elv[x].any() * dist[y].any() * NE_wind[z].any() * C                    
            CF_elv.append(CF)
            

Ee_ps = []
for a in range(len(ΔXCH4)):
    for b in range(len(CF_ps)):
        Ee = ΔXCH4[a].any() * CF_ps[b].any()
        Ee_ps.append(Ee)
        

Ee_elv = []
for a in range(len(ΔXCH4)):
    for b in range(len(CF_elv)):
        Ee = ΔXCH4[a].any() * CF_elv[b].any()
        Ee_elv.append(Ee)
        

for c in range(len(Ee_ps)):
    for d in range(len(Ee_elv)):
        if Ee_ps == Ee_elv:
            print('Emission estimate is the same using a mass correction factor calculated with either surface pressure or surface elevation')
        else:
            print('Emission estimate changes when using surface pressure or surface elevation for the mass correction factor')  
            

