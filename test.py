import sys
import numbers
import pybufrkit
import numpy as np 
import matplotlib.pyplot as plt 

import netCDF4

import gridpp

def get(filename, variable='ta'):
    """ Retrieve observations from a bufr file. Assumes all values are for the same time stamp

    Args:
        filename (str): Filename to parse
    Returns:
        ids (list): List of station IDs
        lats (list): List of latitudes
        lons (list): List of longitudes
        values (list): List of measurements
    """
    ids = list()
    lats = list()
    lons = list()
    values = list()
    lookup_names = ["TEMPERATURE/AIR TEMPERATURE", "TEMPERATURE/DRY-BULB TEMPERATURE"]
    required_names = [
                      "LATITUDE (HIGH ACCURACY)",
                      "LONGITUDE (HIGH ACCURACY)",
                      "WMO STATION NUMBER",
                      "WMO BLOCK NUMBER"]
    with open(filename, 'rb') as ins:
        decoder = pybufrkit.decoder.Decoder()
        for bufr_message in pybufrkit.decoder.generate_bufr_message(decoder, ins.read()):

            # Convert the BUFR message to JSON
            json_data = pybufrkit.renderer.NestedJsonRenderer().render(bufr_message)
            keyvalues = get_all_key_values(json_data)
            has_all = False
            for lookup_name in lookup_names:
                if lookup_name in keyvalues:
                    has_all = True
                    break
            for name in required_names:
                if name not in keyvalues:
                    has_all = False

            if has_all:
                for lookup_name in lookup_names:
                    if lookup_name in keyvalues:
                        values += [keyvalues[lookup_name]]
                        break
                lats += [keyvalues['LATITUDE (HIGH ACCURACY)']]
                lons += [keyvalues['LONGITUDE (HIGH ACCURACY)']]
                # Create the station ID by concatenating the block number and the station number
                id = keyvalues['WMO BLOCK NUMBER'] * 1000 + keyvalues['WMO STATION NUMBER']
                ids += [id]

    #return ids, lats, lons, values
    return np.array(ids), np.array(lats), np.array(lons), np.array(values)

def get_all_key_values(stuff):
    """ Recursively traverse a dictionary tree to find information. It looks for dictionaries with
    'description' and 'value' in them.
    """
    if isinstance(stuff, dict):
        if 'description' in stuff and 'value' in stuff and isinstance(stuff['value'], numbers.Number):
            return {stuff['description']: stuff['value']}
        ret = dict()
        for key, value in stuff.items():
            if isinstance(value, list):
                for key, value in get_all_key_values(value).items():
                    ret[key] = value
        return ret
    elif isinstance(stuff, list):
        ret = dict()
        for s in stuff:
            for key, value in get_all_key_values(s).items():
                ret[key] = value
        return ret
    else:
        return dict()

def print_contents(filename):
    decoder = pybufrkit.decoder.Decoder()
    with open(filename, 'rb') as ins:
        for bufr_message in pybufrkit.decoder.generate_bufr_message(decoder, ins.read()):
            json_data = pybufrkit.renderer.FlatTextRenderer().render(bufr_message)
            print(json_data)
            return

input_file = "syno_2021111200.bufr" #"/lustre/storeB/project/metproduction/products/obs_dec/rdb/syno/syno_2021120709.bufr"
ids, lats, lons, values = get(input_file)
print("Loaded %d stations" % len(ids))

print(len(ids))
print(len(lats))
print(len(lons))
print(len(values))
values = values - 273.15


#lons, lats = np.meshgrid(lons, lats)
opoints = gridpp.Points(lats, lons)


def get_apply_array(laf_array, laf_min, laf_max):
    apply_array = np.zeros(laf.shape)
    #print(apply_array.shape)
    for i in range(0, len(apply_array[:,0])):
        for j in range(0, len(apply_array[0,:])):
            if laf[i,j] >= laf_min and laf[i,j] <= laf_max:
                apply_array[i,j] = 1 
    return apply_array

file_geo = netCDF4.Dataset('/mnt/c/users/jmnip/met2021/yr_global_grid.nc') #ecmwf_global (1).nc     
file_global = netCDF4.Dataset('/mnt/c/users/jmnip/met2021/ec_hres_2021111200.nc')

temp = file_global.variables['air_temperature_2m'][60,0,: ,:] - 273.15
laf = file_geo.variables['land_area_fraction'][0,0,: ,: ]
elev = file_geo.variables['surface_geopotential'][0,0,:,:] / 9.81

lat = file_global.variables["latitude"]#[1230:1250]
lon = file_global.variables["longitude"]#[1980:2000]
lon, lat = np.meshgrid(lon, lat)

igrid = gridpp.Grid(lat, lon)
print(lats.shape)


temperature = file_global.variables["air_temperature_2m"][60, 0, 1230:1250, 1980:2000]

print(np.round(laf[:,0:14],2))

file_global.close()

#laf_min = 0
#laf_max = 0.85
#apply_array[laf >= laf_min & laf <= laf_max] = 1

#apply_array( laf >= 0 & laf <= 0.85).astype(int)

apply_array = get_apply_array(laf, 0, 0.85)

tFinal = gridpp.neighbourhood_search(temp, laf, 1, 0.7, 1, 0.1, apply_array)
print(tFinal.shape)

#gridpp.simple_gradient(grid, points, input, -0.0065)

t_diff = tFinal - temp 
"""
for i in range(len(ids)):
    print(str(i) + ': ' + str(format(ids[i], '05d'))  + ' \t ' +  str(format(lats[i], '0.5f')) + '\t ' +  str(format(lons[i], '03.5f')) + '\t ' +  str(np.round(values[i], 1)))
"""
raw = gridpp.nearest(igrid, opoints, temp)
output = gridpp.nearest(igrid, opoints, tFinal)
olat = gridpp.nearest(igrid, opoints, lat)
olon = gridpp.nearest(igrid, opoints, lon)
    
print('tFinal' + '\t', 'raw', 'lat', 'lon', 'raw_temp', 'lat', 'lon')
for i in range(len(output)):
    if output[i] - raw[i] != 0:
        if lats[i] < -33 and lats[i] > -35 and lons[i] < 20 and lons[i]  > 18:
            print('-------------')
        print(i), print('{:.1f}'.format(output[i]), '{:.1f}'.format(raw[i]), '{:.1f}'.format(output[i] - raw[i]), olat[i], olon[i], "{:.1f}".format(values[i]), lats[i], lons[i])
        if lats[i] < -33 and lats[i] > -35 and lons[i] < 20 and lons[i]  > 18:
            print('-------------')
fig, ax = plt.subplots(figsize= (5,5))

s = ax.scatter(lons, lats, c = values, cmap= 'RdBu_r', vmin = -30, vmax = 30, s = 5)
s.set_clim([-30,30])
cbar = fig.colorbar(s)
plt.savefig('stations_new.png')
plt.grid('on')
plt.close() 

# Hele verden ppå et tidsserie,
# Lag veriffil 
# MAE, obsfcst