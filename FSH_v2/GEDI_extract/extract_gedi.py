# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 20:35:04 2023

@author: 11913
"""

#!/usr/bin/env python3


import os
import h5py
import numpy as np
import pandas as pd
import geopandas as gp
from shapely.geometry import Point
import hdf5storage as h5ps
#import geoviews as gv
#from geoviews import opts, tile_sources as gvts
#import holoviews as hv
#gv.extension('bokeh', 'matplotlib')
import pdb


def runCmd(cmd):
    import subprocess
    out = subprocess.getoutput(cmd)
    return out

# Define a function for visualizing GEDI points
# def pointVisual(features, vdims):
#     return (gvts.EsriImagery * gv.Points(features, vdims=vdims).options(tools=['hover'], height=500, width=900, size=5,
#                                                                         color='yellow', fontsize={'xticks': 10, 'yticks': 10,
#                                                                                                   'xlabel':16, 'ylabel': 16}))

if __name__ == '__main__':
    '''
    Main driver.
    '''
    inDir = os.getcwd() + os.sep  # Set input directory to the current working directory
    
    gediFiles = [g for g in os.listdir() if g.endswith('.h5')]  # List GEDI L2A .h5 files in the inDir

    # Set up lists to store data
    shotNum, dem, zElevation, zHigh, zLat, zLon, rh95, rh98, rh100 ,quality ,degrade, sensitivity ,beamI = ([] for i in range(13))
    
    #define the spatial coverage of the data 
    
    extra_radius = 0.01;
    minLat = 38.445 - extra_radius;
    maxLat = 40.274  + extra_radius;
    minLon = 120.46  - extra_radius;
    maxLon = 122.63 + extra_radius; 
    
    matFile_out = 'gedi_crop.mat'

    lat_out = np.zeros(0)
    lon_out = np.zeros(0)
    rh100_out = np.zeros(0)
    DEM_out = np.zeros(0)
    LiDAR_z0_out = np.zeros(0)
    rh98_out = np.zeros(0)
    rh95_out = np.zeros(0)
    
    #masking the trees higher than this number
    h_threshold  = 55
    
    #    pdb.set_trace()
    for gediFile in gediFiles:
        
        print(gediFile) 
    
        gediL2A = h5py.File(gediFile, 'r')
        
        beamNames = [g for g in gediL2A.keys() if g.startswith('BEAM')]
        
        gediL2A_objs = []
        gediL2A.visit(gediL2A_objs.append)                                           # Retrieve list of datasets
        gediSDS = [o for o in gediL2A_objs if isinstance(gediL2A[o], h5py.Dataset)]  # Search for relevant SDS inside data file
        
        shotNum, dem, zElevation, zHigh, zLat, zLon, rh95, rh98, rh100 ,quality ,degrade, sensitivity ,beamI, solar_elev = ([] for i in range(14))
        
        
        # Loop through each beam and open the SDS needed
        for b in beamNames:
            if [g for g in gediSDS if g.endswith('/shot_number') and b in g].__len__() == 0:
                continue
                
            #if gediL2A[b].attrs['description'] == 'Coverage beam':
            #    continue
            [shotNum.append(h) for h in gediL2A[[g for g in gediSDS if g.endswith('/shot_number') and b in g][0]][()]]
            [dem.append(h) for h in gediL2A[[g for g in gediSDS if g.endswith('/digital_elevation_model') and b in g][0]][()]]
            [zElevation.append(h) for h in gediL2A[[g for g in gediSDS if g.endswith('/elev_lowestmode') and b in g][0]][()]]
            [zHigh.append(h) for h in gediL2A[[g for g in gediSDS if g.endswith('/elev_highestreturn') and b in g][0]][()]]
            [zLat.append(h) for h in gediL2A[[g for g in gediSDS if g.endswith('/lat_lowestmode') and b in g][0]][()]]
            [zLon.append(h) for h in gediL2A[[g for g in gediSDS if g.endswith('/lon_lowestmode') and b in g][0]][()]]
            [rh95.append(h[95]) for h in gediL2A[[g for g in gediSDS if g.endswith('/rh') and b in g][0]][()]]
            [rh98.append(h[98]) for h in gediL2A[[g for g in gediSDS if g.endswith('/rh') and b in g][0]][()]]
            [rh100.append(h[100]) for h in gediL2A[[g for g in gediSDS if g.endswith('/rh') and b in g][0]][()]]
            [quality.append(h) for h in gediL2A[[g for g in gediSDS if g.endswith('/quality_flag') and b in g][0]][()]]
            [degrade.append(h) for h in gediL2A[[g for g in gediSDS if g.endswith('/degrade_flag') and b in g][0]][()]]
            [sensitivity.append(h) for h in gediL2A[[g for g in gediSDS if g.endswith('/sensitivity') and b in g][0]][()]]
            [solar_elev.append(h) for h in gediL2A[[g for g in gediSDS if g.endswith('/solar_elevation') and b in g][0]][()]]
            [beamI.append(h) for h in [b] * len(gediL2A[[g for g in gediSDS if g.endswith('/shot_number') and b in g][0]][()])]

        # Convert lists to Pandas dataframe
        allDF = pd.DataFrame({'Shot Number': shotNum, 'Beam': beamI, 'Latitude': zLat, 'Longitude': zLon, 'radar_dem': dem,
                         'Elevation (m)': zElevation, 'Canopy Elevation (m)': zHigh, 'Canopy Height (rh100)': rh100, 'RH 98': rh98,
                         'RH 95': rh95, 'Quality Flag': quality, 'Degrade Flag': degrade, 'Sensitivity': sensitivity,'solar_elevation':solar_elev})

        #del beamI, degrade, dem, gediSDS, rh100, rh98, rh95, quality, sensitivity, zElevation, zHigh, zLat, zLon, shotNum

    
        allDF = allDF.where(allDF['Latitude'] > minLat)
        allDF = allDF.where(allDF['Latitude'] < maxLat)
        allDF = allDF.where(allDF['Longitude'] > minLon)
        allDF = allDF.where(allDF['Longitude'] < maxLon)
        #allDF = allDF.where(allDF['solar_elevation'] <0)

        allDF = allDF.dropna()  # Drop shots outside of the ROI

        # Set any poor quality returns to NaN
        allDF = allDF.where(allDF['Quality Flag'].ne(0))
        allDF = allDF.where(allDF['Degrade Flag'].ne(1))
        allDF = allDF.where(allDF['Sensitivity'] > 0.95)
        allDF = allDF.where(abs(allDF['radar_dem'] - allDF['Elevation (m)']) <= h_threshold )
    
        allDF = allDF.dropna()
    
        lat_out =  np.append(lat_out,np.array(allDF['Latitude']))
        lon_out = np.append(lon_out,np.array(allDF['Longitude']))
        rh100_out = np.append(rh100_out,np.float64(np.array(allDF['Canopy Height (rh100)'])))
        rh98_out = np.append(rh98_out,np.float64(np.array(allDF['RH 98'])))
        rh95_out = np.append(rh95_out,np.float64(np.array(allDF['RH 95'])))
        DEM_out = np.append(DEM_out,np.float64(np.array(allDF['radar_dem'])))
        LiDAR_z0_out = np.append(LiDAR_z0_out,np.float64(np.array(allDF['Elevation (m)'])))


    # Take the lat/lon dataframe and convert each lat/lon to a shapely point
    #allDF['geometry'] = allDF.apply(lambda row: Point(row.Longitude, row.Latitude), axis=1)

    # Convert to geodataframe
    #allDF = gp.GeoDataFrame(allDF)
    #allDF = allDF.drop(columns=['Latitude','Longitude'])

    #allDF['Shot Number'] = allDF['Shot Number'].astype(str)  # Convert shot number to string
   

    h5ps.savemat(matFile_out, {'lat': lat_out, 'lon': lon_out, 'rh100': rh100_out, 'rh98': rh98_out,
                               'rh95':rh95_out,'radar_dem':DEM_out, 'Lidar_z0':LiDAR_z0_out }, 
                               format='7.3', oned_as='column',store_python_metadata=True)
    
    
    
    #vdims = []
    #for f in allDF:
    #    if f not in ['geometry']:
    #        vdims.append(f)

    #visual = pointVisual(allDF, vdims = vdims)
    #visual
#    visual * gv.Polygons(redwoodNP).opts(line_color='red', color=None)

    # Plot the basemap and geoviews Points, defining the color as the Canopy Height for each shot
    #(gvts.EsriImagery * gv.Points(allDF, vdims=vdims).options(color='Canopy Height (rh100)',cmap='plasma', size=3, tools=['hover'],
    #                                                          clim=(0,102), colorbar=True, clabel='Meters',
    #                                                          title='GEDI Canopy Height over Redwood National Park: June 19, 2019',
    #                                                          fontsize={'xticks': 10, 'yticks': 10, 'xlabel':16, 'clabel':12,
    #                                                          'cticks':10,'title':16,'ylabel':16})).options(height=500,
    #                                                                                                        width=900)

    #(gvts.EsriImagery * gv.Points(allDF, vdims=vdims).options(color='Elevation (m)',cmap='terrain', size=3, tools=['hover'],
    #                                                          clim=(min(allDF['Elevation (m)']), max(allDF['Elevation (m)'])),
    #                                                          colorbar=True, clabel='Meters',
    #                                                          title='GEDI Elevation over Redwood National Park: June 19, 2019',
    #                                                          fontsize={'xticks': 10, 'yticks': 10, 'xlabel':16, 'clabel':12,
    #                                                          'cticks':10,'title':16,'ylabel':16})).options(height=500,
    #                                                                                                        width=900)

    # outName = gediL2A.filename.replace('.h5', '.json')  # Create an output file name using the input file name

    # allDF.to_file(outName, driver='GeoJSON')  # Export to GeoJSON
    del allDF
    

