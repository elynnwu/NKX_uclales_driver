#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Stitch the raw files from ucla-les from different processors
This script has only been tested on nx = 2 and ny = 4, use with caution
@author: elynn
"""
import sys, getopt
import os
from netCDF4 import Dataset
import numpy as np
from scipy.io import netcdf

def main(argv):
    '''Configure here'''
    file_dir = ''
    output_prefix = argv[0]
    output_dir = ''
    nx = int(argv[3])#2 #number of processors in x-direction
    ny = int(argv[4])#4 #number of processors in y-direction
    x_pts = int(argv[5])#48 #number of points in x-direciton for each processor
    y_pts = int(argv[6])#24 #number of points in y-direciton for each processor
    z_pts = int(argv[1]) #number of points in z-direction
    t_pts = int(argv[2]) #number of time index 
    variables = ['l','w','t','q','u','v','rflx','p','lwfd','lfwu'] #variables to be stitched
    '''End configuration'''

    domain_var = {} #dictionary to store all variables
    for i in range(len(variables)): #initialize the stitched variables
        domain_var[variables[i]] = np.zeros((t_pts,y_pts*ny,x_pts*nx,z_pts))
    xm_all = []
    ym_all = []
    for x in range(nx):
        for y in range(ny):
            filename = output_prefix+'.'+str(x).zfill(4)+str(y).zfill(4)+'.nc' #filename from each processor
            current = Dataset(filename)
            for i in range(len(variables)):
                domain_var[variables[i]][:,y*y_pts:y*y_pts+y_pts,x*x_pts:x*x_pts+x_pts,:] = current[variables[i]][:,:,:,:]
            zt = current['zm'][:]
            xm_all.append(current['xm'][:])
            ym_all.append(current['ym'][:])
            current.close()

    x = sorted(np.unique(np.array(xm_all))) #get x grids from processors
    y = sorted(np.unique(np.array(ym_all))) #get y grids from processors
    time = np.arange(1,t_pts+1)

    f = Dataset(output_dir+output_prefix+'_stitched.nc', 'w')
    f.createDimension('t_index', t_pts)
    f.createDimension('x', nx*x_pts)
    f.createDimension('y', ny*y_pts)
    f.createDimension('z', z_pts)
    times = f.createVariable('time', 'f4', ('t_index',))
    times[:] = time[:]
    times.units = 'minutes from hr 3'
    xcoord = f.createVariable('xm', 'f4', ('x',))
    ycoord = f.createVariable('ym', 'f4', ('y',))
    zcoord = f.createVariable('zm', 'f4', ('z',))
    xcoord[:] = x[:]
    ycoord[:] = y[:]
    zcoord[:] = zt[:]
    temps = {} #initialize netCDF4 variables to be written
    for i in range(len(variables)):
        print variables[i]
        temps['temp'+str(i+1)] = f.createVariable(variables[i], 'f4', ('t_index','y', 'x', 'z'))
        temps['temp'+str(i+1)][:,:,:,:] = domain_var[variables[i]][:,:,:,:]
    f.close()
if __name__ == "__main__":
   main(sys.argv[1:])
