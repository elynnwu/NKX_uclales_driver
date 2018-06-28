#!/usr/bin/python

import sys, getopt
import os
import prep_uclales_methods as prep_les
import pandas as pd
import numpy as np
def main(argv):
    if len(argv)!=3:
        sys.exis('Input arguments must be year, month, day!')
    year = argv[0]
    month = argv[1]
    day = argv[2]
    home_dir = '/mnt/lab_45d1/database/Sc_group/'
    '''Input a date'''
    date = pd.date_range(year+'-'+month+'-'+day+' 12:00', year+'-'+month+'-'+day+' 12:00')[0]
    caseName = 'NKX_'+date.strftime('%Y%m%d')
    try: 
        os.makedirs(home_dir+'uclales_output/'+caseName)
        os.makedirs(home_dir+'uclales_output/'+caseName+'/hr_1_3_60min_with_spinup')
        os.makedirs(home_dir+'uclales_output/'+caseName+'/hr_3_4_1min')
        os.makedirs(home_dir+'uclales_output/'+caseName+'/hr_4_18_60min')
    except OSError:
        if not os.path.isdir('./uclales_output/'+caseName):
            raise
    '''Hour minus 2 to 3, saving hourly data'''
    os.chdir(home_dir+'uclales_output/'+caseName+'/hr_1_3_60min_with_spinup')
    '''Make sound_in file'''
    IC, uwind, vwind, z = prep_les.make_sound_in_file(home_dir,date)
    '''Make backrad_in file'''
    prep_les.make_backrad_in_file(home_dir,date,IC,z)
    '''Make zm_grid and zt_grid files'''
    domain_H = np.max([IC['z_inv_base']*2,1000.])
    zm, zt = prep_les.write_z_grid_relax_transition(IC['z_inv_base'],domain_H)
    '''Make NAMELIST'''
    prep_les.write_NAMELIST(nzp=len(zm),timmax=18000.,runtype='INITIAL',frqanl=3600.,\
                   filprf=caseName,hfilin=caseName+'.rst',strtim=float(date.dayofyear),\
                   dthcon=IC['SHF'],drtcon=IC['LHF'],\
                   th00=IC['eq_thetalBL'],umean=uwind.mean(),vmean=vwind.mean(),div=IC['D_SEAarea'],fr0=IC['F0'],fr1=IC['F1'],xka=IC['xka'])
    '''Run UCLALES'''
    os.system('ln -s /home/elw014/radtyp2/uclales/build/uclales .')
    os.system('ln -s /home/elw014/radtyp2/uclales/bin/datafiles .')
    os.system('mpiexec -n 8 ./uclales > log')
    
    '''Hour 3 - 4, saving minute data'''
    os.chdir(home_dir+'uclales_output/'+caseName+'/hr_3_4_1min')
    prep_les.write_NAMELIST(nzp=len(zm),timmax=21600.,runtype='HISTORY',frqanl=60.,\
                   filprf=caseName,hfilin=caseName+'.rst',strtim=float(date.dayofyear),\
                   dthcon=IC['SHF'],drtcon=IC['LHF'],\
                   th00=IC['eq_thetalBL'],umean=uwind.mean(),vmean=vwind.mean(),div=IC['D_SEAarea'],fr0=IC['F0'],fr1=IC['F1'],xka=IC['xka'])
    os.system('cp ../hr_1_3_60min_with_spinup/zm_grid_in .')
    os.system('cp ../hr_1_3_60min_with_spinup/zt_grid_in .')
    os.system('cp ../hr_1_3_60min_with_spinup/sound_in .')
    os.system('cp ../hr_1_3_60min_with_spinup/backrad_in .')
    os.system('cp ../hr_1_3_60min_with_spinup/*.rst .')
    os.system('ln -s /home/elw014/radtyp2/uclales/build/uclales .')
    os.system('ln -s /home/elw014/radtyp2/uclales/bin/datafiles .')
    os.system('mpiexec -n 8 ./uclales > log')
    prep_les.SEND_ALERT_EMAIL(caseName, ' Hour 3 - 4')
    
    '''Hour 4 - 18, saving hourly data'''
    os.chdir(home_dir+'uclales_output/'+caseName+'/hr_4_18_60min')
    prep_les.write_NAMELIST(nzp=len(zm),timmax=72000.,runtype='HISTORY',frqanl=3600.,\
                   filprf=caseName,hfilin=caseName+'.rst',strtim=float(date.dayofyear),\
                   dthcon=IC['SHF'],drtcon=IC['LHF'],\
                   th00=IC['eq_thetalBL'],umean=uwind.mean(),vmean=vwind.mean(),div=IC['D_SEAarea'],fr0=IC['F0'],fr1=IC['F1'],xka=IC['xka'])
    os.system('cp ../hr_3_4_1min/zm_grid_in .')
    os.system('cp ../hr_3_4_1min/zt_grid_in .')
    os.system('cp ../hr_3_4_1min/sound_in .')
    os.system('cp ../hr_3_4_1min/backrad_in .')
    os.system('cp ../hr_3_4_1min/*.rst .')
    os.system('ln -s /home/elw014/radtyp2/uclales/build/uclales .')
    os.system('ln -s /home/elw014/radtyp2/uclales/bin/datafiles .')
    os.system('mpiexec -n 8 ./uclales > log')
    prep_les.SEND_ALERT_EMAIL(caseName, ' Hour 4 - 18')
if __name__ == "__main__":
   main(sys.argv[1:])
