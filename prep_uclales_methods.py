import codecs
import numpy as np
from scipy import interpolate
import pandas as pd
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText
import os 
import smtplib
def write_z_grid(ibh,dzrat=1.05,domain_H=2500.,dz_ct=5.0,dz_bl=15.0):
    zm = []
    ibh = int(round(ibh/5.0)*5.0) #round to closest 5
    zz = 0.
    while (zz<ibh-100.): #below cloud
        zm.append(zz)
        zz = zz + dz_bl
    zz = ibh-100.
    while (zz>=ibh-100.)&(zz<ibh+100.): #fine resolution below/above cloud top
        zm.append(zz)
        zz = zz + dz_ct
    dz = dz_bl
    while(zz<=domain_H):
        dz = dz*dzrat
        zz = zz+dz
        zm.append(zz)
    zt = []
    zt.append(zm[0]-dz_bl/2.)
    for i in range(len(zm)-1):
        zt.append(0.5*zm[i]+0.5*zm[i+1])
    zm = np.array(zm)
    zt = np.array(zt)
    np.savetxt('zm_grid_in',zm,fmt='%5.1f')
    np.savetxt('zt_grid_in',zt,fmt='%5.1f')
    return zm, zt

def write_z_grid_relax_transition(ibh,domain_H,dzrat=1.095,dz_ct=5.0,dz_bl=15.0):
    padding = 0
    counter = 0
    dz = dz_bl
    while (dz > dz_ct):
        dz = dz/dzrat
        padding = padding + dz
    if dz < dz_ct: #shrink one too much
        padding = padding - dz
    zm = []
    zz = 0.
    while (zz<(ibh-50-padding)): #below cloud
        zm.append(zz)
        zz = zz + dz_bl
    dz = dz_bl
    while (zz<ibh-50.):
        zm.append(zz)
        dz = dz/dzrat
        zz = zz + dz
    while (zz>=ibh-50.) & (zz<=ibh+200.): #fine resolution below/above cloud top
        zm.append(zz)
        zz = zz + dz_ct
    dz = dz_ct
    while (dz < dz_bl):
        zm.append(zz)
        dz = dz * dzrat
        zz = zz + dz   
    dz = dz_bl
    while(zz<=domain_H):
        zm.append(zz)
        dz = dz*dzrat
        zz = zz+dz
    zt = []
    zt.append(zm[0]-dz_bl/2.)
    for i in range(len(zm)-1):
        zt.append(0.5*zm[i]+0.5*zm[i+1])
    zm = np.array(zm)
    zt = np.array(zt)
    np.savetxt('zm_grid_in',zm,fmt='%5.1f')
    np.savetxt('zt_grid_in',zt,fmt='%5.1f')
    return zm, zt
def temperature_above_inv(thetaL_BL,dthetaL,z,zi):
    thetaL_above = np.zeros(len(z))
    for i in range(len(z)):
        thetaL_above[i] = thetaL_BL+dthetaL+(z[i]-zi)**(1./3.)
    return thetaL_above
def make_sound_in_file(home_dir,date):
    sounding = pd.read_csv(home_dir+'NKX_sounding/72293_'+date.strftime('%Y_%m_%d')+'_12Z.csv')
    IC = pd.read_csv(home_dir+'NKX_sounding/FinalTable.csv',index_col=0,parse_dates=True).loc[date]
    '''Get vertical coordinates'''
    z_BL = sounding['HGHT [m]'][(sounding['HGHT [m]']>137.) & (sounding['HGHT [m]']<IC['z_inv_base'])]
    z_above = sounding['HGHT [m]'][(sounding['HGHT [m]']>=IC['z_inv_base']) & (sounding['HGHT [m]']<3000.)]
    z = np.array([0.] + list(z_BL) + [IC['z_inv_base']-1] + list(z_above))
    '''Idealize the sounding'''
    thetaL_above = temperature_above_inv(IC['eq_thetalBL'],IC['eq_dthetal'],z[z>=IC['z_inv_base']],IC['z_inv_base'])
    tl = np.concatenate((IC['eq_thetalBL']*np.ones(len(z[z<IC['z_inv_base']])),list(thetaL_above)))
    qt = (IC['eq_qtBL']+IC['eq_dqt'])*np.ones(len(z))
    qt[z<IC['z_inv_base']] = IC['eq_qtBL']
    '''Get wind component'''
    usounding = -(sounding['SKNT [knot]']*0.514)*np.sin(sounding['DRCT [deg]']*np.pi/180.)
    vsounding = -(sounding['SKNT [knot]']*0.514)*np.cos(sounding['DRCT [deg]']*np.pi/180.)
    fz = z<3000.
    uwind = np.array([0.]+list(usounding[(sounding['HGHT [m]']>137.) & (sounding['HGHT [m]']<=IC['z_inv_base'])])+\
                    list(usounding[(sounding['HGHT [m]']>=IC['z_inv_base']) & (sounding['HGHT [m]']<3000.)]))
    vwind = np.array([0.]+list(vsounding[(sounding['HGHT [m]']>137.) & (sounding['HGHT [m]']<=IC['z_inv_base'])])+\
                    list(vsounding[(sounding['HGHT [m]']>=IC['z_inv_base']) & (sounding['HGHT [m]']<3000.)]))
    output = np.stack((z,tl,qt,uwind,vwind),axis=-1)
    output[0,0] = IC['PSFC']/100.
    np.savetxt('sound_in',output,fmt='%10.3f')
    return IC, uwind, vwind, z
def make_backrad_in_file(home_dir,date,IC,z):
    sounding = pd.read_csv(home_dir+'NKX_sounding/72293_'+date.strftime('%Y_%m_%d')+'_12Z.csv')
    sounding = sounding.dropna()
    temp = sounding['TEMP [C]'].dropna().as_matrix()+273.15
    pres = sounding['PRES [hPa]'].dropna().as_matrix()
    p = interpolate.interp1d(pres[0:2],temp[0:2],fill_value='extrapolate')
    tsrf = p(IC['PSFC']/100.) #extrapolate to get surface temperature
    i_pinv = np.where(z==IC['z_inv_base'])[0][0]
    i_3km = np.where(z<3000.)[0][-1]
    temp = np.array([tsrf.tolist()]+list(temp))
    pres = np.array([IC['PSFC']/100.]+list(pres))
    watr = np.array(list(IC['eq_qtBL']*np.ones(i_pinv+1)) + \
                    list((IC['eq_qtBL']+IC['eq_dqt'])*np.ones(i_3km-i_pinv)) + \
                    list(sounding['MIXR [g/kg]'][i_3km:]))/1000.
    ozone_ref = np.genfromtxt(home_dir+'NKX_sounding/backrad_in_ref',skip_header=1)
    pozon = interpolate.interp1d(ozone_ref[:,0],ozone_ref[:,3])
    if pres[0]>ozone_ref[:,0][-1]: #quick hack to make sure interpolation is within reference points
        pres[0] = ozone_ref[:,0][-1]
    ozone = pozon(pres)
    backrad_output = np.stack((pres[::-1],temp[::-1],watr[::-1],ozone[::-1],np.zeros(len(temp))),axis=-1)
    np.savetxt('backrad_in',backrad_output,fmt='%15.3f %8.2f %12.6f %15.10f %8.3f',header=str(tsrf)+' '+str(len(temp)),comments='')
def write_NAMELIST(nzp,timmax,runtype,frqanl,filprf,hfilin,strtim,dthcon,drtcon,th00,umean,vmean,div,fr0,fr1,\
                   nxp=100,nyp=100,igrdtyp=-3,deltax=35.,deltay=35.,deltaz=10.,\
                   nxpart='.true.',dtlong=2.,distim=100.,level=2,CCN=55.0e6,prndtl=-0.33333,\
                   ssam_intvl=15.,savg_intvl=900.,corflg='.true.',iradtyp=2,ubmin=-0.25,):
    '''Write UCLALES NAMELIST
    Default variables are at the end of the input vars
    05/29/18: fixing radiation type to 2, dynamic input for fr0 and fr1
	      subsidence cooling added to keep the tropospheric temperature profile constant
    '''
    with codecs.open('NAMELIST','w',encoding='ascii') as f:
        f.write('&model\n')
        f.write(' nxp = '+str(nxp)+'\n')
        f.write(' nyp = '+str(nyp)+'\n')
        f.write(' nzp = '+str(nzp)+'\n')
        f.write(' deltax = '+str(deltax)+'\n')
        f.write(' deltay = '+str(deltay)+'\n')
        f.write(' deltaz = '+str(deltaz)+'\n')
        f.write(' nxpart = '+nxpart+'\n')
        f.write(' igrdtyp = 3\n')
        f.write(' dtlong = '+str(dtlong)+'\n')
        f.write(' distim = '+str(distim)+'\n')
        f.write(' timmax = '+str(timmax)+'\n')
        f.write(' runtype = "'+str(runtype)+'"\n')
        f.write(' frqanl = '+str(frqanl)+'\n')
        f.write(' level = '+str(level)+'\n')
        f.write(' CCN = '+str(CCN)+'\n')
        f.write(' prndtl = '+str(prndtl)+'\n')
        f.write(' filprf = \''+str(filprf)+'\'\n')
        f.write(' hfilin = \''+str(hfilin)+'\'\n')
        f.write(' ssam_intvl = '+str(ssam_intvl)+'\n')
        f.write(' savg_intvl = '+str(savg_intvl)+'\n')
        f.write(' corflg = '+corflg+'\n')
        f.write(' strtim = '+str(strtim-0.083)+'\n') #start simulation 2 hours before midnight to account for spinup time
        f.write(' iradtyp = '+str(iradtyp)+'\n')
        f.write(' dthcon = '+str(dthcon)+'\n')
        f.write(' drtcon = '+str(drtcon)+'\n')
        f.write(' ubmin = '+str(ubmin)+'\n')
        f.write(' th00 = '+str(th00)+'\n')
        f.write(' umean = '+str(umean)+'\n')
        f.write(' vmean = '+str(vmean)+'\n')
        f.write(' div = '+str(np.abs(div))+'\n')
	f.write(' fr0 = '+str(fr0)+'\n')
	f.write(' fr1 = '+str(fr1)+'\n')        
	f.write('/\n')
        f.close()

def SEND_ALERT_EMAIL(caseName,timePeriod):
    '''Send an email once the model completes
    '''
    fromaddr = "ucsdwrf@gmail.com"
    toaddr = ["elw014@eng.ucsd.edu"]
    msg = MIMEMultipart()
    msg['From'] = fromaddr
    msg['To'] = ", ".join(toaddr)
    msg['Subject'] = 'UCLALES COMPLETED - ' + caseName + timePeriod
    hostname = os.popen("uname -n").read()
    body = "LES completed on " + hostname   
    msg.attach(MIMEText(body, 'plain'))
    server = smtplib.SMTP('smtp.gmail.com', 587)
    server.starttls()
    server.login(fromaddr, "ucsdwrfforecast")
    text = msg.as_string()
    server.sendmail(fromaddr, toaddr, text)
    server.quit()
