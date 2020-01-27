#!/usr/bin/env python
# radiative kernels
from netCDF4 import Dataset
from matplotlib import pyplot as plt
from matplotlib import rc
import matplotlib.colors as mc
import matplotlib as mpl
import numpy as np
import climlab
from climlab import domain
from climlab.domain import field
import matplotlib.colors as mc
from matplotlib import ticker, cm

# default color cycle
# print (plt.rcParams['axes.prop_cycle'].by_key()['color'])
#['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
def_colors=plt.rcParams['axes.prop_cycle'].by_key()['color'][:]

H=7.e3 # [m]
pr=1000. # [mb]
a=6371.e3 # [m]
#omg=7.292e-5 # [s-1]
Rd=287.04 # [m2s-2K-1]
Rv=461.52 #J/(kgK)
kappa=2./7.
g=9.8 # [m/s2]
#Tr=g*H/R
#rho=100.*pr/R/Tr*exp(-z/H)
##N2=5.e-4 # [s-2]
##S=H*N2/R
Mair=29. #[g] molar mass of air
Mo3=48. #[g] molar mass of ozone
Lv=2.50000e6 # [Jkg-1]
Ls=3.34e5 # [Jkg-1]
cp=Rd/kappa # [Jkg-1K-1]
stefan=5.670373e-8 # [Wm-2K-4]
spd=86400.

# formula 3: Goff-Gratch
def esat(TK):
    log10_ew =  -7.90298*(373.16/TK-1.) + 5.02808*np.log10(373.16/TK) \
    -1.3816e-7*(10**(11.344*(1.-TK/373.16))-1.) \
    +8.1328e-3*(10**(-3.49149*(373.16/TK-1))-1.) + np.log10(1013.246)
    return 10**log10_ew*100. # Pa

data_dir='/mnt/storage/Research/data/HIRAM/'
time_domain='month'
directory=data_dir+'slab-C400-cre0-uw-a30.2'+'/history/'
filename=directory+'19820101.atmos_month.nc'
ncfile = Dataset(filename,'r')
var_list=ncfile.variables.keys()
p=ncfile.variables['pfull'][:]
phalf=ncfile.variables['phalf'][:]
nr=0 # global, tropics, extra-tropics
lat=ncfile.variables['lat'][:]
rlat=np.pi*lat/180. # deg to rad
z=H*log(pr/p) # [km] log pressure coordinate
zhalf=H*log(pr/phalf) # [km] log pressure coordinate
nz = np.size(p)
ny = np.size(lat)
if nr==0:
    ny1=0
    ny2=ny
ifact=1./np.sum(np.cos(rlat[ny1:ny2]))

directory=data_dir+'slab-C400-cre0-uw-a30.7-S1365-oh1e8-diur'+'/history/'
filename=directory+'lon_mean_'+'19990101'+'_'+'month'+'.nc'
ncfile = Dataset(filename,'r')
qo3=np.mean(ncfile.variables['qo3'][:,:,:,0],axis=0)
ch4=1.557889297535934e-6
n2o=3.007494866529774e-7
o2=0.
ccl4=0.
cfc11=1.680e-10
cfc12=2.850e-10
cfc113=1.737268993839836e-11
cfc22=4.8743326488706363-11

#Gas value is taken from timeseries at time: 1980 Jan 01 00:00:00
 #ch4 value is   1.557889297535934E-006
#NOTE from PE    0: radiative_gases_mod: PROCESSING TIMESERIES FOR n2o
#Gas value is taken from timeseries at time: 1980 Jan 01 00:00:00
 #n2o value is   3.007494866529774E-007
#NOTE from PE    0: radiative_gases_mod: PROCESSING TIMESERIES FOR f11
#Gas value is taken from timeseries at time: 1980 Jan 01 00:00:00
 #f11 value is   1.680000000000000E-010
#NOTE from PE    0: radiative_gases_mod: PROCESSING TIMESERIES FOR f12
#Gas value is taken from timeseries at time: 1980 Jan 01 00:00:00
 #f12 value is   2.850000000000000E-010
#NOTE from PE    0: radiative_gases_mod: PROCESSING TIMESERIES FOR f113
#Gas value is taken from timeseries at time: 1980 Jan 01 00:00:00
 #f113 value is   1.737268993839836E-011
#NOTE from PE    0: radiative_gases_mod: PROCESSING TIMESERIES FOR f22
#Gas value is taken from timeseries at time: 1980 Jan 01 00:00:00
 #f22 value is   4.874332648870636E-011
 #Ozone data is obtained from a clim_zonal ozone file for year 1990

# define tropopause
#dT=np.zeros((nz,ny))
#for k in range(0,ncase):
    #for j in range(0,ny):
        ##nn=temp[k,:,j].argmin(0) # should not use this, problems at poles sometimes
        ##ptrop[j]=p[nn]
        #dT[:,j]=np.where(temp_2xco2[k,:,j]-temp[k,:,j]>0.,temp_2xco2[k,:,j]-temp[k,:,j],np.nan*np.zeros(nz))

#plt.figure()
#plt.contourf(lat,p,dT)
#plt.ylim(1000.,0.)

lat0=np.linspace(-90.,90.,19)
pt0=np.array([300.,285.,280.,270.,240.,195.,100.,100.,100.,100.,
    100.,100.,100.,195.,240.,270.,280.,285.,300.])
pt=numpy.interp(lat, lat0, pt0)

ntrop=np.zeros(ny, dtype=int16)
for j in range(0,ny):
    ntrop[j]=np.abs(p[:]-pt[j]).argmin(0)

#plt.figure()
#plt.plot(lat0,pt0)
#plt.plot(lat,pt)
#plt.plot(lat,p[ntrop])
#plt.ylim(290.,95.)

dp=np.zeros(nz)
for j in range(0,nz):
    dp[j]=phalf[j+1]-phalf[j]

# function to get global value of variable
def gb(var):
    if len(np.shape(var))==2:
        var_iz=np.zeros(ny)
        for j in range(0,ny):
            var_iz[j]=sum(var[ntrop[j]:,j]*dp[ntrop[j]:])
            #var_iz[j]=sum(var[:,j]*dp[:])            
    elif len(np.shape(var))==1:
        var_iz=var
    return sum(var_iz*np.cos(rlat))*ifact


# function to get symmetric value of variable
def sym_lat(var):
    var_sym=np.zeros_like(var)
    if len(np.shape(var))==4:
        for j in range(0,ny):
            var_sym[:,:,:,j] = 0.5*(var[:,:,:,j]+var[:,:,:,ny-1-j])
    elif len(np.shape(var))==3:
        for j in range(0,ny):
            var_sym[:,:,j] = 0.5*(var[:,:,j]+var[:,:,ny-1-j])
    elif len(np.shape(var))==2:
        for j in range(0,ny):
            var_sym[:,j] = 0.5*(var[:,j]+var[:,ny-1-j])
    elif len(np.shape(var))==1:
        for j in range(0,ny):
            var_sym[j] = 0.5*(var[j]+var[ny-1-j])
    return var_sym


###oh1e7###

#casename='cre0-uw-a30.1-S1365-oh1e7'
#casename='cre0-uw-a30.5-S1365-oh1e7-dail'
#casename='cre0-uw-a31.0-S1365-oh1e7-diur'

#casename='cre1-uw-a04.8-S1365-oh1e7'
#casename='cre1-uw-a00.6-S1365-oh1e7-dail'
#casename='cre1-uw-a12.6-S1365-oh1e7-diur'


###oh1e8###
#annual-mean
#casename='cre0-uw-a30.3-S1365-oh1e8'
#name='cre0-uw'

#casename='cre0-uw-a30.3-S1365-oh1e8-inst'
#name='cre0-uw'

#casename='cre1-uw-a04.8-S1365-oh1e8'
#name='2x-cre1-cins' #constant insolation (no time-variations)
#name='4x-cre1-cins' #constant insolation (no time-variations)
#name='8x-cre1-cins' #constant insolation (no time-variations)

#casename='cre1-uw-a32.1-S1365-oh1e8-l0'
#casename='cre1-uw-a00.0-S1377-oh1e8-i0'

#diur
#casename='cre0-uw-a30.7-S1365-oh1e8-diur'
#name='4x-cre0-uw'
#name='2x-cre0-uw'
#name='8x-cre0-uw'

#casename='cre0-uwl-a28.3-S1365-oh1e8-diur'
#name='cre0-lin'

#casename='cre0-uwe0-a30.6-S1365-oh1e8-diur'
#name='cre0-ce0'

#casename='cre0-uwve0-a30.5-S1365-oh1e8-diur'
#name='cre0-ve0'

#casename='cre1-uw-a17.3-S1365-oh1e8-diur'
#name='4x-cre1-uw'
#name='2x-cre1-uw'
#name='8x-cre1-uw'

#casename='cre1-uwl-a13.9-S1365-oh1e8-diur'
#name='4x-cre1-uwl'
#name='2x-cre1-uwl'
#name='8x-cre1-uwl'

casename='cre1-uwe0-a02.5-S1365-oh1e8-diur'
name='4x-cre1-ce0'
#name='2x-cre1-ce0'
#name='8x-cre1-ce0'

#casename='cre1-uwve0-a12.8-S1365-oh1e8-diur'
#name='cre1-ve0'

#casename='cre1-uw-a34.1-S1365-oh1e8-l0-diur'
#name='lrad-ctl'

#casename='cre1-uw-a10.7-S1365-oh1e8-i0-diur'
#name='irad-ctl'

linestyles=['-','-','-','-','-','-','-','-']
#colors=['b','b','r','r',
#'b','b','r','r',
#'g','g','m','m',
#'g','g','m','m']

# default color cycle
# print (plt.rcParams['axes.prop_cycle'].by_key()['color'])
#['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', 
#'#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
def_colors=plt.rcParams['axes.prop_cycle'].by_key()['color'][:]
colors=[def_colors[0],def_colors[1],def_colors[2],def_colors[3],
        def_colors[0],def_colors[1],def_colors[2],def_colors[3]]

ncase=2 # T400-C400, T400-C800

# domain creation for rad calc
##sfc, atm = domain.zonal_mean_column(num_lat=ny, num_lev=nz)
sfc, atm = domain.zonal_mean_column(lat=lat, lev=p)
# radiative gases
o3 = field.Field(qo3.T*Mair/Mo3, domain=atm)
absorber_C100={'O3': o3, 'CO2': 100.e-6, 'CH4':ch4, 'N2O':n2o, 'O2': o2,'CCL4':ccl4, 
    'CFC11':cfc11, 'CFC12':cfc12, 'CFC113':cfc113, 'CFC22':cfc22}
absorber_C200={'O3': o3, 'CO2': 200.e-6, 'CH4':ch4, 'N2O':n2o, 'O2': o2,'CCL4':ccl4, 
    'CFC11':cfc11, 'CFC12':cfc12, 'CFC113':cfc113, 'CFC22':cfc22}
absorber_C400={'O3': o3, 'CO2': 400.e-6, 'CH4':ch4, 'N2O':n2o, 'O2': o2,'CCL4':ccl4, 
    'CFC11':cfc11, 'CFC12':cfc12, 'CFC113':cfc113, 'CFC22':cfc22}
absorber_C800={'O3': o3, 'CO2': 800.e-6, 'CH4':ch4, 'N2O':n2o, 'O2': o2,'CCL4':ccl4, 
    'CFC11':cfc11, 'CFC12':cfc12, 'CFC113':cfc113, 'CFC22':cfc22}

# get surface albedo
nn=casename.find('-a')
alb_sfc= float(casename[nn+2:nn+6])*0.01 # fraction not percent
albedo = field.Field(np.ones(ny)*alb_sfc, domain=sfc)            

# hiram output variables
olr=np.zeros((ncase,ny))
olr_clr=np.zeros((ncase,ny))
swdn_toa=np.zeros((ncase,ny))
swup_toa=np.zeros((ncase,ny))
swup_toa_clr=np.zeros((ncase,ny))
swdn_sfc=np.zeros((ncase,ny))
swup_sfc=np.zeros((ncase,ny))
lwdn_sfc=np.zeros((ncase,ny))
lwup_sfc=np.zeros((ncase,ny))

swdn_sfc_clr=np.zeros((ncase,ny))
swup_sfc_clr=np.zeros((ncase,ny))
lwdn_sfc_clr=np.zeros((ncase,ny))
lwup_sfc_clr=np.zeros((ncase,ny))

precip=np.zeros((ncase,ny))

#year=['19820101']
#year=['19850101']
year=['19820101','19830101','19840101','19850101','19860101',
    '19870101','19880101','19890101','19900101','19910101','19920101']
nyear=len(year)
nt=12*nyear

# vars needed for rad calc
tsurf=np.zeros((ncase,nt,ny))
insol=np.zeros((ncase,nt,ny))
temp=np.zeros((ncase,nt,nz,ny))
rh=np.zeros((ncase,nt,nz,ny))
sphum=np.zeros((ncase,nt,nz,ny))
qs=np.zeros((ncase,nt,nz,ny))
stoch_cld=np.zeros((ncase,nt,nz,ny))
stoch_lwp=np.zeros((ncase,nt,nz,ny))
stoch_drop_size=np.zeros((ncase,nt,nz,ny))
stoch_iwp=np.zeros((ncase,nt,nz,ny))
stoch_ice_size=np.zeros((ncase,nt,nz,ny))

# read in data
for k in range(0,ncase): #T400-C400, T400-C800
    if k==0:
        if '8x' in name:
            fname='T400-C100-'
        elif '4x' in name:
            fname='T400-C200-'
        elif '2x' in name:
            fname='T400-C400-'
    elif k==1:
        fname='T400-C800-'
    directory=data_dir+fname+casename+'/history/'

    for i in range(0,nyear):
        filename=directory+'lon_mean_'+year[i]+'_'+time_domain+'.nc'
        ncfile = Dataset(filename,'r')
        
        #some hiram output vars
        olr[k,:]=olr[k,:]+np.mean(ncfile.variables['olr'][:,0,:,0],axis=0)/nyear
        olr_clr[k,:]=olr_clr[k,:]+np.mean(ncfile.variables['olr_clr'][:,0,:,0],axis=0)/nyear
        swup_toa[k,:]=swup_toa[k,:]+np.mean(ncfile.variables['swup_toa'][:,0,:,0],axis=0)/nyear
        swup_toa_clr[k,:]=swup_toa_clr[k,:]+np.mean(ncfile.variables['swup_toa_clr'][:,0,:,0],axis=0)/nyear
        swdn_toa[k,:]=swdn_toa[k,:]+np.mean(ncfile.variables['swdn_toa'][:,0,:,0],axis=0)/nyear
        
        swup_sfc[k,:]=swup_sfc[k,:]+np.mean(ncfile.variables['swup_sfc'][:,0,:,0],axis=0)/nyear
        swdn_sfc[k,:]=swdn_sfc[k,:]+np.mean(ncfile.variables['swdn_sfc'][:,0,:,0],axis=0)/nyear
        lwup_sfc[k,:]=lwup_sfc[k,:]+np.mean(ncfile.variables['lwup_sfc'][:,0,:,0],axis=0)/nyear
        lwdn_sfc[k,:]=lwdn_sfc[k,:]+np.mean(ncfile.variables['lwdn_sfc'][:,0,:,0],axis=0)/nyear

        swup_sfc_clr[k,:]=swup_sfc_clr[k,:]+np.mean(ncfile.variables['swup_sfc_clr'][:,0,:,0],axis=0)/nyear
        swdn_sfc_clr[k,:]=swdn_sfc_clr[k,:]+np.mean(ncfile.variables['swdn_sfc_clr'][:,0,:,0],axis=0)/nyear
        lwup_sfc_clr[k,:]=lwup_sfc_clr[k,:]+np.mean(ncfile.variables['lwup_sfc_clr'][:,0,:,0],axis=0)/nyear
        lwdn_sfc_clr[k,:]=lwdn_sfc_clr[k,:]+np.mean(ncfile.variables['lwdn_sfc_clr'][:,0,:,0],axis=0)/nyear

        precip[k,:]=precip[k,:]+np.mean(ncfile.variables['precip'][:,0,:,0],axis=0)/nyear

        # vars needed for rad calc
        tsurf[k,i*12:(i+1)*12,:]=ncfile.variables['t_surf'][0:12,0,:,0]
        insol[k,i*12:(i+1)*12,:]=ncfile.variables['swdn_toa'][0:12,0,:,0]
        temp[k,i*12:(i+1)*12,:,:]=ncfile.variables['temp'][0:12,:,:,0]
        rh[k,i*12:(i+1)*12,:,:]=0.01*ncfile.variables['rh'][0:12,:,:,0]
        sphum[k,i*12:(i+1)*12,:,:]=ncfile.variables['sphum'][0:12,:,:,0]
        if 'cre1' in casename:
            stoch_cld[k,i*12:(i+1)*12,:,:]=np.maximum(ncfile.variables['stoch_cld_ave'][0:12,:,:,0],0.) #3d
            if 'l0' not in casename:
                stoch_lwp[k,i*12:(i+1)*12,:,:]=np.maximum(ncfile.variables['stoch_drop_conc_ave'][0:12,:,:,0],0.) #3d, g/m3
                var=0.5*ncfile.variables['stoch_drop_size_ave'][0:12,:,:,0] #3d, diameter --> radius
                var[var<2.5]=2.5
                var[var>60.]=60.
                stoch_drop_size[k,i*12:(i+1)*12,:,:]=var
            if 'i0' not in casename:
                stoch_iwp[k,i*12:(i+1)*12,:,:]=np.maximum(ncfile.variables['stoch_ice_conc_ave'][0:12,:,:,0],0.) #3d, g/m3
                var=0.5*ncfile.variables['stoch_ice_size_ave'][0:12,:,:,0] #3d, diameter --> radius
                var[var<13.]=13.
                var[var>130.]=130.
                stoch_drop_size[k,i*12:(i+1)*12,:,:]=var

    for j in range(0,nz):
        stoch_iwp[k,:,j,:]=stoch_iwp[k,:,j,:]*np.log(phalf[j+1]/phalf[j])*Rd*temp[k,:,j,:]/g
        stoch_lwp[k,:,j,:]=stoch_lwp[k,:,j,:]*np.log(phalf[j+1]/phalf[j])*Rd*temp[k,:,j,:]/g
        qs[k,:,j,:]=esat(temp[k,:,j,:])/(100.*p[j])*Rd/Rv
    #sphum[k,:]=rh[k,:]*qs[k,:]
    

# radiative calculation here:

rad_lw_up=np.zeros((6,nz+1,ny))
rad_lw_dn=np.zeros((6,nz+1,ny))
rad_lw_up_clr=np.zeros((6,nz+1,ny))
rad_lw_dn_clr=np.zeros((6,nz+1,ny))

rad_sw_up=np.zeros((6,nz+1,ny))
rad_sw_dn=np.zeros((6,nz+1,ny))
rad_sw_up_clr=np.zeros((6,nz+1,ny))
rad_sw_dn_clr=np.zeros((6,nz+1,ny))

for k in range(0,6):
    if k==0:
        k0=0
        kt=0 #temperature
        kc=0 #cloud
        kw=0 #water vapour
        if '2x' in name:
            absorber=absorber_C400
        elif '4x' in name:
            absorber=absorber_C200
        elif '8x' in name:
            absorber=absorber_C100
    elif k==1:
        k0=1
        kt=1 #temperature
        kc=1 #cloud
        kw=1 #water vapour
        absorber=absorber_C800
    elif k==2: #instant forcing
        k0=0
        kt=0 #temperature
        kc=0 #cloud
        kw=0 #water vapour
        absorber=absorber_C800
    elif k==3: #temperature adjustment
        k0=0
        kt=1 #temperature
        kc=0 #cloud
        kw=0 #water vapour
        absorber=absorber_C800
    elif k==4: #cloud adjustment
        k0=0
        kt=0 #temperature
        kc=1 #cloud
        kw=0 #water vapour
        absorber=absorber_C800
    elif k==5: #water vapour adjustment
        k0=0
        kt=0 #temperature
        kc=0 #cloud
        kw=1 #water vapour
        absorber=absorber_C800
    for l in range(0,nt):
        insolation = field.Field(insol[k0,l,:], domain=sfc)
        ts = field.Field(tsurf[k0,l,:], domain=sfc)        
        temperature = field.Field(temp[kt,l,:,:].T, domain=atm)
        state= {'Tatm': temperature, 'Ts': ts}
        specific_humidity = field.Field(sphum[kw,l,:,:].T, domain=atm)
        
        cldfrac = field.Field(stoch_cld[kc,l,:,:].T, domain=atm)
        clwp = field.Field(stoch_lwp[kc,l,:,:].T, domain=atm)
        ciwp = field.Field(stoch_iwp[kc,l,:,:].T, domain=atm)
        r_liq = field.Field(stoch_drop_size[kc,l,:,:].T, domain=atm)
        r_ice = field.Field(stoch_ice_size[kc,l,:,:].T, domain=atm)        

        #insolation = field.Field(sym_lat(insol[k0,l,:]), domain=sfc)
        #ts = field.Field(sym_lat(tsurf[k0,l,:]), domain=sfc)        
        #temperature = field.Field(sym_lat(temp[kt,l,:,:]).T, domain=atm)
        #state= {'Tatm': temperature, 'Ts': ts}
        #specific_humidity = field.Field(sym_lat(sphum[kw,l,:,:]).T, domain=atm)
        
        #cldfrac = field.Field(sym_lat(stoch_cld[kc,l,:,:]).T, domain=atm)
        #clwp = field.Field(sym_lat(stoch_lwp[kc,l,:,:]).T, domain=atm)
        #ciwp = field.Field(sym_lat(stoch_iwp[kc,l,:,:]).T, domain=atm)
        #r_liq = field.Field(sym_lat(stoch_drop_size[kc,l,:,:]).T, domain=atm)
        #r_ice = field.Field(sym_lat(stoch_ice_size[kc,l,:,:]).T, domain=atm)

        r_liq[r_liq<2.5]=2.5
        r_liq[r_liq>60.]=60.
        r_ice[r_ice<13.]=13.
        r_ice[r_ice>130.]=130.
    
        rad = climlab.radiation.RRTMG(name='Radiation', state=state, 
            specific_humidity=specific_humidity, absorber_vmr=absorber, 
            albedo=albedo, insolation = insolation, cldfrac=cldfrac, 
            r_liq=r_liq, r_ice=r_ice, clwp=clwp, ciwp=ciwp)
        rad.compute_diagnostics()
        rad_lw_up[k,:,:]=rad_lw_up[k,:,:]+rad.LW_flux_up.T/nt
        rad_lw_dn[k,:,:]=rad_lw_dn[k,:,:]+rad.LW_flux_down.T/nt
        rad_lw_up_clr[k,:,:]=rad_lw_up_clr[k,:,:]+rad.LW_flux_up_clr.T/nt
        rad_lw_dn_clr[k,:,:]=rad_lw_dn_clr[k,:,:]+rad.LW_flux_down_clr.T/nt
        rad_sw_up[k,:,:]=rad_sw_up[k,:,:]+rad.SW_flux_up.T/nt
        rad_sw_dn[k,:,:]=rad_sw_dn[k,:,:]+rad.SW_flux_down.T/nt
        rad_sw_up_clr[k,:,:]=rad_sw_up_clr[k,:,:]+rad.SW_flux_up_clr.T/nt
        rad_sw_dn_clr[k,:,:]=rad_sw_dn_clr[k,:,:]+rad.SW_flux_down_clr.T/nt

# instant forcing at toa
flw_inst=rad_lw_up[0,0,:]-rad_lw_up[2,0,:]
flw_inst_clr=rad_lw_up_clr[0,0,:]-rad_lw_up_clr[2,0,:]
fsw_inst=rad_sw_up[0,0,:]-rad_sw_up[2,0,:]
fsw_inst_clr=rad_sw_up_clr[0,0,:]-rad_sw_up_clr[2,0,:]

# instant forcing at surface
flw_srf_inst=-rad_lw_up[0,nz,:]+rad_lw_up[2,nz,:]+rad_lw_dn[0,nz,:]-rad_lw_dn[2,nz,:]
flw_srf_inst_clr=-rad_lw_up_clr[0,nz,:]+rad_lw_up_clr[2,nz,:]+rad_lw_dn_clr[0,nz,:]-rad_lw_dn_clr[2,nz,:]
fsw_srf_inst=-rad_sw_up[0,nz,:]+rad_sw_up[2,nz,:]+rad_sw_dn[0,nz,:]-rad_sw_dn[2,nz,:]
fsw_srf_inst_clr=-rad_sw_up_clr[0,nz,:]+rad_sw_up_clr[2,nz,:]+rad_sw_dn_clr[0,nz,:]-rad_sw_dn_clr[2,nz,:]

# instant forcing at tropopause
flw_itrop=np.zeros(ny)
flw_itrop_clr=np.zeros(ny)
fsw_itrop=np.zeros(ny)
fsw_itrop_clr=np.zeros(ny)
for j in range(0,ny):
    flw_itrop[j]=rad_lw_up[0,ntrop[j],j]-rad_lw_up[2,ntrop[j],j]
    flw_itrop_clr[j]=rad_lw_up_clr[0,ntrop[j],j]-rad_lw_up_clr[2,ntrop[j],j]
    fsw_itrop[j]=rad_sw_up[0,ntrop[j],j]-rad_sw_up[2,ntrop[j],j]
    fsw_itrop_clr[j]=rad_sw_up_clr[0,ntrop[j],j]-rad_sw_up_clr[2,ntrop[j],j]

# co2 forcing after adjustment, rrtm calc
flw_eff=rad_lw_up[0,0,:]-rad_lw_up[1,0,:]
flw_eff_clr=rad_lw_up_clr[0,0,:]-rad_lw_up_clr[1,0,:]
fsw_eff=rad_sw_up[0,0,:]-rad_sw_up[1,0,:]
fsw_eff_clr=rad_sw_up_clr[0,0,:]-rad_sw_up_clr[1,0,:]

# co2 forcing at surface after adjustment, rrtm calc
flw_srf_eff=-rad_lw_up[0,nz,:]+rad_lw_up[1,nz,:]+rad_lw_dn[0,nz,:]-rad_lw_dn[1,nz,:]
flw_srf_eff_clr=-rad_lw_up_clr[0,nz,:]+rad_lw_up_clr[1,nz,:]+rad_lw_dn_clr[0,nz,:]-rad_lw_dn_clr[1,nz,:]
fsw_srf_eff=-rad_sw_up[0,nz,:]+rad_sw_up[1,nz,:]+rad_sw_dn[0,nz,:]-rad_sw_dn[1,nz,:]
fsw_srf_eff_clr=-rad_sw_up_clr[0,nz,:]+rad_sw_up_clr[1,nz,:]+rad_sw_dn_clr[0,nz,:]-rad_sw_dn_clr[1,nz,:]

# co2 forcing after adjustment, hiram calc
flw_hiram=olr[0,:]-olr[1,:]
flw_hiram_clr=olr_clr[0,:]-olr_clr[1,:]
fsw_hiram=swup_toa[0,:]-swup_toa[1,:]
fsw_hiram_clr=swup_toa_clr[0,:]-swup_toa_clr[1,:]

# co2 forcing at surface after adjustment, hiram calc
flw_srf_hiram=-lwup_sfc[0,:]+lwup_sfc[1,:]+lwdn_sfc[0,:]-lwdn_sfc[1,:]
flw_srf_hiram_clr=-lwup_sfc_clr[0,:]+lwup_sfc_clr[1,:]+lwdn_sfc_clr[0,:]-lwdn_sfc_clr[1,:]
fsw_srf_hiram=-swup_sfc[0,:]+swup_sfc[1,:]+swdn_sfc[0,:]-swdn_sfc[1,:]
fsw_srf_hiram_clr=-swup_sfc_clr[0,:]+swup_sfc_clr[1,:]+swdn_sfc_clr[0,:]-swdn_sfc_clr[1,:]

#fast adjustment due to temperature only
flw_temp=rad_lw_up[0,0,:]-rad_lw_up[3,0,:]
flw_temp_clr=rad_lw_up_clr[0,0,:]-rad_lw_up_clr[3,0,:]
fsw_temp=rad_sw_up[0,0,:]-rad_sw_up[3,0,:]
fsw_temp_clr=rad_sw_up_clr[0,0,:]-rad_sw_up_clr[3,0,:]
flw_temp=flw_temp-flw_inst
flw_temp_clr=flw_temp_clr-flw_inst_clr
fsw_temp=fsw_temp-fsw_inst
fsw_temp_clr=fsw_temp_clr-fsw_inst_clr

#fast adjustment due to water vapour only
flw_wv=rad_lw_up[0,0,:]-rad_lw_up[5,0,:]
flw_wv_clr=rad_lw_up_clr[0,0,:]-rad_lw_up_clr[5,0,:]
fsw_wv=rad_sw_up[0,0,:]-rad_sw_up[5,0,:]
fsw_wv_clr=rad_sw_up_clr[0,0,:]-rad_sw_up_clr[5,0,:]
flw_wv=flw_wv-flw_inst
flw_wv_clr=flw_wv_clr-flw_inst_clr
fsw_wv=fsw_wv-fsw_inst
fsw_wv_clr=fsw_wv_clr-fsw_inst_clr

#fast adjustment due to cloud only
flw_cld=rad_lw_up[0,0,:]-rad_lw_up[4,0,:]
flw_cld_clr=rad_lw_up_clr[0,0,:]-rad_lw_up_clr[4,0,:]
fsw_cld=rad_sw_up[0,0,:]-rad_sw_up[4,0,:]
fsw_cld_clr=rad_sw_up_clr[0,0,:]-rad_sw_up_clr[4,0,:]
flw_cld=flw_cld-flw_inst
flw_cld_clr=flw_cld_clr-flw_inst_clr
fsw_cld=fsw_cld-fsw_inst
fsw_cld_clr=fsw_cld_clr-fsw_inst_clr

flw_cld_hiram=flw_hiram-flw_inst-flw_temp-flw_wv
flw_cld_hiram_clr=flw_hiram_clr-flw_inst_clr-flw_temp_clr-flw_wv_clr #some error here
fsw_cld_hiram=fsw_hiram-fsw_inst-fsw_temp-fsw_wv
fsw_cld_hiram_clr=fsw_hiram_clr-fsw_inst_clr-fsw_temp_clr-fsw_wv_clr #some error here


#plt.figure(1)
#plt.title(r'hiram/rrtm olr diff')
#plt.plot(lat,rad_lw_up[0,0,:]-olr[0,:],'-')
#plt.plot(lat,rad_lw_up[1,0,:]-olr[1,:],'-')

#plt.figure(2)
#plt.title(r'hiram/rrtm sw diff')
#plt.plot(lat,rad_sw_up[0,0,:]-swup_toa[0,:],'-')
#plt.plot(lat,rad_sw_up[1,0,:]-swup_toa[1,:],'-')

#plt.figure(3) #okay
#plt.title(r'hiram/rrtm sw diff')
#plt.plot(lat,rad_sw_dn[0,0,:]-swdn_toa[0,:],'--')
#plt.plot(lat,rad_sw_dn[1,0,:]-swdn_toa[1,:],'--')

#plt.figure(4)
#plt.title(r'hiram/rrtm sw sfc')
#plt.plot(lat,swup_sfc[0,:]-rad_sw_up[0,nz,:],'-')
#plt.plot(lat,swdn_sfc[0,:]-rad_sw_dn[0,nz,:],'-')


print('***forcing lw+sw***')
print('eff hiram',gb(flw_hiram)+gb(fsw_hiram))
print('eff rrtm',gb(flw_eff)+gb(fsw_eff))
print('inst toa',gb(flw_inst)+gb(fsw_inst))
print('inst trop',gb(flw_itrop)+gb(fsw_itrop))
print('adj temp',gb(flw_temp)+gb(fsw_temp))

print('***adjustment lw+sw***')
print('total adj',gb(flw_eff+fsw_eff-(flw_inst+fsw_inst)))
print('temp',gb(flw_temp+fsw_temp))
print('cloud rrtm/hiram',gb(flw_cld+fsw_cld),gb(flw_cld_hiram+fsw_cld_hiram))
print('water vapour',gb(flw_wv+fsw_wv))

print('flw_cld rrtm/hiram',gb(flw_cld),gb(flw_cld_hiram))
print('fsw_cld rrtm/hiram',gb(fsw_cld),gb(fsw_cld_hiram))
print('flw_cld_clr rrtm/hiram',gb(flw_cld_clr),gb(flw_cld_hiram_clr))
print('fsw_cld_clr rrtm/hiram',gb(fsw_cld_clr),gb(fsw_cld_hiram_clr))

font_size=20
from matplotlib import rc
rc('text', usetex=False)  #very slow
#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('font', **{'size' : 20})

width=11.
height=6.
left=0.1
bottom=0.15
right=0.95
top=.9

figs_directory='../../../papers/2018-ECS/figs-201903/'
txt=casename.replace(".","")



plt.figure()
fig = matplotlib.pyplot.gcf()
ax = plt.gca()
fig.set_size_inches(width, height)
plt.plot(lat,sym_lat(flw_inst+fsw_inst),label=r'IRF',color=def_colors[0])
#plt.plot(lat,flw_itrop+fsw_itrop,label=r'inst trop rrtm',color=def_colors[2])
#plt.plot(lat,flw_hiram+fsw_hiram,label=r'ERF',color=def_colors[0])
plt.plot(lat,sym_lat(flw_eff+fsw_eff),label=r'ERF',color=def_colors[1])
plt.plot(lat,sym_lat(flw_temp+fsw_temp),label=r'temp. adj.',color=def_colors[3])
plt.plot(lat,sym_lat(flw_wv+fsw_wv),label=r'water vapour adj.',color='k',linestyle='--')
if 'cre1' in casename:
    plt.plot(lat,sym_lat(flw_cld+fsw_cld),label=r'cld. adj.',color=def_colors[2])
#plt.plot(lat,flw_cld_hiram+fsw_cld_hiram,label=r'cld adj hiram',color=def_colors[4])
#plt.plot(lat,0.5*Lv*(precip[0,:]-precip[1,:]),color=def_colors[3],label=r'$-$0.5*Lv*dprecip')
#plt.plot([-90.,90.],[0.,0.],'k:')
plt.legend()
plt.xlabel(r'Latitude')
plt.xlim(-90.,90.)
plt.xticks(np.arange(-90., 120., step=30))
plt.ylabel(r'(W$\,$m$^{-2}$)')
plt.title(r'TOA radiative flux change',fontsize=font_size)
plt.subplots_adjust(left=left,bottom=bottom,right=right,top=top)
#filename=figs_directory+'toa-radflx-lat-'+name+'.pdf'
#plt.savefig(filename, bbox_inches=0)





plt.figure()
fig = matplotlib.pyplot.gcf()
ax = plt.gca()
fig.set_size_inches(6, 6.)
plt.plot(lat,sym_lat(flw_cld),label=r'adj. lw',color=def_colors[0])
plt.plot(lat,sym_lat(fsw_cld),label=r'adj. sw',color=def_colors[1])
plt.plot(lat,sym_lat(flw_inst-flw_inst_clr),
    label=r'mask. lw',color=def_colors[2],linestyle='-')
plt.plot(lat,sym_lat(fsw_inst-fsw_inst_clr),
    label=r'mask. sw',color=def_colors[3],linestyle='-')
#plt.plot(lat,sym_lat(flw_cld+fsw_cld+flw_inst-flw_inst_clr+fsw_inst-fsw_inst_clr),
#    label=r'cld. total',color=def_colors[2])
#plt.plot(lat,sym_lat(flw_cld),label=r'cld. adj. lw',color=def_colors[0])
#plt.plot(lat,sym_lat(fsw_cld),label=r'cld. adj. sw',color=def_colors[1])
#plt.plot([-90.,90.],[0.,0.],'k:')
plt.legend(ncol=2,loc='lower center')
plt.xlabel(r'Latitude')
plt.xlim(0.,90.)
plt.xticks(np.arange(0., 120., step=30))
plt.ylim(-5.,5.)
plt.ylabel(r'(W$\,$m$^{-2}$)')
plt.title(r'Cloud Radiative Effect at TOA',fontsize=font_size)
#plt.subplots_adjust(left=left,bottom=bottom,right=right,top=top)
plt.subplots_adjust(left=0.15,bottom=0.13,right=0.95,top=0.93)
#filename=figs_directory+'toa-radflxcld-lat-'+name+'.pdf'
#plt.savefig(filename, bbox_inches=0)


#plt.figure()
#fig = matplotlib.pyplot.gcf()
#ax = plt.gca()
#fig.set_size_inches(6, 6.)
#plt.plot(lat,sym_lat(flw_cld),color=def_colors[0],linestyle='--')
#plt.plot(lat,sym_lat(fsw_cld),color=def_colors[1],linestyle='--')
#plt.plot(lat,sym_lat(flw_inst-flw_inst_clr),color=def_colors[2],linestyle='--')
#plt.plot(lat,sym_lat(fsw_inst-fsw_inst_clr),color=def_colors[3],linestyle='--')
##plt.legend(ncol=2,loc='lower center')
#plt.xlabel(r'Latitude')
#plt.xlim(0.,90.)
#plt.xticks(np.arange(0., 120., step=30))
#plt.ylim(-6.,7.5)
#plt.ylabel(r'(W$\,$m$^{-2}$)')
#plt.title(r'Cloud Radiative Effect at TOA',fontsize=font_size)
##plt.subplots_adjust(left=left,bottom=bottom,right=right,top=top)
#plt.subplots_adjust(left=0.15,bottom=0.13,right=0.95,top=0.93)
#filename=figs_directory+'toa-radflxcld-lat.pdf'
#plt.savefig(filename, bbox_inches=0)


## save data to file
directory='radforc-out/'


## save latitudinal data to file
root_grp = Dataset(directory+name+'_lat.nc', 'w', format='NETCDF4')

# y-dimension
root_grp.createDimension('lat', ny)

#forcing
var=root_grp.createVariable('flw_inst', 'f8', ('lat',))
var[:]=flw_inst
var=root_grp.createVariable('flw_inst_clr', 'f8', ('lat',))
var[:]=flw_inst_clr

var=root_grp.createVariable('fsw_inst', 'f8', ('lat',))
var[:]=fsw_inst
var=root_grp.createVariable('fsw_inst_clr', 'f8', ('lat',))
var[:]=fsw_inst_clr

var=root_grp.createVariable('flw_srf_inst', 'f8', ('lat',))
var[:]=flw_srf_inst
var=root_grp.createVariable('flw_srf_inst_clr', 'f8', ('lat',))
var[:]=flw_srf_inst_clr

var=root_grp.createVariable('fsw_srf_inst', 'f8', ('lat',))
var[:]=fsw_srf_inst
var=root_grp.createVariable('fsw_srf_inst_clr', 'f8', ('lat',))
var[:]=fsw_srf_inst_clr

var=root_grp.createVariable('flw_itrop', 'f8', ('lat',))
var[:]=flw_itrop
var=root_grp.createVariable('flw_itrop_clr', 'f8', ('lat',))
var[:]=flw_itrop_clr

var=root_grp.createVariable('fsw_itrop', 'f8', ('lat',))
var[:]=fsw_itrop
var=root_grp.createVariable('fsw_itrop_clr', 'f8', ('lat',))
var[:]=fsw_itrop_clr

var=root_grp.createVariable('flw_temp', 'f8', ('lat',))
var[:]=flw_temp
var=root_grp.createVariable('flw_temp_clr', 'f8', ('lat',))
var[:]=flw_temp_clr

var=root_grp.createVariable('fsw_temp', 'f8', ('lat',))
var[:]=fsw_temp
var=root_grp.createVariable('fsw_temp_clr', 'f8', ('lat',))
var[:]=fsw_temp_clr

var=root_grp.createVariable('flw_cld', 'f8', ('lat',))
var[:]=flw_cld
var=root_grp.createVariable('flw_cld_clr', 'f8', ('lat',))
var[:]=flw_cld_clr

var=root_grp.createVariable('fsw_cld', 'f8', ('lat',))
var[:]=fsw_cld
var=root_grp.createVariable('fsw_cld_clr', 'f8', ('lat',))
var[:]=fsw_cld_clr

var=root_grp.createVariable('flw_wv', 'f8', ('lat',))
var[:]=flw_wv
var=root_grp.createVariable('flw_wv_clr', 'f8', ('lat',))
var[:]=flw_wv_clr

var=root_grp.createVariable('fsw_wv', 'f8', ('lat',))
var[:]=fsw_wv
var=root_grp.createVariable('fsw_wv_clr', 'f8', ('lat',))
var[:]=fsw_wv_clr

var=root_grp.createVariable('flw_eff', 'f8', ('lat',))
var[:]=flw_eff
var=root_grp.createVariable('flw_eff_clr', 'f8', ('lat',))
var[:]=flw_eff_clr

var=root_grp.createVariable('fsw_eff', 'f8', ('lat',))
var[:]=fsw_eff
var=root_grp.createVariable('fsw_eff_clr', 'f8', ('lat',))
var[:]=fsw_eff_clr

var=root_grp.createVariable('flw_srf_eff', 'f8', ('lat',))
var[:]=flw_srf_eff
var=root_grp.createVariable('flw_srf_eff_clr', 'f8', ('lat',))
var[:]=flw_srf_eff_clr

var=root_grp.createVariable('fsw_srf_eff', 'f8', ('lat',))
var[:]=fsw_srf_eff
var=root_grp.createVariable('fsw_srf_eff_clr', 'f8', ('lat',))
var[:]=fsw_srf_eff_clr

var=root_grp.createVariable('flw_hiram', 'f8', ('lat',))
var[:]=flw_hiram
var=root_grp.createVariable('flw_hiram_clr', 'f8', ('lat',))
var[:]=flw_hiram_clr

var=root_grp.createVariable('fsw_hiram', 'f8', ('lat',))
var[:]=gb(fsw_hiram)
var=root_grp.createVariable('fsw_hiram_clr', 'f8', ('lat',))
var[:]=gb(fsw_hiram_clr)

var=root_grp.createVariable('flw_srf_hiram', 'f8', ('lat',))
var[:]=flw_srf_hiram
var=root_grp.createVariable('flw_srf_hiram_clr', 'f8', ('lat',))
var[:]=flw_srf_hiram_clr

var=root_grp.createVariable('fsw_srf_hiram', 'f8', ('lat',))
var[:]=gb(fsw_srf_hiram)
var=root_grp.createVariable('fsw_srf_hiram_clr', 'f8', ('lat',))
var[:]=gb(fsw_srf_hiram_clr)

var=root_grp.createVariable('flw_cld_hiram', 'f8', ('lat',))
var[:]=flw_cld_hiram
var=root_grp.createVariable('flw_cld_hiram_clr', 'f8', ('lat',))
var[:]=flw_cld_hiram_clr

var=root_grp.createVariable('fsw_cld_hiram', 'f8', ('lat',))
var[:]=fsw_cld_hiram
var=root_grp.createVariable('fsw_cld_hiram_clr', 'f8', ('lat',))
var[:]=fsw_cld_hiram_clr

root_grp.close()



### save global data to file
root_grp = Dataset(directory+name+'.nc', 'w', format='NETCDF4')

#forcing
var=root_grp.createVariable('flw_inst', 'f8')
var[:]=gb(flw_inst)
var=root_grp.createVariable('flw_inst_clr', 'f8')
var[:]=gb(flw_inst_clr)

var=root_grp.createVariable('fsw_inst', 'f8')
var[:]=gb(fsw_inst)
var=root_grp.createVariable('fsw_inst_clr', 'f8')
var[:]=gb(fsw_inst_clr)

var=root_grp.createVariable('flw_srf_inst', 'f8')
var[:]=gb(flw_srf_inst)
var=root_grp.createVariable('flw_srf_inst_clr', 'f8')
var[:]=gb(flw_srf_inst_clr)

var=root_grp.createVariable('fsw_srf_inst', 'f8')
var[:]=gb(fsw_srf_inst)
var=root_grp.createVariable('fsw_srf_inst_clr', 'f8')
var[:]=gb(fsw_srf_inst_clr)

var=root_grp.createVariable('flw_itrop', 'f8')
var[:]=gb(flw_itrop)
var=root_grp.createVariable('flw_itrop_clr', 'f8')
var[:]=gb(flw_itrop_clr)

var=root_grp.createVariable('fsw_itrop', 'f8')
var[:]=gb(fsw_itrop)
var=root_grp.createVariable('fsw_itrop_clr', 'f8')
var[:]=gb(fsw_itrop_clr)

var=root_grp.createVariable('flw_temp', 'f8')
var[:]=gb(flw_temp)
var=root_grp.createVariable('flw_temp_clr', 'f8')
var[:]=gb(flw_temp_clr)

var=root_grp.createVariable('fsw_temp', 'f8')
var[:]=gb(fsw_temp)
var=root_grp.createVariable('fsw_temp_clr', 'f8')
var[:]=gb(fsw_temp_clr)

var=root_grp.createVariable('flw_cld', 'f8')
var[:]=gb(flw_cld)
var=root_grp.createVariable('flw_cld_clr', 'f8')
var[:]=gb(flw_cld_clr)

var=root_grp.createVariable('fsw_cld', 'f8')
var[:]=gb(fsw_cld)
var=root_grp.createVariable('fsw_cld_clr', 'f8')
var[:]=gb(fsw_cld_clr)

var=root_grp.createVariable('flw_wv', 'f8')
var[:]=gb(flw_wv)
var=root_grp.createVariable('flw_wv_clr', 'f8')
var[:]=gb(flw_wv_clr)

var=root_grp.createVariable('fsw_wv', 'f8')
var[:]=gb(fsw_wv)
var=root_grp.createVariable('fsw_wv_clr', 'f8')
var[:]=gb(fsw_wv_clr)

var=root_grp.createVariable('flw_eff', 'f8')
var[:]=gb(flw_eff)
var=root_grp.createVariable('flw_eff_clr', 'f8')
var[:]=gb(flw_eff_clr)

var=root_grp.createVariable('fsw_eff', 'f8')
var[:]=gb(fsw_eff)
var=root_grp.createVariable('fsw_eff_clr', 'f8')
var[:]=gb(fsw_eff_clr)

var=root_grp.createVariable('flw_srf_eff', 'f8')
var[:]=gb(flw_srf_eff)
var=root_grp.createVariable('flw_srf_eff_clr', 'f8')
var[:]=gb(flw_srf_eff_clr)

var=root_grp.createVariable('fsw_srf_eff', 'f8')
var[:]=gb(fsw_srf_eff)
var=root_grp.createVariable('fsw_srf_eff_clr', 'f8')
var[:]=gb(fsw_srf_eff_clr)

var=root_grp.createVariable('flw_hiram', 'f8')
var[:]=gb(flw_hiram)
var=root_grp.createVariable('flw_hiram_clr', 'f8')
var[:]=gb(flw_hiram_clr)

var=root_grp.createVariable('fsw_hiram', 'f8')
var[:]=gb(fsw_hiram)
var=root_grp.createVariable('fsw_hiram_clr', 'f8')
var[:]=gb(fsw_hiram_clr)

var=root_grp.createVariable('flw_srf_hiram', 'f8')
var[:]=gb(flw_srf_hiram)
var=root_grp.createVariable('flw_srf_hiram_clr', 'f8')
var[:]=gb(flw_srf_hiram_clr)

var=root_grp.createVariable('fsw_srf_hiram', 'f8')
var[:]=gb(fsw_srf_hiram)
var=root_grp.createVariable('fsw_srf_hiram_clr', 'f8')
var[:]=gb(fsw_srf_hiram_clr)

var=root_grp.createVariable('flw_cld_hiram', 'f8')
var[:]=gb(flw_cld_hiram)
var=root_grp.createVariable('flw_cld_hiram_clr', 'f8')
var[:]=gb(flw_cld_hiram_clr)

var=root_grp.createVariable('fsw_cld_hiram', 'f8')
var[:]=gb(fsw_cld_hiram)
var=root_grp.createVariable('fsw_cld_hiram_clr', 'f8')
var[:]=gb(fsw_cld_hiram_clr)

root_grp.close()
