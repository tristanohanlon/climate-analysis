import numpy as np#code to analyse the complete cloud fraction dataset
"""Created on the 26th of March"""
def cloud_fraction(ff):
    array1=[]
    array11=[]
    array2=[]
    array22=[]
    array3=[]
    array33=[]
    ii=0
    undefined=[]
    if len(ff)>0:
        for i in ff:
            try:
                print ii
                #print i
                from pyhdf import SD
    
                import pylab as py
                import numpy, scipy
                from scipy import interpolate
                import numpy as numpy
                import numpy as np
    
                filename=i
                f = SD.SD(filename)
                subdataset_name='Latitude'
                sds =f.select(subdataset_name)
                lat=sds.get()
                subdataset_name='Longitude'
                sds =f.select(subdataset_name)
                lon=sds.get()
                subdataset_name ='Cloud_Top_Height_Nadir'
                sds = f.select(subdataset_name)
                height = sds.get()
                subdataset_name='Cloud_Fraction_Nadir'
                sds=f.select(subdataset_name)
                cloud_fraction=sds.get()
                z1=np.where(cloud_fraction<>sds.attributes()['_FillValue'])
                cloud_fraction=cloud_fraction[z1].ravel()*sds.attributes()['scale_factor']
                lat=lat[z1].ravel()
                height=height[z1].ravel()
                lon=lon[z1].ravel()
                y=numpy.where(height<2500) and numpy.where(height>200)
                array1=array1+cloud_fraction[y].tolist()
                array11=array11+cloud_fraction.ravel().tolist()
                array2=array2+(np.round(lat[y]*2)/2).tolist()
                array22=array22+(np.round(lat.ravel()*2)/2).tolist()
                array3=array3+(np.around(lon[y]*2)/2).tolist()
                array33=array33+(np.around(lon.ravel()*2)/2).tolist()
                ii=ii+1
            except:
                print i
                undefined.append(i)
                
            #instead you need to invoke a condition that simply considers the clouds that are not low clouds in a grid point and consider their effect on the average.
        cf=np.array(array1)
        cf_high=np.array(array11)
        lat_high=np.array(array22)
        lon_high=np.array(array33)
        lat1=np.array(array2)
        lon1=np.array(array3)
        lat_bins=np.arange(-30,5.5,0.5)[::-1]
        lon_bins=np.arange(-130,-74.5,0.5)
        lon2=[lon1[np.where(lat1==k)] for k in lat_bins]
        cf1=[cf[np.where(lat1==k)] for k in lat_bins]
        cf22=[cf_high[np.where(lat_high==k)] for k in lat_bins]
        lon1_high=[lon_high[np.where(lat_high==k)] for k in lat_bins]

        std=np.ones([71,111])*np.nan
        grid=np.ones([71,111])*np.nan

        counter=0

        for i in lon2:
            if len(lon2)==len(cf1):
            #if counter<71:
                array_i=lon2[counter]
                cf2=cf1[counter]
                cf_all=cf22[counter]
                lon_high_count=lon1_high[counter]
                #lon_bins=np.arange(-130,-74.5,0.5)
                
                if len(array_i[array_i<-75])>0 and len(array_i[array_i>-130]):
                    for j in range(len(lon_bins)):

                        grid[counter,j]=np.nansum(cf2[np.where(array_i==lon_bins[j])])/len(cf_all[np.where(lon_high_count==lon_bins[j])][cf_all[np.where(lon_high_count==lon_bins[j])]>-0.01])
                        std[counter,j]=np.nanstd(cf2[np.where(i==lon_bins[j])])
                    counter=counter+1

                else:
                    counter=counter+1
                    


        return grid,std,undefined

    else:
        return [],[],[]
def convert(year,day):
    if day<100 and day>=10:
        return str(int(year))+'0'+str(int(day))
    elif day<10:
        return str(int(year))+'00'+str(int(day))
    else:
        return str(int(year))+str(int(day))