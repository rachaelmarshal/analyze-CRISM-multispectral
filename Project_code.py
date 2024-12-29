#!/usr/bin/env python
# coding: utf-8

# # Programming Solutions Manual
# 
# This manual will explain the code produced to manipulate summary product data, in specific do the following operations: 
# * Read summary products
# * Mask Outliers
# * Display Maps of the Summary Products
# * Carry out Principal Component Analysis 
# * Calculate percentage contribution of each variable to each principal component

# Step 1 : The code snippet below will be utilised to read, mask outliers and save summary products as georeferenced images/maps. 

# In[ ]:


#Import all necessary packages
#import skimage
#import skimage.transform  #use this if you would want to scale down your dataset
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import numpy as np
import gdal
import osr
import pandas as pd
import sys

#set directory and specify file name
directory= r'/home/s6039723/Documents/AcidaliaMosaics/Subsets/subsets' #specify directory
file= '/band' #specify file name 

#outliers defined using the code from Appendix
outliers_viviano = [[0.11655,0.25224],
                    [-0.01910,0.05693],
                    [-0.01200,0.23752],
                    [-0.027382,0.03615],
                    [-0.03663,0.02860],
                    [-0.07869,-0.00267],
                    [-0.03558,0.06528],
                    [-0.01650,0.06663],
                    [-0.03741,0.02718],
                    [-0.00952,0.03849],
                    [-0.01029,0.04600],
                    [-0.19906,0.23912],
                    [-0.10895,0.33268]]
band_nr_viviano = [4,5,14,16,17,18,20,26,32,33,42,50,51] #selected bands after going through all summary product maps
band_names_viviano=['SH600_2','SH770','OLINDEX3','BD1300',
                     'LCPINDEX2','HCPINDEX2','ISLOPE1',
                     'BD1900_2','MIN2200','BD2210_2','SINDEX2','BD3400_2','CINDEX2']
MaskList1=list()

#for loop for reading summary products, masking outliers, saving them as maps
for k in range (60):
            l=k+1
            if l in [1,2,3,6,7,8,9,10,11,12,13,15,19,21,22,23,24,25,27,28,29,30,31,34,35,36,37,38,39,
                     40,41,43,44,45,46,47,48,49,52,53,54,55,56,57,58,59,60]:
                continue
            else:
                fp = directory+file+str(l)+'.tiff'    
                ds = gdal.Open(fp)
                data = ds.ReadAsArray()
                outlier_min=outliers_viviano[band_nr_viviano.index(l)][0]
                outlier_max=outliers_viviano[band_nr_viviano.index(l)][1]
                product=band_names_viviano[band_nr_viviano.index(l)]
                mask_band= np.logical_or(
                                data == 65535,
                                np.logical_or(
                                data < outlier_min,
                                data > outlier_max))
                maskeddata1=np.ma.array(data,mask=mask_band)
                x=np.ma.filled(maskeddata1,0)
                format = "GTiff"
                driver = gdal.GetDriverByName( format )
                dst_ds1 = driver.Create(str(product)+'tr.tiff', 30047,14034, 1, gdal.GDT_Float64)
                dst_ds1.SetGeoTransform( [ -65.06, 0.004, 0, 67.67, 0, -0.0039])   
                srs = osr.SpatialReference()                                        
                f=open(r'/home/s6039723/Mars_2000.prj')     #locate the projection file                                        
                inproj = f.read()
                srs.ImportFromWkt(inproj)
                dst_ds1.SetProjection( srs.ExportToWkt() )           
                dst_ds1.GetRasterBand(1).WriteArray(x)                   
                dst_ds1.FlushCache()    
                maskeddata1=np.ma.array(data,mask=mask_band)
                MaskList1.append(maskeddata1)
                if k is 3:                                        
                    mask_total=mask_band
                else:
                    mask_total= np.logical_or(mask_band,mask_total) #total summed up mask


# Step 2: Apply the summed up mask to the data

# In[ ]:


FinalMaskList=list() 
for i in range(13):
    maskeddata2=np.ma.array(MaskList1[i],mask=mask_total)   
    FinalMaskList.append(maskeddata2)
#denote all the masked values as NaNs    
NoOut=list()     
for m in range(13):
    noout=np.ma.filled(FinalMaskList[m],fill_value=np.nan)
    NoOut.append(noout)


# Step 3: Empty the list of arrays into a 2D array and exclude the NaN values from the analysis

# In[ ]:


RowArray=list()
ColArray=list()
NonNaNArray=list()
for j in range (13):
    values=NoOut[j]
    array_nonNaN=values[np.logical_not(np.isnan(values))]
    [row,col]=np.where(np.logical_not(np.isnan(values)))
    NonNaNArray.append(array_nonNaN)
    RowArray.append(row)
    ColArray.append(col)


# Step 4: Standardise the Data 

# In[ ]:


MB_matrix=np.zeros((NonNaNArray[0].size,13))
for n in range(13):
    MB_array=NonNaNArray[n]
    MB_arrayStd=(MB_array-MB_array.mean())/MB_array.std()
    MB_matrix[:,n]=MB_arrayStd  
#fp2=np.memmap('MB.npy',dtype='float64',shape=(14034,30047),mode='w+') - creating a memmap array
#fp2[:,:]=MB_matrixw[:,:]    - copying values 
#del fp2 
#fp2=np.memmap('MB.npy',dtype='float64',shape=(14034,30047))
#sys.getsizeof(fp2) # to check the size of the memmap object 


# Step 5: Carry out PCA

# In[ ]:


import sklearn
from sklearn import decomposition
from sklearn.decomposition import PCA

pca=PCA(n_components=5) #change according to required components
pca.fit(MB_matrix)     #pca.fit(fp2)
x_pca=pca.transform(MB_matrix) #x_pca=pca.transform(fp2)
array_zeros=np.zeros((14034,30047))
x_comp=pca.components_
plt.plot(np.cumsum(pca.explained_variance_ratio_))
plt.xlabel('No. of Components')
plt.ylabel('Cumulative explained variance')


# Step 6: Save PCA components as geotiff images/maps

# In[ ]:


arrayfinal=np.zeros((14034,30047,5))
for i in range(5):
    array_zeros[RowNew[i],ColNew[i]]=x_pca[:,i]
    arrayfinal[:,:,i]=array_zeros
    format = "GTiff"
    driver = gdal.GetDriverByName( format )

    dst_ds1 = driver.Create('pc'+str(i+0)+'.tiff', 30047,14034, 1, gdal.GDT_Float64)

 

    dst_ds1.SetGeoTransform( [ -65.06, 0.004, 0, 67.67, 0, -0.0039])   
    srs = osr.SpatialReference()                                       
    f=open(r'/home/s6039723/Mars_2000.prj')                                              
    inproj = f.read()
    srs.ImportFromWkt(inproj) 
    dst_ds1.SetProjection( srs.ExportToWkt() )          
    dst_ds1.GetRasterBand(1).WriteArray(arrayfinal[:,:,i])                   
    dst_ds1.FlushCache()    


# Step 7 : Compute variance captured per variable per component (adapted from R implementation) 

# In[ ]:


def varcapr(x,y):
    contrib=x*100/y
    return contrib
                  

comp2=np.transpose(x_comp)
sdev=np.std(comp2, axis=0)
compsdev=np.std(comp2,axis=0)
varcoord=compon2*compsdev
varcos2=varcoord**2
compcos2=np.sum(varcos2,axis=0)
varcontr=np.transpose(np.apply_along_axis(varcapr,1,varcos2,compcos2))
varcontr2=np.transpose(varcontr)

column_names=['PC1','PC2','PC3','PC4','PC5']
row_names=['SH600_2','SH770','OLINDEX3','BD1300',
                     'LCPINDEX2','HCPINDEX2','ISLOPE1',
                     'BD1900_2','MIN2200','BD2210_2','SINDEX2','BD3400_2','CINDEX2']
df=pd.DataFrame(np.transpose(x_comp),columns=column_names)
df.index=row_names
df.insert(0,'Product',value=df.index)
df.plot(kind='bar', stacked=True)
plt.legend(loc='left', bbox_to_anchor=(2.0, 0.5))

# Heat map to visualise loadings 
ax=sns.heatmap(pca.components_,cmap='YlGnBu',xticklabels=row_names,yticklabels=column_names,cbar_kws={"orientation":"horizontal"})
ax.set_aspect("equal")                     


# Appendix: Calculate Outliers 

# In[ ]:


import numpy as np
import gdal
import os
import osr
import scipy
import scipy.ndimage
import pandas as pd

directory= r'/home/s6039723/Documents/AcidaliaMosaics/Subsets/subsets'
file= '/band'
OutlierMin=list()
OutlierMax=list()
for i in range(60):
        fp = directory+file+str(i+1)+'.tiff'    
        ds = gdal.Open(fp)
        data = ds.ReadAsArray()
        if i==5:
            continue
        else:
            data[np.isnan(data)]=65535   
                
            masked_array = np.ma.MaskedArray(
                                data, np.logical_or(data == 0.0, data == 65535)
                                )
                
            Q1 = np.percentile(masked_array.compressed(),25)
            Q3 = np.percentile(masked_array.compressed(),75)
            IQD = 1.5*(Q3-Q1)
                
            outlier_min1 = Q1-IQD
            outlier_max1 = Q3+IQD
            OutlierMin.append(outlier_min1)
            OutlierMax.append(outlier_max1)
df=pd.DataFrame()
df['Min']=OutlierMin
df['Max']=OutlierMax
df.to_csv('Outliers.csv')

