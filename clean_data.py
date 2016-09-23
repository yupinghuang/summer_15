# filter out data points in the time series files according to the following criteria
# analog data (scan file names that start with 'a')
# data that are uncalibrated (psredit -c "scale" show units other than Jansky)
# data with a noise baseline that is not consistent with the radiometer equation
import glob
import subprocess
from edit_time_series import time_series_file
import numpy as np
import matplotlib.pyplot as plt
import math
import psrchive

# get the names of the pulsars and save to file because it will be very useful
def incorrect_noise_scale(fn):
    a=psrchive.Archive_load(fn)
    tInt=a.integration_length()
    #rmsOff=a.rms_baseline()
    rmsOff=subprocess.check_output('psrstat -c "off:rms" '+fn+' -Q',shell=True).split()[1]
    product=float(rmsOff)*math.sqrt(tInt) #this quantity should be the same for all 1024chan pulsars
    prodList.append(product)
    print fn+' '+str(tInt)+' '+str(rmsOff)+' '+str(product)


def uncalibrated(fn):
    a=psrchive.Archive_load(fn)
    unit=a.get_scale()
    if unit=='Jansky':
        return False
    else:
        return True

class outliers_id:
    def __init__(self,fn):
        #plt.clf()
        #plt.ion()
        fluxList=[]
        mjdList=[]
        productList=[]
        self.scanFileList=[]
        self.t=time_series_file(fn)
        for line in self.t.fileBuffer:
            if line[0]!='#':
                l=line.split()
                scanFile=l[6]
                mjd=float(l[0])#used to look up and delete entries in the time_series_file interface
                flux=float(l[3])
                #scanFileName=scanFile.split('/')[-1]
                tInt=float(l[1])
                #rmsOff=a.rms_baseline()
                rmsOff=subprocess.check_output('psrstat -c "off:rms" '+scanFile+' -Q',shell=True).split()[1]
                product=float(rmsOff)*math.sqrt(tInt)
                self.scanFileList.append(scanFile)
                fluxList.append(flux)
                mjdList.append(mjd)
                productList.append(product)
        self.data=np.array([mjdList,fluxList,productList])

    def median_absolute_deviation(self,array):
        med=np.median(array)
        return (med,np.median(np.abs(array-med)))

    def identify_outlier(self):
        # Metric: if the median absolute deviation of a point's product(which is the rms:off*sqrt(length)) is more than 10 times the median median absolute deviation then the point will be removed.
        # or if the product is more than 6 times the median median absolute deviation, with the flux being more than 4 times the median median absolute deviation, the point will be removed as well.
        (fluxMed,fluxMad)=self.median_absolute_deviation(self.data[1])
        (productMed,productMad)=self.median_absolute_deviation(self.data[2])
        fluxDev=abs(self.data[1]-fluxMed)/fluxMad
        productDev=abs(self.data[2]-productMed)/productMad
        
        for i in range(0,productDev.shape[0]):
            if productDev[i]>6:
                if productDev[i]>10 or fluxDev[i]>4:
                    dummy=self.t.get_index_from_mjd(self.data[0][i]) 
                    self.t.comment_out()
                    self.write_log(self.scanFileList[i]+' '+str(productDev[i])+' '+str(fluxDev[i]))
        self.t.unload()
        self.write_log('#file '+self.t.filename+' saved')
    def write_log(self,info):
        with open('outlierCleaning.log','a') as logFile:
            logFile.write(info+'\n')
def main():
    timeSeriesList=glob.glob('/DATA/RUBICON_2/joh414/yuping/test/*/*.time_series')
    count=len(timeSeriesList)
    for i in range(0,count):
        try:
            o=outliers_id(timeSeriesList[i])
            o.identify_outlier()
        except ValueError:
            print 'ValueError '+str(o.data.shape)

    
    
