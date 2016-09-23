import subprocess
import psrchive
import math
import os
import glob
import numpy as np
from scipy import interpolate
from scipy import fftpack
import matplotlib.pyplot as plt
from matplotlib import gridspec

class time_series:
    # approximate central freq corresponding to each wavelength
    stdCmToFreq={10:3098,20:1370,50:651}

    def __init__(self,directory,cm,nbin,stdProfDir,overwrite=False):
        self.directory=directory
        # cm is the wavelength in cm, we have 10 20 and 50
        self.cm=cm
        self.nbin=nbin
        self.stdProfDir=stdProfDir
        self.filteredFileNames=[]
        self.sortedData=[]
        # name will be set when psredit is run in the function freqNnbin
        self.name=''
        self.timeSeriesFileName=''
        #search for existing time series file in the current directory. If found, just load it
        potentialTimeSeriesFile=glob.glob(directory+'/'+'*.time_series')
        if (not overwrite) and potentialTimeSeriesFile!=[]:
            self.load(directory,potentialTimeSeriesFile[0])
        else:
            self.freqNnbin()
            self.group_psrflux()

    def load(self,directory,timeSeriesFileName):
        #load time series file, for debugging and testing mostly
        self.directory=directory
        with open(timeSeriesFileName,'r') as f:
            self.name=f.readline().split()[1]
            for line in f:
                if line[0]!='#':
                    self.sortedData.append(line.rstrip('\n').split())
        print 'found time series for '+self.name 
        #tweak timeSeriesFileName such that it does not include the directory for consistency

        self.timeSeriesFileName=timeSeriesFileName.split('/')[-1]


    def freqNnbin(self):
        # freq and nbin needs to be integers
        # list of all *.rf.dzTFp files in the directory
        filenames=glob.glob(self.directory+'/*.rf.dzTFp')
        stdFreq=time_series.stdCmToFreq[self.cm]
        self.name=psrchive.Archive_load(filenames[0]).get_source()
        for filename in filenames:
            archive=psrchive.Archive_load(filename)
            fileFreq=archive.get_centre_frequency()
            fileNbin=archive.get_nbin()
            #Check if the file is calibrated, 1024 bin, and at the right freq
            if abs(fileFreq-stdFreq)<300 and fileNbin==self.nbin and archive.get_poln_calibrated():
                self.filteredFileNames.append(filename)

    def group_psrflux(self):
        self.timeSeriesFileName=self.name+'.'+str(self.cm)+'.time_series'
        # return a list of the data points in tuples
        unsortedData=[]

        for fileName in self.filteredFileNames:
            stdFileName=self.stdProfDir+'/'+self.name+'_'+str(self.cm)+'cm_paas.std'
            try:
                print subprocess.check_output('psrflux '+fileName+' -s '+stdFileName,shell=True)
            except subprocess.CalledProcessError,e:
                with open('timeSeriesErrorLog','a') as errorLog:
                    z=e
                    errorLog.write(str(z)+'\n')
            #integration time from vap. Simon note that psrflux does not give the correct integration time
#            vapOut=subprocess.check_output('vap -nc "length" '+fileName,shell=True).split()
            archive=psrchive.Archive_load(fileName)
            time=archive.integration_length()
            fileNameWithoutDir=archive.get_filename()
            
            try:
                with open(fileName+'.ds','r') as f:
                    lines=f.readlines()
                    mjd=float(lines[4][8:])
                    dataLineList=lines[7].split()
                    freq=float(dataLineList[3])
                    flux=float(dataLineList[4])
                    flux_err=float(dataLineList[5])
                    toa_err=subprocess.check_output('pat -a"'+stdFileName+'" '+fileName,shell=True).split()[4]
                    if time!=0.0:
                        scaledToa_err=float(toa_err)/math.sqrt(time)
                        unsortedData.append((mjd,time,freq,flux,flux_err,scaledToa_err,fileNameWithoutDir))
            except IOError,r:
                with open('SfErrorLog','a') as errorLog:
                    errorLog.write(str(r)+'\n')
                break
        self.sortedData=sorted(unsortedData, key=lambda dat:dat[0])

        with open(self.directory+'/'+self.timeSeriesFileName,'w') as of:
            of.write('#Name '+self.name+' '+str(self.cm)+'cm'+'\n')
            of.write('#MJD0 time freq flux flux_err toa_err/sqrt(intTime) filename\n')
            for dat in self.sortedData:
                for entry in dat:
                    of.write(str(entry)+' ')
                of.write('\n')
        print 'write time series file '+self.timeSeriesFileName


class structure_function:
   # timeSeriesFile
   # header
   # mjd
   # flux
    def __init__(self,timeSeriesFile):
        self.mjd=[]
        self.flux=[]
        self.fluxErr=[]
        self.timeSeriesFile=timeSeriesFile
        self.valid=True
        self.sf=[]
        with open(timeSeriesFile,'r') as f:
            firstLine=True
            for line in f:
                if firstLine:
                    self.header=line
                    firstLine=False
                else:
                    if line[0]!='#':
                        self.mjd.append(float(line.split()[0]))
                        self.flux.append(float(line.split()[3]))
                        self.fluxErr.append(float(line.split()[4]))

    def create_equal_weight_sf(self,outputDir,info,binNum):
        # Note that this method INDEPENDENTLY creates a sf file. No need to call create_sf
        #create equal weight sf. Each bin has almost the same number of points
        # binNum is the number of bins
        numDat=len(self.flux)
        if numDat==0:
            self.valid=False
            return
        numPairs=numDat*(numDat-1)/2
        unbinSf=np.zeros((3,numPairs))#(lag,sf Value,sfErr)
        i=0
        k=0
        outFileName=outputDir+'/'+info+'.sf'
        while i<numDat:
            j=i+1
            while j<numDat:
                unbinSf[0][k]=self.mjd[j]-self.mjd[i]
                fluxI=self.flux[i]
                fluxJ=self.flux[j]
                errI=self.fluxErr[i]
                errJ=self.fluxErr[j]
                # sf Value for the point
                unbinSf[1][k]=(fluxI-fluxJ)**2-(errI**2+errJ**2)
                # sf uncertainty (in stadard deviation)
                unbinSf[2][k]=math.sqrt(
                        (
                            2*(errI**4+errJ**4)+
                            4*(errI**2+errJ**2)*((fluxI-fluxJ)**2)+
                                4*(errI**2)*(errJ**2)
                                )
                        )
                #print unbinSf[2][k]
                        
                k+=1
                j+=1
            i+=1
        assert k==numPairs,'Check expression for numPairs,the calc isnt right'
        #NEED TO SORT unbinsf BEFORE MOVING ON
        sortIndex=np.argsort(unbinSf[0])

        if binNum<numPairs and binNum!=0:
            self.binned=True
        else:
            self.binned=False
            binNum=numPairs#so that things are unbinned

        binSize=numPairs/binNum
        self.sf=np.zeros((5,binNum))#[average lag,average sf value, sf uncertainty, SD of points in bin, points within bin]
        #cumsum
        j=0
        self.binArr=[]
        for i in sortIndex:
            binIndex=j/binSize
            if binIndex>=binNum:
                binIndex=binNum-1#the last bin will have at least as many data points as all other
            self.sf[0][binIndex]+=unbinSf[0][i]
            self.sf[1][binIndex]+=unbinSf[1][i]
            self.sf[2][binIndex]+=unbinSf[2][i]**2
            self.sf[3][binIndex]+=unbinSf[1][i]**2
            self.sf[4][binIndex]+=1
            j+=1
        """    if binIndex==50:
                self.binArr.append(unbinSf[1][i])"""
        #average lag and sf value
        #print [len(self.binArr),self.sf[0][4]/self.sf[4][4],np.mean(self.binArr),sum(self.binArr)/self.sf[4][4]]
        self.sf[0]=np.divide(self.sf[0],self.sf[4])
        self.sf[1]=np.divide(self.sf[1],self.sf[4])
        # See You et al 2007. I derived and confirmed. One should not use an rms.
        self.sf[2]=np.divide(np.sqrt(self.sf[2]),self.sf[4])
        #print self.sf[1][2]
        # Finally, estimate the standard deviation within each bin
        self.sf[3]=np.sqrt(np.divide(self.sf[3],self.sf[4])-np.square(self.sf[1]))
        #meanFluxSquare=np.mean(self.flux)**2
        #self.sf[1]=self.sf[1]/meanFluxSquare
        if self.binned:
            outFileName=outFileName+'.binned'
        #I really don't like this....
        np.savetxt(outFileName,np.transpose(self.sf))

class acf_fft:
    #should use interpolate.interp1d
    #np.correlate
    def __init__(self,timeSeriesFile):
        self.timeSeriesFile=timeSeriesFile
        mjd=[]
        flux=[]
        with open(timeSeriesFile,'r') as f:
            for line in f:
                if line[0]!='#':
                    l=line.split()
                    mjd.append(float(l[0]))
                    flux.append(float(l[3]))
        self.data=np.array([mjd,flux])

    def create_acf(self,interval=33):
        # first interpolate the time series
        self.interval=interval
        interpTimeSeries=interpolate.interp1d(self.data[0],self.data[1])
        x=np.arange(self.data[0].min(),self.data[0].max(),self.interval)
        self.equalSampleTimeSeries=interpTimeSeries(x)
        m=np.mean(self.equalSampleTimeSeries)
        acf=np.correlate(self.equalSampleTimeSeries-m,self.equalSampleTimeSeries-m,"full")
        self.untruncAcf=np.divide(acf,acf.max())
        #Truncate the acf
        truncAcf=self.untruncAcf[(self.untruncAcf.shape[0]-1)/2:]
        self.truncAcf=np.array([np.arange(0,truncAcf.size*self.interval,self.interval),truncAcf])

    def get_tau(self,threshold=0.2):
        for i in range(0,self.truncAcf[1].size-1,1):
            if self.truncAcf[1][i]>=threshold and self.truncAcf[1][i+1]<=threshold:
                # linear interpolate between the two points and solve the quation
                y1=self.truncAcf[1][i]
                y2=self.truncAcf[1][i+1]
                x1=self.truncAcf[0][i]
                x2=self.truncAcf[0][i+1]
                print [x1,x2,y1,y2]
                return (threshold*x1-threshold*x2+x2*y1-x1*y2)/(y1-y2) 
        # Can't find a solution
        return 0.

    def time_series_power_spectrum(self):
        # FFT the interpolated time series to get the power spectrum of time series
        timeSeriesPowerSpec=np.absolute(fftpack.fft(self.equalSampleTimeSeries-np.mean(self.equalSampleTimeSeries)))**2
        timeSeriesFreqs=fftpack.fftfreq(self.equalSampleTimeSeries.size,float(self.interval))
        if timeSeriesPowerSpec.size % 2==0:
            self.truncPowerSpec=timeSeriesPowerSpec[:timeSeriesPowerSpec.size/2]
            self.truncFreqs=timeSeriesFreqs[:timeSeriesPowerSpec.size/2]
        else:
            self.truncPowerSpec=timeSeriesPowerSpec[:(timeSeriesPowerSpec.size-1)/2:]
            self.truncFreqs=timeSeriesFreqs[:(timeSeriesPowerSpec.size-1)/2]
        self.maxPeakFreq=self.truncFreqs[np.argmax(self.truncPowerSpec)]

    def acf_linear_spectrum(self):
        # FFT the untruncated ACF to get the linear spectrum of ACF
        self.fft=fftpack.fft(self.untruncAcf)
        self.powerSpec=np.absolute(self.fft)
        fftFreqs=fftpack.fftfreq(self.fft.size,float(self.interval))
        print fftFreqs[np.argmax(self.powerSpec)]
        plt.figure()
        plt.plot(fftFreqs,self.powerSpec,marker='o')
        plt.show()

def main(allDirs):
    firstDir='/DATA/RUBICON_2/joh414/yuping/J1048-5832/'
    tempdir='/DATA/NEWTON_1/pulsar/ker14a/P574/test/old_p574_pipeline_templates'
    plotdir='/DATA/RUBICON_2/joh414/yuping/time_series_plots'
    dmList={}
    plt.ioff()
    with open('DMList.txt','r') as f:
        for line in f:
            l=line.split()
            dmList.update({l[0]:l[1]})

    if allDirs=='':
        allDirs=glob.glob('/DATA/RUBICON_2/joh414/yuping/test/J*')
    for testdir in allDirs:
        t=time_series(testdir,20,1024,tempdir,overwrite=False)
        s=structure_function(t.directory+'/'+t.timeSeriesFileName)
        s.create_equal_weight_sf(t.directory,t.name+'.'+str(t.cm)+'cm',100)
        if s.valid:
            fig,ax=plt.subplots(2,3)
            fig.set_size_inches(16,11)
            fig.suptitle(t.name+' DM= '+dmList[t.name])
            # Structure Function Plot
            [ax1,ax2,ax5],[ax3,ax4,ax6]=ax
            ax1.set_title('Binned Structure Function')
            ax1.set_yscale('log')
            ax1.set_xscale('log')
            ax1.set_xlabel('lag(Days)')
            ax1.set_ylabel('Structure Function D')
            ax1.errorbar(s.sf[0],s.sf[1],yerr=s.sf[2],marker='o',ms=4)

            # Time Series
            ax2.set_title('Flux Time Series')
            ax2.set_xlabel('MJD-54000')
            ax2.set_ylabel('Flux(mJy)')
            ax2.errorbar([float(mjd)-54000 for [mjd,time,freq,flux,flux_err,toa_err,fn,gof,scaledNoise] in t.sortedData],
                    [float(flux) for [mjd,time,freq,flux,flux_err,toa_err,fn,gof,scaledNoise] in t.sortedData],
                    yerr=[float(flux_err) for [mjd,time,freq,flux,flux_err,toa_err,fn,gof,scaledNoise] in t.sortedData],
                    ecolor='r',marker='o',ms=4)

            # Autocorrelation Function
            acfRef=acf_fft(t.directory+'/'+t.timeSeriesFileName)
            acfRef.create_acf()
            ax3.set_title('acf'+' decor time= '+str(acfRef.get_tau()))
            ax3.set_xlabel('Lag/(Days)')
            ax3.set_ylabel('Autocorrelation')
            ax3.plot(acfRef.truncAcf[0],acfRef.truncAcf[1],marker='o',ms=5)

            # Goodness of fit
            ax4.set_title('Goodness of Fit\n(0 means no corresponding TOA)')
            ax4.set_xlabel('MJD-54000')
            ax4.plot([float(mjd)-54000 for [mjd,time,freq,flux,flux_err,toa_err,fn,gof,scaledNoise] in t.sortedData],
                    [float(gof) for [mjd,time,freq,flux,flux_err,toa_err,fn,gof,scaledNoise] in t.sortedData],marker='o',ms=4)

            #Power spectrum of the time series
            acfRef.time_series_power_spectrum()
            ax5.set_title('Power Spectrum of the Time Series\n'+'peak:  '+str(1./acfRef.maxPeakFreq)+' days')
            ax5.set_xlabel('1/(1000 Days)')
            ax5.plot(acfRef.truncFreqs*1000,acfRef.truncPowerSpec,marker='o')

            # Noise baseline of each scan multiplied by srt of integration time. Should be only dependent on the setup
            ax6.set_title('off:rms*sqrt(integration time)')
            ax6.set_xlabel('MJD-54000')
            ax6.set_ylabel('mJy*sqrt(t)')
            ax6.plot([float(mjd)-54000 for [mjd,time,freq,flux,flux_err,toa_err,fn,gof,scaledNoise] in t.sortedData],
                    [float(scaledNoise) for [mjd,time,freq,flux,flux_err,toa_err,fn,gof,scaledNoise] in t.sortedData],marker='o',ms=4)

            fig.savefig(plotdir+'/'+t.name+str(t.cm)+'.cm.time_series.png')
            fig.clf()

def time_series_all():
    allDirs=glob.glob('/DATA/RUBICON_2/joh414/yuping/test/J*')
    tempdir='/DATA/NEWTON_1/pulsar/ker14a/P574/test/old_p574_pipeline_templates'
    #allDirs=[firstDir]
    for testdir in allDirs:
        try:
            t=time_series(testdir,20,1024,tempdir,overwrite=False)
        except subprocess.CalledProcessError,e:
            with open('timeSeriesErrorLog','a') as errorLog:
                z=e
                errorLog.write(str(z)+'\n')
