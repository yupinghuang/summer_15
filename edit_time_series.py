from matplotlib import pyplot as plt
import math
import glob
import subprocess
import psrchive

class time_series_file:
    def __init__(self,filename):
        self.filename=filename
        # fileBuffer is the place to store the data from the file. All entries are in string
        # eventually fileBuffer is written to the output file
        self.fileBuffer=[]
        self.currentInd=-1#index pointing to the current line being viewed&to be deleted
        self.mjdList=[]
        with open(self.filename,'r') as f:
            for line in f:
                    self.fileBuffer.append(line)
                    if line[0]=='#':
                        # mjdList is for lookup purpose
                        self.mjdList.append(-99)
                    else:
                        self.mjdList.append(float(line.split()[0]))

    def get_index_from_mjd(self,mjd):
        # float:MJD
        # returns the list index of the record with a given mjd in self.fileBuffer
        # if can't find, return -1
        index=-1
        for i in range(0,len(self.mjdList)):
            if mjd==self.mjdList[i]:
                index=i
                self.currentInd=i
                break
        return index

    def view_scan(self,mjd):
        # float:MJD
        # given a MJD in float, view the scan with pav -D
        index=self.get_index_from_mjd(mjd)
        if index==-1:
            print "Can't find corresponding entry for MJD "+str(mjd)
            return
        scanFileName=self.fileBuffer[index].rstrip('\n').split()[6]
        subprocess.check_output('pav -CD '+scanFileName,shell=True)

    def comment_out(self):
        print 'deleting data point: '+self.fileBuffer[self.currentInd]
        if self.currentInd>-1:
            self.fileBuffer[self.currentInd]='#'+self.fileBuffer[self.currentInd]
            self.currentInd=-1

    def comment_out_mjd(self,mjd):
        self.get_index_from_mjd(mjd)
        self.comment_out()

    def unload(self):
        #ehhhhhhhhhhhhh 
        with open(self.filename,'w') as of:
            for line in self.fileBuffer:
                of.write(line)
        print 'file saved'


class plot_interface:
    def __init__(self):
        self.mjd=[]
        self.flux=[]
        self.flux_err=[]
        self.refresh=False
        self.filename=''

    def add_file_list(self,fileList):
        self.fileList=fileList
        self.currentFileIndex=0
        self.load_from_file(self.fileList[self.currentFileIndex])

    def load_from_file(self,filename):
        print 'loading '+filename
        self.mjd=[]
        self.flux=[]
        self.flux_err=[]
        self.filename=filename
        # read the columns
        # pass the file to time_series_file for manipulation
        self.timeSeries=time_series_file(filename)
        for entry in self.timeSeries.fileBuffer:
            #READ FROM DATA FILE 
            if entry[0]!='#':
                splittedLine=entry.rstrip('\n').split()
                self.mjd.append(float(splittedLine[0]))
                self.flux.append(float(splittedLine[3]))
                self.flux_err.append(float(splittedLine[4]))

    def reload_data(self):
        #reload after changes made to the file buffer
        self.mjd=[]
        self.flux=[]
        self.flux_err=[]
        #reload_data from the existing buffer
        for entry in self.timeSeries.fileBuffer:
            #READ FROM DATA FILE 
            if entry[0]!='#':
                splittedLine=entry.rstrip('\n').split()
                self.mjd.append(float(splittedLine[0]))
                self.flux.append(float(splittedLine[3]))
                self.flux_err.append(float(splittedLine[4]))
        self.refresh=True

    def on_pick(self,event):
        self.timeSeries.view_scan(event.artist.name)

    def next_file(self):
        self.currentFileIndex+=1
        fn=self.fileList[self.currentFileIndex]
        self.load_from_file(fn)

    def prev_file(self):
        self.currentFileIndex-=1
        fn=self.fileList[self.currentFileIndex]
        self.load_from_file(fn)

    def on_press(self,event):
        if event.key=='n':
            self.next_file()
            self.create_plot()
        if event.key=='d':
            self.timeSeries.comment_out()
            self.reload_data()
            self.create_plot()
        if event.key=='w':
            self.timeSeries.unload()
        if event.key=='q':
            plt.close()
            self.refresh=False
        if event.key=='p':
            self.prev_file()
            self.create_plot()

    def create_plot(self):
        plt.ion()
        if self.refresh:
            self.fig.clear()
            plt.figure(1)
        else:
            self.fig=plt.figure(1)
            self.fig.canvas.callbacks.connect('pick_event',self.on_pick)
            self.fig.canvas.mpl_connect('key_press_event',self.on_press)
        for x,y in zip(self.mjd,self.flux):
            pt,=plt.plot(x,y,'.b',picker=3)
            pt.name=x
        plt.errorbar(self.mjd,self.flux,yerr=self.flux_err)
        if not self.refresh:
            plt.show()
        self.refresh=True


#################################################
class add_par:
    # Note that all methods in this class only add parameters to the LAST column of the  time series file.
    # Hence add_gof should be done before add_scaledOffRms for a file without any of these additional columns.
    # `
    def __init__(self,timeSeriesFile,toaFile):
        self.t=time_series_file(timeSeriesFile)
        self.toaFile=toaFile

    def add_gof(self):
        gofList={}
        with open(self.toaFile,'r') as f:
            for line in f:
                l=line.rstrip('\n').split()
                gofList.update({l[0]:l[-1]})

        if self.t.fileBuffer[1].rstrip('\n').split()[-1]!='gof':
            self.t.fileBuffer[1]=self.t.fileBuffer[1].rstrip('\n')+' gof\n'
        for i in range(2,len(self.t.fileBuffer)):
            scanName=self.t.fileBuffer[i].rstrip('\n').split()[-1].split('/')[-1]
            try:
                gof=gofList[scanName]
                self.t.fileBuffer[i]=self.t.fileBuffer[i].rstrip('\n')+' '+gof+'\n'
            except KeyError:
                print "Can't find TOA for "+self.t.fileBuffer[i]
                self.t.fileBuffer[i]=self.t.fileBuffer[i].rstrip('\n')+' '+'0.0'+'\n'
        self.t.unload()

    def add_scaledOffRms(self):
        if self.t.fileBuffer[1].rstrip('\n').split()[-1]!='off:rms*sqrt(tInt)':
            self.t.fileBuffer[1]=self.t.fileBuffer[1].rstrip('\n')+' off:rms*sqrt(tInt)\n'
        for i in range(2,len(self.t.fileBuffer)):
            scanFile=self.t.fileBuffer[i].rstrip('\n').split()[-2]
            [fn,tInt,offRms]=subprocess.check_output('psrstat -c "length","off:rms" '+scanFile+' -Q',shell=True).split()
            self.t.fileBuffer[i]=self.t.fileBuffer[i].rstrip('\n')+' '+str( float(offRms)*math.sqrt(float(tInt)))+'\n'
        self.t.unload()

def add_gof_to_all():
    psrList=[]
    with open('DMList.txt','r') as f:
        for line in f:
            psrList.append(line.split()[0])
    for psr in psrList:
        g=add_par('/DATA/RUBICON_2/joh414/yuping/test/'+psr+'/'+psr+'.20.time_series','/DATA/NEWTON_1/pulsar/ker14a/P574/test/'+psr+'/'+psr+'.toa')
        g.add_gof()
        g=add_par('/DATA/RUBICON_2/joh414/yuping/test/'+psr+'/'+psr+'.20.time_series','/DATA/NEWTON_1/pulsar/ker14a/P574/test/'+psr+'/'+psr+'.toa')
        g.add_scaledOffRms()
