# summer_15
# Below are the useful interfaces of various class
# Yuping Huang huangy@carleton.edu Summer 2015

time_series.py
"Main program that deal with generating time series, structure function, acf and power spec. Also does plotting.

--------------------------
    main(arg):
        " Do the batch processing for all data. Generate new/load existing time series files, create structure functions, acf, power spec, and plot and output. Good place to see the example usage of various classes in this file.
        when arg=='' It will generate the plots for all files.
        when arg is an array of full paths of data directories, then the processing will only be done to the directories mentioned.
---------------------------
    class time_series
        self.sortedData: a splitted array of strings of the entries in the time series. Available after initializing.
        self.name: name of the pulsar. Available after initializing
        self.directory: PATH of the time series file without file name
        self.timeSeriesFileName: NAME of the time series file

        __init__(self,directory,cm,nbin,stdProfDir,overwrite=False):
            "Generate the time series file if it does not exist. If it exist, load the file into self.sortedData.
            @str directory: the pulsar scans directory. The output time series file will be in this directory as well
            @int cm: observing wavelength in cm. Can choose 10,20,50
            @int nbin: number of bins, should match the nbin of template. Usually 1024
            @str stdProfDir:Directory of the std profiles. Concatenated to get the template file name in method group_psrflux
            @boolean overWrite=False: overwrite existing time_series file. Default to False
-----------------------------
    class structure_function
        self.sf: 4*binNum np array. The four dimensions are respectively mean lag of the bin, mean structure function value, structure function value std, and number of points in the bin. Available after running 

        __init__(self,timeSeriesFile):
            "Initialize and load the time series file for processing
            @str timeSeriesFile: full path of the time series file.
            
        create_equal_weight_sf(self,outputDir,info,binNum):
            " Compute the binned structure function and write to file.
            @str outputDir: the output folder of the structure function file
            @str info: extensions to the output filename. Output file will be in the format $info.sf.binned or $info.sf.unbinned (you may just want to parse in the pulsar name here)
            @int binNum: number of bins wanted. If set to 0 or larger than the number of sf values, then the structure function will be an unbinned one and the file name will be $info.sf.unbinned.

-----------------------------
    class acf_fft
        self.truncAcf: 2D np array. [lag,acfValue] The acf with not negative lag. Available after create_acf()
        self.truncPowerSpec: 1D np array. The power spectrum of the time series with non-negative frequency.
        self.truncFreqs: 1D np array. The frequencies corresponding to the power in self.truncPowerSpec. In unit 1/Days. Hence you want to plot truncFreqs on x and truncPowerSpec on y. Both available after time_series_power_spectrum.
        self.maxPeakFreq: The frequency corresponding to the max peak in the time series power spectrum. In unit 1/Days.

        __init__(self,timeSeriesFile):
            @str timeSeriesFile the time series file full path

        create_acf(self[,interval=33]):
            "create the acf and write the truncated version to self.truncAcf
            @int interval: resampling interval. Default to 33.

        float get_tau(self,[threshold=0.2]):
            "Return the timescale when the acf drops below the threshold. Return 0. if not found.
            @float threshold: threshold value for the time scale, default to 0.2.

        time_series_power_spectrum(self):
            " Create the power spectrum for the time series and store the positive frequency part in [self.truncFreqs,self.truncPowerSpec].

=================================================================
edit_time_series.py
"very useful interface to visually examine data points and clean up.
"Also time_series_file is my interface for editting time series files.

--------------------
    class add_par
        "very incomplete and ugly procedures to add more column entries to the time series files. However it is a good example of using the time_series_file interface.

--------------------
    class time_series_file
        " Main interface if you want to edit a existing time series file, including adding entries and especially cleaning up data.
        self.fileBuffer: the file buffer. An array of all lines in the time series file (including the headers and the \n character at the end of each line) IN STRINGS. Any changes to the time series file should first be applied to this buffer. Then use unload() to write it to file.

        __init__(self,filename):
            " load the time series file in to self.fileBuffer
            @str filename: full path of the time series file.

        view_scan(self,mjd):
            " run pav to view the scan corresponding to the MJD time
            @float mjd: mjd of the data point to view. Should be exactly the same as the mjd in the time series file. It is easy to get the MJD value from the fileBuffer

        comment_out_mjd(self,mjd)
        @float mjd
        " comment out the data entry in the time series file corresponding to the input mjd. Note that mjd should be exactly the same as the mjd in the file. Do nothing if the mjd is not found. The reason for using mjd instead of index is that mjd is much less confusing especially when one has another data array with floats for some other purpose. The indices may be different from that array and the file buffer because the file buffer contains commented out lines and the headers. THIS WILL NOT WRITE TO THE FILE. IT ONLY CHANGE THE BUFFER

        unload()
        "write the buffer to the file.

-------------------
    class plot_interface
        "COOL STUFF. RECOMMENDED INTERFACE TO MANUALLY CLEAN UP DATA. plot a time series file and allow interactive viewing and commenting out of data points.
        __init__(self):
            initializer that require no arguments.

        add_file_list(self,fileList):
            "add the LIST of files to edit.
            @fileList: a list of full paths of the files to edit. If editing a single file. Just all the square brackets around it.

       create_plot(self):
           create a time series plot using matplotlib on which you can do the following:
            1. click on a point: view the scan with pav (after viewing you want to reactivate the pyplot window by clicking on it to issue other commands).
            2. press key 'n' or 'p': go to next or previous file. Python will complain when at the end of the list.
            3. press key 'd':delete the point just viewed. WILL NEED TO PRESS 'w' TO SAVE THE FILE.
            4. press key 'w': save the changes.


================================
clean_data.py
" identifying outliers in the dataset according to set metrics. The time_series_file interface is used. See the main() method for the process. Writes to the log file outliersCleaning.log.

class outliers_id: Do the job. initializer takes in the full path of a time series file.
