#plotData.py

""" This program will read in a data file from Gnu Radio, apply 
    temperature corrections, and plot it.
    Ellie White, Summer 2019
"""

import pylab as plt
import os
import numpy as np
import numpy.polynomial.polynomial as poly
from tempCorrections import tempCorrections
#from assessCorrections import assessCorrs
from astropy.time import Time
import time
import datetime

def main():

    # Filenames
    infile_a = 'DIAG_DATA_28AUG19_1A.obs' 
    infile_b = 'DIAG_DATA_28AUG19_1B.obs'
    #in_tfilename = 'TIMEFILE_28AUG19_1.txt'
    out_tfilename = 'timestamps.txt'
    rpi_filename = 'RPI_DATA_08_28_2019_20-41-35.txt'

    #ofname_a = 'outfile_6a.txt'
    #ofname_b = 'outfile_6b.txt'

    # Define some constants
    fftsize = 8192
    tot_time = 86400 #observation duration, in seconds
    srate = 4000000.0 
    decimation = 3276.0
    cfreq = 169010000

    rpi_samptime = 300 # number of seconds between sensor readings

    K_CONST = 1.3806*pow(10,-23) #Boltzmann constant

    cable_gain = -7.5
    filter_gain = -0.5
    ettus_gain = 78.2824
    balun_gain = 30
    rx_gain = 57.5



    """ Read in observation data """

    size_a = os.path.getsize(infile_a) / 4
    length_a = size_a/fftsize
    shape_a = (int(length_a), fftsize)
    #print(str(length_a))
    xa = np.memmap(infile_a, dtype='float32', mode = 'r', shape=shape_a) 

    size_b = os.path.getsize(infile_b) / 4
    length_b = size_b/fftsize
    shape_b = (int(length_b), fftsize)
    xb = np.memmap(infile_b, dtype='float32', mode = 'r', shape=shape_b)



    """ Create bandpass plot """

    freqPlotA = np.mean(xa, axis=0)
    freqPlotB = np.mean(xb, axis=0)
    fmin = (169010000-(srate/2))/1000000
    fmax = (169010000+(srate/2))/1000000
    fidx = np.linspace(fmin, fmax, freqPlotA.size)

    for f in range(freqPlotA.size):
        freqPlotA[f] = (freqPlotA[f]/decimation)#pow(freqPlot[f]*-17.8398111277, 2)/50.0
        freqPlotB[f] = (freqPlotB[f]/decimation)



    """ Read in data from RPI file as array
        Balun temps and Rx temps:                         """

    tc_b_chanA = tempCorrections('balun', 1, rpi_filename, 0)
    tc_rx_chanA = tempCorrections('rx', 1, rpi_filename, 1)

    tc_b_chanB = tempCorrections('balun', 2, rpi_filename, 0)
    tc_rx_chanB = tempCorrections('rx', 2, rpi_filename, 1)
    
    #convert temperature array to temperature list
    #balun temps:
    balun_temps_A = []
    btemps_A = tc_b_chanA.getTemps()

    balun_temps_B = []
    btemps_B = tc_b_chanB.getTemps()

    for ba in btemps_A:
        balun_temps_A.append(ba)

    for bb in btemps_B:
        balun_temps_B.append(bb)

    #Rx temps:
    rx_temps_A = []
    rxtemps_A = tc_rx_chanA.getTemps()

    rx_temps_B = []
    rxtemps_B = tc_rx_chanB.getTemps()

    for ra in rxtemps_A:
        rx_temps_A.append(ra)

    for rb in rxtemps_B:
        rx_temps_B.append(rb)



    """ Create time arrays """
    
    #amount of time between data bursts...
    delta_t = tot_time / length_a

    #read times from RPI file...
    times = tc_b_chanA.getTimes() #times = RPI times
    start_time = times[0] - rpi_samptime
 
    #create array of times in usual Unix format from RPI time
    t = 0
    unix_times = []

    for t in times:
        time = t + 654690.1978 #this offset is an arbitrary value inserted to correct for RPI time
        unix_times.append(time)

    unix_times = np.array(unix_times)

    #Create array of LST timestamps from RPI times 
    t = Time(unix_times, format='unix', scale='utc')
    lst = t.sidereal_time('apparent', '-79.8d')


    """ Assign temperature reading to each spectra """

    data_btemps_a = []
    data_rxtemps_a = []
    ant_power_a = []

    data_btemps_b = []
    data_rxtemps_b = []
    ant_power_b = []

    for i in range(int(size_a/fftsize)):

        #assign sensor row to antenna data row
        row_time = start_time + (i+1)*delta_t #antenna data time

        #channel A
        ta_val = min(times, key=lambda xa:abs(xa-row_time))
        ta_index = times.index(ta_val)

        #get temperatures for balun and rx from sensor row
        btemp_ia = balun_temps_A[ta_index]
        data_btemps_a.append(btemp_ia)

        rxtemp_ia = rx_temps_A[ta_index] 
        data_rxtemps_a.append(rxtemp_ia)

        #channel B
        tb_val = min(times, key=lambda xb:abs(xb-row_time))
        tb_index = times.index(tb_val)

        #get temperatures for balun and rx from sensor row
        btemp_ib = balun_temps_B[tb_index]
        data_btemps_b.append(btemp_ib)

        rxtemp_ib = rx_temps_B[tb_index] 
        data_rxtemps_b.append(rxtemp_ib)


    ''' Create temperature model 2-D array '''

    s12mod_mag_ba, s12mod_ph_ba = tc_b_chanA.getS12MagPhase() 
    s12mod_mag_bb, s12mod_ph_bb = tc_b_chanB.getS12MagPhase()

    s12mod_mag_ra, s12mod_ph_ra = tc_rx_chanA.getS12MagPhase() 
    s12mod_mag_rb, s12mod_ph_rb = tc_rx_chanB.getS12MagPhase()


    """ Take a 4 MHz slice from the temp model """

    start_freq = cfreq - (srate/2)
    freqs = []

    for q in range(fftsize):
        freq = start_freq + ((srate/fftsize)*q)
        freqs.append(freq)
    
    freq_array = np.array(freqs)

    fvals = tc_b_chanA.gc.returnFreqs()
    fvals_slice = fvals[47:53]*1000000

    s12mag_ba_slice = s12mod_mag_ba[:, 47:53]
    s12mag_bb_slice = s12mod_mag_bb[:, 47:53]
    s12mag_ra_slice = s12mod_mag_ra[:, 47:53]
    s12mag_rb_slice = s12mod_mag_rb[:, 47:53]

    mod_curves_ba = np.empty([int(size_a/fftsize), freqPlotA.size])
    mod_curves_bb = np.empty([int(size_b/fftsize), freqPlotA.size])
    mod_curves_ra = np.empty([int(size_a/fftsize), freqPlotA.size])
    mod_curves_rb = np.empty([int(size_b/fftsize), freqPlotA.size])

    bgain_vals_a = []
    bgain_vals_b = []
    rxgain_vals_a = []
    rxgain_vals_b = []

    #call tempCorrections with given temperature to get model curve
    for k in range(len(data_btemps_a)):
        btemp_a = data_btemps_a[k]
        index_ba = balun_temps_A.index(btemp_a)
        mod_curves_ba[k] = np.interp(freq_array, fvals_slice, s12mag_ba_slice[index_ba]) 

        btemp_b = data_btemps_b[k]
        index_bb = balun_temps_B.index(btemp_b)
        mod_curves_bb[k] = np.interp(freq_array, fvals_slice, s12mag_bb_slice[index_bb]) 

        rxtemp_a = data_rxtemps_a[k]
        index_ra = rx_temps_A.index(rxtemp_a)
        mod_curves_ra[k] = np.interp(freq_array, fvals_slice, s12mag_ra_slice[index_ra]) 

        rxtemp_b = data_rxtemps_b[k]
        index_rb = rx_temps_B.index(rxtemp_b)
        mod_curves_rb[k] = np.interp(freq_array, fvals_slice, s12mag_rb_slice[index_rb]) 

        #convert model curve to power
        #apply any offsets to slice of row i data, same slice as used for timeseries plot
        # i.e. ...

        corrlist_ba = []
        corrlist_bb = []
        corrlist_ra = []
        corrlist_rb = []

        for c in range(100):
          corrlist_ba.append(mod_curves_ba[k, (4400+c)])
          corrlist_bb.append(mod_curves_bb[k, (4400+c)])
          corrlist_ra.append(mod_curves_ra[k, (4400+c)])
          corrlist_rb.append(mod_curves_rb[k, (4400+c)])

        carray_ba = np.array(corrlist_ba)
        cval_ba = np.mean(carray_ba)
        bgain_vals_a.append(cval_ba) 

        carray_bb = np.array(corrlist_bb)
        cval_bb = np.mean(carray_bb) 
        bgain_vals_b.append(cval_bb) 

        carray_ra = np.array(corrlist_ra)
        cval_ra = np.mean(carray_ra) 
        rxgain_vals_a.append(cval_ra)  

        carray_rb = np.array(corrlist_rb)
        cval_rb = np.mean(carray_rb)
        rxgain_vals_b.append(cval_rb) 


    ''' Apply corrections and create timeseries plot '''

    ant_temps_a = [] 
    ant_temps_b = []

    final_vals_a = []
    final_vals_b = []

    nominal_gain = balun_gain + cable_gain + rx_gain + filter_gain + ettus_gain


    for j in range(int(size_a/fftsize)):

        alist_a = []
        alist_b = []
        flist_a = []
        flist_b = []

        for w in range(100):
            atemp_a = (xa[j, (2400+w)])/(K_CONST*(srate/float(fftsize))*100*pow(10, nominal_gain/10.0))
            atemp_b = (xb[j, (2400+w)])/(K_CONST*(srate/float(fftsize))*100*pow(10, nominal_gain/10.0))
            alist_a.append(atemp_a)
            alist_b.append(atemp_b)

            dB_A = ((bgain_vals_a[j]+ 3) + cable_gain + rxgain_vals_a[j] + filter_gain + ettus_gain)
            ftemp_a = (xa[j, (2400+w)])/(K_CONST*(srate/float(fftsize))*100*pow(10, dB_A/10.0))
            flist_a.append(ftemp_a)

            dB_B = ((bgain_vals_b[j] + 3) + cable_gain + rxgain_vals_b[j] + filter_gain + ettus_gain)
            ftemp_b = (xb[j, (2400+w)])/(K_CONST*(srate/float(fftsize))*100*pow(10, dB_B/10.0))
            flist_b.append(ftemp_b)

        a_array_a = np.array(alist_a)
        a_array_b = np.array(alist_b)

        f_array_a = np.array(flist_a)
        f_array_b = np.array(flist_b)
 
        aval_a = np.sum(a_array_a)
        ant_temps_a.append((aval_a/decimation) - 488.0)

        aval_b = np.sum(a_array_b)
        ant_temps_b.append((aval_b/decimation) - 488.0)

        fval_a = np.sum(f_array_a)
        final_vals_a.append((fval_a/decimation) - 488.0)

        fval_b = np.sum(f_array_b)
        final_vals_b.append((fval_b/decimation) - 488.0)

    tempslength_A = len(balun_temps_A)
    print("Length: {0}".format(tempslength_A))

    
    ''' Plot timeseries graphs '''

    #print(final_vals_a)

    at_len = len(ant_temps_a)
    tlinespace = np.linspace(0, (tempslength_A - 1), tempslength_A)
    alinespace = np.linspace(0, (tempslength_A - 1), at_len)

    # create Local Sidereal Time labels for the x-axis
    linespace_locs = []
    lst_labels = []
    num_of_labels = 15
    

    for l in range(num_of_labels):
        lindex = int((tempslength_A/num_of_labels)*l)
        lst_labels.append(lst[lindex])
        linespace_locs.append(tlinespace[lindex])

    lst_labels.append(lst[-1])
    linespace_locs.append(tlinespace[-1])

    '''#Radiometer analysis
    radiometer_data = np.array(ant_temps_a[110:121])
    rad_coeffs = poly.polyfit(alinespace[110:121], radiometer_data, 1)
    rad_fit = alinespace[110:121]*rad_coeffs[1] + rad_coeffs[0]

    rad_diff = radiometer_data - rad_fit

    stderr = np.std(rad_diff)
    print(stderr)

    plt.plot(rad_fit)
    plt.plot(radiometer_data)
    #plt.plot(rad_diff)
    plt.show()'''


    ''' Plot resulting data graphs '''

    # Plot the timeseries graphs ant_temps_a
    #plt.plot(alinespace, ant_temps_a, label="Antenna Temp A", color='b')
    plt.plot(alinespace, final_vals_a, label="Corrected Ant. Temps A", color='g')
    #plt.plot(tlinespace, balun_temps_A, label = "Physical Temp", color='r')
    #plt.plot(alinespace, final_vals, label = "Temp. Corrections", color='r')
    plt.xticks(linespace_locs, lst_labels, rotation=90)
    plt.legend(loc='upper right')
    plt.title(infile_a[:-4])
    plt.subplots_adjust(bottom=0.28)
    plt.ylabel("Temperature [K]")
    plt.xlabel("Local Sidereal Time")
    plt.show()

    plt.plot(alinespace, final_vals_b, label="Corrected Ant. Temps B", color='b')
    #plt.plot(alinespace, final_vals_b, label="Physical Temp", color='g')
    #plt.plot(tlinespace, balun_temps_A, label = "Physical Temp", color='g')
    #plt.plot(alinespace, final_vals, label = "Temp. Corrections", color='r')
    plt.xticks(linespace_locs, lst_labels, rotation=90)
    plt.legend(loc='upper right')
    plt.title(infile_b[:-4])
    plt.subplots_adjust(bottom=0.28)
    plt.ylabel("Temperature [K]")
    plt.xlabel("Local Sidereal Time")
    plt.show()

    '''humidities = tc_b_chanA.getHumidity()
    plt.plot(humidities)
    plt.show()

    # Plot the timeseries graphs, with temp. corrections applied
    #plt.plot(alinespace, ant_temps_a, label='Channel A', color='b')
    #plt.plot(alinespace, ant_temps_b, label='Channel B', color='g')
    plt.plot(tlinespace, balun_temps_A, label = "Physical Temp", color='r')
    plt.xticks(linespace_locs, lst_labels, rotation=90)
    plt.legend(loc='upper right')
    plt.title(infile_a[:-7])
    plt.subplots_adjust(bottom=0.28)
    plt.ylabel("Temperature [K]")
    plt.xlabel("Local Sidereal Time")
    plt.show()'''

    # Show the averaged bandpass plot
    plt.plot(freqPlotA)
    #plt.plot(freqPlotB)
    plt.show()

    '''outfile_a = open(ofname_a, 'w')
    outfile_b = open(ofname_b, 'w') 


    for ind in range(len(ant_temps_a)):
        outfile_a.write("{0}\t{1}\n".format(ant_temps_a[ind], final_vals_a[ind]))
        outfile_b.write("{0}\t{1}\n".format(ant_temps_b[ind], final_vals_b[ind]))

    outfile_a.close()
    outfile_b.close()

    assessCorrs(linespace_locs, lst_labels)'''
    

if __name__ == "__main__":
    main()
