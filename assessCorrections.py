#assessCorrections.py

""" This program will be used to plot uncorrected and temperature-
    adjusted antenna temperatures on top of each other to show 
    improvements (or lack thereof) in agreement.

    Ellie White 25 Aug. 2019

"""

import pylab as plt
from astropy.time import Time
import numpy as np
import numpy.polynomial.polynomial as poly

def assessCorrs(llocs, lstlabs):
    ifname1 = 'outfile_1a.txt'
    ifname2 = 'outfile_2a.txt'
    ifname3 = 'outfile_3a.txt'
    ifname4 = 'outfile_4a.txt'
    ifname5 = 'outfile_5a.txt'
    ifname6 = 'outfile_6a.txt'

    infile1 = open(ifname1, 'r')
    infile2 = open(ifname2, 'r')
    infile3 = open(ifname3, 'r')
    infile4 = open(ifname4, 'r')
    infile5 = open(ifname5, 'r')
    infile6 = open(ifname6, 'r')

    atemps1 = []
    atemps2 = []
    atemps3 = []
    atemps4 = []
    atemps5 = []
    atemps6 = []

    ctemps1 = []
    ctemps2 = []
    ctemps3 = []
    ctemps4 = []
    ctemps5 = []
    ctemps6 = []

    lst1 = []
    lst2 = []
    lst3 = []
    lst4 = []

    for line1 in infile1.readlines():
        line1list = line1.split()
        atemps1.append(float(line1list[0]))
        ctemps1.append(float(line1list[1])) 
        #lst1.append(float(line1list[2])) 

    for line2 in infile2.readlines():
        line2list = line2.split()
        atemps2.append(float(line2list[0]))
        ctemps2.append(float(line2list[1])) 
        #lst2.append(float(line2list[2])) 

    for line3 in infile3.readlines():
        line3list = line3.split()
        atemps3.append(float(line3list[0]))
        ctemps3.append(float(line3list[1])) 
        #lst3.append(float(line3list[2])) 

    for line4 in infile4.readlines():
        line4list = line4.split()
        atemps4.append(float(line4list[0]))
        ctemps4.append(float(line4list[1])) 
        #lst4.append(float(line4list[2])) 

    for line5 in infile5.readlines():
        line5list = line5.split()
        atemps5.append(float(line5list[0]))
        ctemps5.append(float(line5list[1])) 
        #lst4.append(float(line4list[2]))

    for line6 in infile6.readlines():
        line6list = line6.split()
        atemps6.append(float(line6list[0]))
        ctemps6.append(float(line6list[1])) 
     

    atemps1 = atemps1[20:]+atemps1[0:20]
    atemps3 = atemps3[152:]+atemps3[0:152]
    atemps4 = atemps4[130:]+atemps4[0:130]
    atemps5 = atemps5[100:]+atemps5[0:100]
    atemps6 = atemps6[-15:] + atemps6[:-15]

    ctemps1 = ctemps1[20:]+ctemps1[0:20]
    ctemps3 = ctemps3[152:]+ctemps3[0:152]
    ctemps4 = ctemps4[130:]+ctemps4[0:130]
    ctemps5 = ctemps5[100:]+ctemps5[0:100]
    ctemps6 = ctemps6[-15:]+ctemps6[:-15]

    tlinespace = np.linspace(0, (len(atemps1) - 1), len(atemps1))

    linespace_locs = []
    for l in range(15):
        lindex = int((len(atemps2)/15)*l)
        linespace_locs.append(tlinespace[lindex])

    linespace_locs.append(tlinespace[-1])

    plt.plot(atemps1, label='AUG 13')
    plt.plot(atemps2, label='AUG 16')
    plt.plot(atemps3, label='AUG 22')
    plt.plot(atemps4, label='AUG 24')
    plt.plot(atemps5, label='AUG 25')
    plt.plot(atemps6, label='AUG 28')
    plt.legend()
    plt.subplots_adjust(bottom=0.32)
    plt.xticks(linespace_locs, lstlabs, rotation=90)
    plt.xlabel("Local Sidereal Time")
    plt.ylabel("Antenna Temperature [K]")
    plt.title("Uncorrected timeseries graphs")
    #plt.xticks(linespace_locs, lst_labels, rotation=90)
    plt.show()

    plt.plot(ctemps1, label='AUG 13')
    plt.plot(ctemps2, label='AUG 16')
    plt.plot(ctemps3, label='AUG 22')
    plt.plot(ctemps4, label='AUG 24')
    plt.plot(ctemps5, label='AUG 25')
    plt.plot(ctemps6, label='AUG 28')
    plt.legend()
    plt.xticks(linespace_locs, lstlabs, rotation=90)
    plt.subplots_adjust(bottom=0.32)
    plt.xlabel("Local Sidereal Time")
    plt.ylabel("Antenna Temperature [K]")
    plt.title("Corrected timeseries graphs")
    plt.show()

    ''' Statistical comparison of uncorrected vs. corrected plots '''
    '''si = 220
    fi = 230

    uncorr_slice = np.array(atemps1[si:fi]) #sum_uncorr5 / 6.0
    ulinespace = np.linspace(0, (uncorr_slice.size - 1), uncorr_slice.size)

    uncorr_coeffs = poly.polyfit(ulinespace, uncorr_slice, 1)
    uncorr_fit = ulinespace*uncorr_coeffs[1] + uncorr_coeffs[0]

    diff_uc_1 = np.array(atemps1[si:fi]) - uncorr_fit
    diff_uc_2 = np.array(atemps2[si:fi]) - uncorr_fit
    diff_uc_3 = np.array(atemps3[si:fi]) - uncorr_fit
    diff_uc_4 = np.array(atemps4[si:fi]) - uncorr_fit
    diff_uc_5 = np.array(atemps5[si:fi]) - uncorr_fit
    diff_uc_6 = np.array(atemps6[si:fi]) - uncorr_fit

    uncorr = np.concatenate((diff_uc_1, diff_uc_2, diff_uc_3, diff_uc_4, diff_uc_5, diff_uc_6))

    stderr = np.std(uncorr)
    print(stderr)

    corr_slice = np.array(ctemps1[si:fi]) #sum_uncorr5 / 6.0

    corr_coeffs = poly.polyfit(ulinespace, corr_slice, 1)
    corr_fit = ulinespace*corr_coeffs[1] + corr_coeffs[0]

    diff_c_1 = np.array(ctemps1[si:fi]) - corr_fit
    diff_c_2 = np.array(ctemps2[si:fi]) - corr_fit
    diff_c_3 = np.array(ctemps3[si:fi]) - corr_fit
    diff_c_4 = np.array(ctemps4[si:fi]) - corr_fit
    diff_c_5 = np.array(ctemps5[si:fi]) - corr_fit
    diff_c_6 = np.array(ctemps6[si:fi]) - corr_fit

    corr = np.concatenate((diff_c_1, diff_c_2, diff_c_3, diff_c_4, diff_c_5, diff_c_6))

    stderrc = np.std(corr)
    print(stderrc)

    #plt.plot(uncorr_fit)
    #plt.plot(uncorr_slice)
    #plt.plot(corr_slice)
    #plt.plot(corr_fit)
    #plt.plot(rad_diff)
    plt.plot(uncorr_avg)
    plt.plot(corr_avg)
    plt.show()
    
    plt.hist(uncorr)
    plt.show()

    plt.hist(corr)
    plt.show()'''


    # find error bar for uncorrected data
    alinespace = np.linspace(0, (len(atemps1) - 1), len(atemps1))

    uncorr_avg = (np.array(atemps1) + np.array(atemps2) + np.array(atemps3) + np.array(atemps4) + np.array(atemps5) + np.array(atemps6))/6.0

    a1_coeffs = poly.polyfit(alinespace, uncorr_avg, 9)
    uc_fit = (alinespace**9)*a1_coeffs[9] + (alinespace**8)*a1_coeffs[8] + (alinespace**7)*a1_coeffs[7] + (alinespace**6)*a1_coeffs[6] + (alinespace**5)*a1_coeffs[5] + (alinespace**4)*a1_coeffs[4] + (alinespace**3)*a1_coeffs[3] + (alinespace**2)*a1_coeffs[2] + (alinespace)*a1_coeffs[1] + a1_coeffs[0]

    diff_uc_1 = np.array(atemps1) - uc_fit
    diff_uc_2 = np.array(atemps2) - uc_fit
    diff_uc_3 = np.array(atemps3) - uc_fit
    diff_uc_4 = np.array(atemps4) - uc_fit
    diff_uc_5 = np.array(atemps5) - uc_fit
    diff_uc_6 = np.array(atemps6) - uc_fit

    uncorr = np.concatenate((diff_uc_1, diff_uc_2, diff_uc_3, diff_uc_4, diff_uc_5, diff_uc_6))

    stderr = np.std(uncorr)
    print(stderr)

    # find error bar for corrected data
    corr_avg = (np.array(ctemps1) + np.array(ctemps2) + np.array(ctemps3) + np.array(ctemps4) + np.array(ctemps5) + np.array(ctemps6))/6.0

    c_coeffs = poly.polyfit(alinespace, corr_avg, 9)
    c_fit = (alinespace**9)*c_coeffs[9] + (alinespace**8)*c_coeffs[8] + (alinespace**7)*c_coeffs[7] + (alinespace**6)*c_coeffs[6] + (alinespace**5)*c_coeffs[5] + (alinespace**4)*c_coeffs[4] + (alinespace**3)*c_coeffs[3] + (alinespace**2)*c_coeffs[2] + (alinespace)*c_coeffs[1] + c_coeffs[0]

    diff_c_1 = np.array(ctemps1) - c_fit
    diff_c_2 = np.array(ctemps2) - c_fit
    diff_c_3 = np.array(ctemps3) - c_fit
    diff_c_4 = np.array(ctemps4) - c_fit
    diff_c_5 = np.array(ctemps5) - c_fit
    diff_c_6 = np.array(ctemps6) - c_fit

    corr = np.concatenate((diff_c_1, diff_c_2, diff_c_3, diff_c_4, diff_c_5, diff_c_6))

    stderrc = np.std(corr)
    print(stderrc)

if __name__ == "__main__":
    main()

