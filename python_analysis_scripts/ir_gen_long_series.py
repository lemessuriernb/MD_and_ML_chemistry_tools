import numpy as np
#import pandas as pd 
import scipy.fft as ff
import scipy.signal as sig
import matplotlib.pyplot as plt

# This file calculates the IR spectrum from a dipole derivative file. 
# It is designed to handle long time series data by partitioning the data into smaller segments, calculating the correlation for each segment, and then averaging the results.
# If you want to save the data, uncomment the save_data() function call in the main() function and provide the appropriate file paths.
# The other IR code file (ir_chop.py) does not partition the data, and instead calculates the correlation for the entire time series at once. The dipole autocorrelation is truncated before the IR spectrum is calculated.
def main():
    global n_avgs, time_span

    n_avgs = 2
    time_span = 500000
    file_name = '/path/to/your/datafile.txt'  # Replace with your actual dipole derivative file path

    load_data(file_name)

    partition_dipoles()

    calc_padding()

    big_loop_setup(dipole)

    avg_ir()

    #save_data(m_ir_avg,m_freq)

    plot_ir(m_ir_avg, m_freq)


def load_data(filename):
    global dip_data
    dip_data = np.loadtxt(filename, skiprows=2, usecols=(1, 2, 3))

def partition_dipoles():
    global dipole
    dipole = np.zeros((time_span, 3, n_avgs))
    for i in range(n_avgs):
        #dipole[:, :, i] = dip_data[(time_span+25000)*i:((time_span+25000)*i+time_span), :]
        dipole[:, :, i] = dip_data[(time_span)*i:((time_span)*i+time_span), :]

def big_loop_setup(dipole):
    global m_freq
    for j in range(n_avgs):
        m_pad = np.zeros((pad, 3))
        m_corr = np.zeros((2*pad - 1, 5))
        m_unpad = np.zeros((pad, 3))

        m_pad = np.concatenate((dipole[:, :, j], np.zeros((added_zeros, 3))), axis=0) #pad the data with zeros

        calc_correlation(m_pad, m_corr, m_unpad)

        windowing(m_tot_s)

        fft_ir(m_hann, j)

def calc_padding():
    global pad, added_zeros, m_ir
    pad = int(2**np.ceil(np.log2(time_span)))
    added_zeros = pad - time_span
    m_ir = np.zeros((2*pad, n_avgs))

def calc_correlation(padded, correlated, unpadded):
    global m_tot_s
    correlated = sig.correlate(padded, padded, 'full', 'fft')
    unpadded = correlated[pad-1:, 2:]
    m_tot = np.sum(unpadded, axis=1)
    m_tot_s = m_tot / np.flip(range(1, pad + 1))
    
def windowing(m_total):
    global m_hann
    hann = np.hanning(2*len(m_total))
    m_hann = m_total*hann[len(m_total):]

def fft_ir(windowed, index):
    global m_freq
    m_double = np.hstack((windowed,np.zeros(1), np.flip(windowed)))
    m_double = m_double[:-1]
    m_ir[:, index] = np.real(ff.fft(m_double))
    m_freq = ff.fftfreq(len(m_double), 2E-16)/2.998E10 #in units of cm^-1 (wavemumber)
    
    return m_freq

def avg_ir():
    global m_ir_avg
    m_ir_avg = np.mean(m_ir, axis=1)

def save_data(all_ir, avg_ir, frequencies):
    np.savetxt('/path/to/your/output/iravg.txt', avg_ir[0:7000], fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
    np.savetxt('/path/to/your/output/irfreq.txt', frequencies[0:7000], fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)
    np.savetxt('/path/to/your/output/irall.txt', all_ir[0:7000,:], fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ', encoding=None)

def plot_ir(avg_ir, frequencies):
    plt.plot(frequencies[:], avg_ir[:])
    plt.show()

if __name__ == '__main__':
    main()