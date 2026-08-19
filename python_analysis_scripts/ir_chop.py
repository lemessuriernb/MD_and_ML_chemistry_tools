import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import scipy.fft as ff
import scipy.signal as sig



def dipole_ir(dip_center, time_span_f, time_span, length):
    pad1 = int(2**np.ceil(np.log2(time_span_f)))
    added_zeros1 = pad1 - time_span_f
    pad2 = int(2**np.ceil(np.log2(time_span)))
    added_zeros2 = pad2 - time_span
    n_avgs = int(length/time_span_f)

    m_ir_all = np.zeros((2*pad2,n_avgs))

    for i in range(n_avgs):
        m_pad = np.concatenate((dip_center[i*time_span_f:(i*time_span_f+time_span_f),:], np.zeros((added_zeros1, 3))), axis=0)
        m_corr = sig.correlate(m_pad, m_pad, 'full', 'fft')
        m_unpad = m_corr[pad1-1:, 2:]
        m_tot = np.sum(m_unpad, axis=1)
        m_tot_s = m_tot / np.flip(range(1, pad1 + 1))
        #plt.plot(m_tot_s)
        #plt.show()
    
        m_tot_up = m_tot_s[0:time_span]
        hann = np.hanning(2*time_span)
        m_hann = m_tot_up *hann[time_span:]
        m_hann_rep = np.concatenate((m_hann[:], np.zeros((added_zeros2))), axis=0)
        #plt.plot(m_hann_rep)
        #plt.show()

        m_double = np.hstack((m_hann_rep,np.zeros(1), np.flip(m_hann_rep)))
        m_double = m_double[:-1]
        m_ir_all[:,i] = np.abs(ff.fft(m_double))
        m_freq = ff.fftfreq(len(m_double), 40E-16)/2.998E10

    return m_ir_all, m_freq
    
def dipole_derivative(dipole, length): #This function calculates the time derivative of the dipole moment using a central difference method. 
    dip_grad = np.zeros((length,3))
    for i in range(3):
        dip_grad[:,i] = np.gradient(dipole[:,i], edge_order=2)
    dip_grad /= 4.0 # 4 fs is the time between each data point in the input file.

    return dip_grad

def dipole(dipole, length): #This function centers the dipole moment around its mean value.
    avg = np.mean(dipole,axis=0)
    print(avg)
    dip_centered = np.zeros((length,3))
    dip_centered[:,0] = dipole[:,0] - avg[0]
    dip_centered[:,1] = dipole[:,1] - avg[1]
    dip_centered[:,2] = dipole[:,2] - avg[2]

    return dip_centered

def plot_all(ir_all, freqs):
    first = 100
    last = int(0.5*len(ir_all[:,0]))

    avg_ir = np.mean(ir_all[:,:],axis=1)
    ir_std = np.std(ir_all[:,:],axis=1)
    plt.plot(freqs[first:last:1], avg_ir[first:last]) 
    plt.fill_between(freqs[first:last:1], avg_ir[first:last] - ir_std[first:last], avg_ir[first:last] + ir_std[first:last], alpha=0.5)

    #plt.savefig('/path/to/your/plot.eps', format='eps', dpi=600)
    plt.show()
    return avg_ir


def save_data(avg_ir, freqs, typec):
    np.savetxt(typec+'path/to/your/output_ir_data.txt', avg_ir[:], fmt='%.5e', delimiter=',', newline='\n', header='', footer='', comments='# ', encoding=None)
    np.savetxt(typec+'path/to/your/freq_data.txt', freqs[0:int(0.5*len(avg_ir))], fmt='%.5e', delimiter=',', newline='\n', header='', footer='', comments='# ', encoding=None)




if __name__ == '__main__':
    mac = '/Users/xxx/Documents/' #File path headers for managing file paths on different operating systems.
    pc = '/mnt/c/Users/xxx/Documents/'
    dip = np.loadtxt(pc+'path/to/your/dipole_data.txt')#, skiprows=2, usecols=(1,2,3), max_rows=2500000)

    dip_len = len(dip)

    dip_c = dipole(dip, dip_len)

    #dip_d = dipole_derivative(dip, dip_len)

    d_ir_all, d_freqs = dipole_ir(dip_c, 125000, 2500, dip_len)

    avg_d_ir = plot_all(d_ir_all, d_freqs) 

    save_data(avg_d_ir, d_freqs, pc)
