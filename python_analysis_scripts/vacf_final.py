import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.fft as ff
import scipy.signal as sig


def main():
    global n_atoms, filename, time_span, trajectories_max
    n_atoms = 256
    trajectories_max = 2000
    time_span = 2000
    filename = 'Documents/research/waters/cluster-results/spce_m_vxOx_5.txt'

    import_data()

    #summation_method()

    correlate_vacf()

    plot_vacfs()
    #np.savetxt('Documents/research/waters/cluster-results/vacf_spce_1.txt', vacf_tot_c_scaled, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ',encoding=None)


def import_data():
    global v_list
    valid_f = lambda x: 1 <= (x + 1) % (n_atoms + 9) <= 9
    df = pd.read_csv(filename, skiprows=valid_f, delim_whitespace=True, names=[
        "id", "type", "vx", "vy", "vz"
    ])  # "id","type","xu","yu","zu","vx","vy","vz","q"
    v_list = df.values


def summation_method():
    global VACFtot

    VACFtot = np.zeros(time_span)

    for Y in range(trajectories_max):  # loop over each new time zero to represent a new simulation
        line = Y * n_atoms
        VACFav = np.zeros(time_span)
        v1 = v_list[line:(line + n_atoms), :]
        v1 = v1[np.argsort(v1[:, 0])]  # sort indices

        if np.round(Y / 100) == Y / 100:
            print(Y)

        for Q in range(Y, Y + time_span):  # loop over later times
            q_line = Q * n_atoms
            v2 = v_list[q_line:(q_line + n_atoms), :]
            v2 = v2[np.argsort(v2[:, 0])]

            V = np.zeros(n_atoms)
            for k in range(n_atoms):
                v_atom = np.dot(v1[k, 5:8], v2[k, 5:8])
                # V = V + v_atom
                V[k] = v_atom

            VACFav[Q - Y] = np.average(V)  # average over all particles in the system

        VACFtot += VACFav  # average over all simulations
    VACFtot /= trajectories_max


def compute_correlation():
    global vacf_avg
    for j in range(n_atoms):
        v_pad2 = np.concatenate((vtot[j, :, :], np.zeros((added_zeros, 3))), axis=0)
        va2[j, :, :] = sig.correlate(v_pad2, v_pad2, 'full', 'fft')
    v_unpad = va2[:, pad - 1:-added_zeros, 2:]
    v_unpad_1d = np.sum(v_unpad, axis=2)
    vacf_avg = np.mean(v_unpad_1d, axis=0)


def organize_data():
    for i in range(k, k + time_span):
        line = i * n_atoms
        v1 = v_list[line:(line + n_atoms), :]
        v1 = v1[np.argsort(v1[:, 0])]

        vtot[:, i - k, :] = v1[:, 2:5]


def compute_padding():
    global pad, added_zeros
    pad = int(2 ** np.ceil(np.log2(time_span)))
    added_zeros = pad - time_span


def correlate_vacf():
    global k, vtot, vacf_avg, va2, vacf_tot_c_scaled
    compute_padding()
    vacf_tot_c = np.zeros(time_span)
    for k in range(trajectories_max):
        if np.round(k / 100) == k / 100:
            print(k)

        vtot = np.zeros((n_atoms, time_span, 3))
        vacf_avg = np.zeros(2 * pad - 1)
        va2 = np.zeros((n_atoms, 2 * pad - 1, 5))

        organize_data()

        compute_correlation()

        vacf_tot_c += vacf_avg
    vacf_tot_c /= trajectories_max
    vacf_tot_c_scaled = vacf_tot_c / np.flip(range(1, time_span + 1))


def fft_ifft_data():
    global va
    va = np.zeros((n_atoms, 2 * pad, 3))
    for j in range(n_atoms):
        v_corr = np.zeros((2 * pad, 3))
        v_fft = np.zeros((2 * pad, 3))
        v2tot = np.zeros((2 * pad, 3))
        v_pad = np.concatenate((vtot[j, :, :], np.zeros((added_zeros, 3))), axis=0)
        v2tot = np.concatenate((v_pad, np.zeros((1, 3)), np.flip(v_pad, axis=0)), axis=0)
        v2tot = v2tot[:-1, :]

        v_fft = np.real(ff.fft2(v2tot))

        for c in range(0, 3):
            v_corr[:, c] = v_fft[:, c] * v_fft[:, c]

        va[j, :, :] = np.real(ff.ifft2(v_corr))


def plot_vacfs():
    #plt.plot(VACFtot[:700], label='sum VACF')
    plt.plot(vacf_tot_c_scaled[:700], label='correlate VACF')
    plt.title('VACF for TIP4P/2005f Water Oxygens')
    plt.xlabel('Time (fs)')
    plt.ylabel('VACF (A^2/fs)')
    #plt.legend()
    plt.show()


if __name__ == '__main__':
    main()




