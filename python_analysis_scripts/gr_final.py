from turtle import end_fill
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    global points, g, n_atoms, L, density, dr, filename, start, stop, step_skip
    points = 100
    #g = np.zeros((points, stop-start))
    n_atoms = 256  # 864 LJ 256 H2O
    L = 10  # 10.229 LJ 19.72 H2O
    header = 9
    Linv = 1 / L
    density = n_atoms / (L ** 3)
    dr = L / points
    filename = 'Documents/research/spce/dump_spce_oxy.lammpstrj'
    start = 1
    stop = 1000
    step_skip = 1
    g = np.zeros((points, stop-start))

    import_data()

    calc_g()

    plot_g(g, rx)

    integrate_g(g, rx)

def import_data():
    global rlist
    valid_f = lambda x: 1 <= (x + 1) % (n_atoms + 9) <= 9
    df = pd.read_csv(filename, skiprows=valid_f, delim_whitespace=True, names=[
        "id", "type", "xu", "yu", "zu", "vx", "vy", "vz", "q"
    ])
    rlist = df.values


def bin_particles(n_atoms, r, L, dr, h):
    for i in range(n_atoms):
        for j in range(i):
            d2 = 0.0
            rij = r[j, :] - r[i, :]
            rij -= L * np.round(rij / L)

            d2 = np.dot(rij, rij)
            ibin = int(np.sqrt(d2) / dr)
            h[ibin] += 2
    return h


def find_denominator(h, dr, density, n_atoms):
    for i in range(1, len(h)):
        rx[i] = i * dr
        h[i] /= 4 * np.pi * density * n_atoms * rx[i] * rx[i] * dr

    return h, rx


def calc_g():
    global g, rx, dr
    for t in range(start, stop, step_skip):
        tc = t - start + 1
        #g *= (tc - 1) / tc  # t / (t + 1)
        h = np.zeros(points)
        rx = np.zeros(points)

        if np.round(t / 100) == t / 100:
            print("Step =", t)

        l = t * n_atoms
        r = rlist[l:(l + n_atoms), 2:5]

        dr = L / points

        bin_particles(n_atoms, r, L, dr, h)

        find_denominator(h, dr, density, n_atoms)

        g[:,tc-1] = h

    return g, rx


def plot_g(g, rx):
    g_avg = np.mean(g, axis=1)
    g_std = np.std(g, axis=1)/(stop-start)

    plt.figure()
    plt.fill_between(rx, g_avg+g_std, g_avg-g_std, color='y')
    plt.plot(rx, g_avg)
    plt.show()
    #plt.savefig('gr_h2o_good.png', dpi=300, transparent=True)


def save_output(g):
    np.savetxt('gr_h2o_tip4p.txt', g, fmt='%.18e', delimiter=' ', newline='\n', header='', footer='', comments='# ',
               encoding=None)


def integrate_g(g, rx):
    g_avg = np.mean(g, axis=1)
    integ = 4 * np.pi * density * np.trapz(g_avg * rx * rx, dx=L / points)
    print("integral = ", integ)


if __name__ == '__main__':
    main()