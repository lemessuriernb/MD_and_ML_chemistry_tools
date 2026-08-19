import numpy as np
import pandas as pd
import os

def main():
    global n_atoms, n_full_atoms, file_path, file_pos, file_vel, file_frc, file_charge, output_starter
    n_atoms = 29
    n_full_atoms = 528
    file_path = 'C:/Users/lemes/OneDrive - UCB-O365/Documents/2024exps/cp2k/'
    file_pos = 'PYR_1W_QMMM-pos-10ps.xyz'
    file_vel = 'PYR_1W_QMMM-vel-10ps.xyz'
    file_frc = 'PYR_1W_QMMM-frc-1.xyz'
    file_charge = 'hirshfeld_py1w-10ps.txt'
    output_starter = 'pyr_1w_qm_'

    sort_positions()

    sort_velocities()

    sort_forces()

    sort_charges()



def sort_positions():
    valid_pos = lambda x: 0 <= (x + 1) % (n_full_atoms+2) <= 2 or (n_atoms+3) <= (x + 1) % (n_full_atoms+2) <= (n_full_atoms+2) 
    dpos = pd.read_csv(file_path+file_pos, skiprows=valid_pos, delim_whitespace=True, names=["id", "x", "y", "z"])

    valid_comment = lambda x: 0 <= (x + 1) % (n_full_atoms+2) <= 1 or 3 <= (x + 1) % (n_full_atoms+2) <= (n_full_atoms+2) 
    top_info = pd.read_csv(file_path+file_pos, skiprows=valid_comment, delim_whitespace=True, names=['a','b', "step", 'c', 'd', "time", 'e', 'f', "energy"])

    for i in range(len(top_info)):
        comment_line = [f"{n_atoms} \n", f"{top_info.iloc[[i]].to_string(header=False, index=False)} \n"]
        with open(file_path + output_starter + 'pos.xyz', 'a') as file:
            file.writelines(comment_line)
            file.write(dpos.iloc[i*n_atoms:(i*n_atoms+n_atoms),:].to_string(header=False, index=False))
            file.write("\n")

def sort_velocities():
    valid_vel = lambda x: 0 <= (x + 1) % (n_full_atoms+2) <= 2 or (n_atoms+3) <= (x + 1) % (n_full_atoms+2) <= (n_full_atoms+2) 
    dvel = pd.read_csv(file_path+file_vel, skiprows=valid_vel, delim_whitespace=True, names=["id", "x", "y", "z"])

    valid_comment = lambda x: 0 <= (x + 1) % (n_full_atoms+2) <= 1 or 3 <= (x + 1) % (n_full_atoms+2) <= (n_full_atoms+2) 
    top_info_v = pd.read_csv(file_path+file_vel, skiprows=valid_comment, delim_whitespace=True, names=['a','b', "step", 'c', 'd', "time", 'e', 'f', "energy"])

    for i in range(len(top_info_v)):
        comment_line = [f"{n_atoms} \n", f"{top_info_v.iloc[[i]].to_string(header=False, index=False)} \n"]
        with open(file_path + output_starter + 'vel.xyz', 'a') as file:
            file.writelines(comment_line)
            file.write(dvel.iloc[i*n_atoms:(i*n_atoms+n_atoms),:].to_string(header=False, index=False))
            file.write("\n")

def sort_forces():
    valid_frc = lambda x: 0 <= (x + 1) % (n_full_atoms+2) <= 2 or (n_atoms+3) <= (x + 1) % (n_full_atoms+2) <= (n_full_atoms+2) 
    dfrc = pd.read_csv(file_path+file_frc, skiprows=valid_frc, delim_whitespace=True, names=["id", "x", "y", "z"])

    valid_comment = lambda x: 0 <= (x + 1) % (n_full_atoms+2) <= 1 or 3 <= (x + 1) % (n_full_atoms+2) <= (n_full_atoms+2) 
    top_info_f = pd.read_csv(file_path+file_frc, skiprows=valid_comment, delim_whitespace=True, names=['a','b', "step", 'c', 'd', "time", 'e', 'f', "energy"])

    for i in range(len(top_info_f)):
        comment_line = [f"{n_atoms} \n", f"{top_info_f.iloc[[i]].to_string(header=False, index=False)} \n"]
        with open(file_path + output_starter + 'frc.xyz', 'a') as file:
            file.writelines(comment_line)
            file.write(dfrc.iloc[i*n_atoms:(i*n_atoms+n_atoms),:].to_string(header=False, index=False))
            file.write("\n")

def sort_charges():
    valid_charge = lambda x: 0 <= (x + 1) % (n_atoms+8) <= 5 or (n_atoms+6) <= (x + 1) % (n_atoms+8) <= (n_atoms+8)  
    charge = pd.read_csv(file_path+file_charge, skiprows=valid_charge, delim_whitespace=True, names=["Atom", "Element", "Kind", "Ref Charge", "pop 1", "pop 2", "Spin moment", "Net Charge"])

    for i in range(int(len(charge)/n_atoms)-1):
        comment_line = [f"{n_atoms}, step {i} \n"]
        with open(file_path+'pyr_1w_qm_charge.txt', 'a') as file:
            file.writelines(comment_line)
            
            if n_atoms == 29:
                file.write(charge.iloc[(i*n_atoms+1):(i*n_atoms+12),:].to_string(header=True, index=False))
                file.write("\n")
                file.write(charge.iloc[[(i*n_atoms+17)]].to_string(header=False, index=False))
                file.write("\n")
                file.write(charge.iloc[(i*n_atoms+12):(i*n_atoms+17),:].to_string(header=False, index=False))
                file.write("\n")
                file.write(charge.iloc[(i*n_atoms+18):(i*n_atoms+27),:].to_string(header=False, index=False))
                file.write("\n")
                file.write(charge.iloc[[(i*n_atoms)]].to_string(header=False, index=False))
                file.write("\n")
                file.write(charge.iloc[(i*n_atoms+27):(i*n_atoms+29),:].to_string(header=False, index=False))

            elif n_atoms == 32:
                file.write(charge.iloc[(i*n_atoms+2):(i*n_atoms+13),:].to_string(header=True, index=False))
                file.write("\n")
                file.write(charge.iloc[[(i*n_atoms+18)]].to_string(header=False, index=False))
                file.write("\n")
                file.write(charge.iloc[(i*n_atoms+13):(i*n_atoms+18),:].to_string(header=False, index=False))
                file.write("\n")
                file.write(charge.iloc[(i*n_atoms+19):(i*n_atoms+28),:].to_string(header=False, index=False))
                file.write("\n")
                file.write(charge.iloc[[(i*n_atoms)]].to_string(header=False, index=False))
                file.write("\n")
                file.write(charge.iloc[(i*n_atoms+28):(i*n_atoms+30),:].to_string(header=False, index=False))
                file.write("\n")
                file.write(charge.iloc[[(i*n_atoms+1)]].to_string(header=False, index=False))
                file.write("\n")
                file.write(charge.iloc[(i*n_atoms+30):(i*n_atoms+32),:].to_string(header=False, index=False))

            file.write("\n")


if __name__ == '__main__':
    main()