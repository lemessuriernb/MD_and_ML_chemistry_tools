import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import scipy.stats as st

def gaussian_bash_parsing(num_files, num_atoms, config_tag, output_tag, pos=True, forces=True, energies=True, dipoles=True, charges=True, pos_standard=True):
    gauss_file = 'gauss_' + config_tag + '_$i.log'
    chunk_start = "for i in $(seq 0 "+str(num_files-1)+"); do "

    if pos==True:
        chunk_pos = """awk '/Input orientation:/{{getline;getline;getline;getline;} for(i;i<="""+str(num_atoms-1)+""";i++) {getline;print $2 " " $4 " " $5 " " $6}}' """ + gauss_file + """ >> ../positions_""" + output_tag + """.txt; """
    elif pos==False:
        chunk_pos = ""

    if forces==True:
        chunk_forces = """awk '/Forces \(Hartrees/{{getline;getline;} for(i;i<="""+str(num_atoms-1)+""";i++) {getline;print $3 " " $4 " " $5}}' """ + gauss_file + """ >> ../forces_""" + output_tag + """.txt; """
    elif forces==False:
        chunk_forces = ""

    if energies==True:
        chunk_energies = """awk '/SCF Done/{print $5}' """ + gauss_file + """ >> ../energies_""" + output_tag + """.txt; """
    elif energies==False:
        chunk_energies = ""

    if dipoles==True:
        chunk_dipoles = """awk  '/Dipole moment/{getline;print $2 " " $4 " " $6}'  """ + gauss_file + """ >> ../dipole_""" + output_tag + """.txt; """
    elif dipoles==False:
        chunk_dipoles = ""

    if charges==True:
        chunk_charges = """awk '/Charges from ESP fit,/{{getline;getline;} for (i;i<="""+str(num_atoms-1)+""";i++) {getline;print $0}}' """ + gauss_file + """ >> ../charges_""" + output_tag + """.txt;  """
    elif charges==False:
        chunk_charges = ""

    if pos_standard==True:
        chunk_pos_standard = """awk '/Standard orientation:/{{getline;getline;getline;getline;} for(i;i<="""+str(num_atoms-1)+""";i++) {getline;print $2 " " $4 " " $5 " " $6}}' """ + gauss_file + """ >> ../positions_standard_""" + output_tag + """.txt; """
    elif pos_standard==False:
        chunk_pos_standard = ""

    chunk_end = "done"
    bash_command = chunk_start + chunk_pos + chunk_forces + chunk_energies + chunk_dipoles + chunk_charges + chunk_pos_standard + chunk_end
    print(bash_command)
    return bash_command

def load_forces_energy_data(num_files, num_atoms, output_tag, filepath=str):
    #hartree_kcalmol = 627.509
    bohr_angstrom = 1.88973 
    #debye_eA = 0.208194 
    hartree_eV = 27.211386
    
    dp = pd.read_csv(filepath + 'positions_'+output_tag+'.txt', sep=r"\s+", names=["id", "x", "y", "z"]) 
    #Convert the Atomic numbers into the atomic symbols:
    
    dp["id"] = dp["id"].map({6: 'C', 1: 'H', 8: 'O'})
    #dp["id"].replace(6,'C', inplace=True)
    #dp["id"].replace(1,'H', inplace=True)
    #dp["id"].replace(8,'O', inplace=True)

    forces = pd.read_csv(filepath + 'forces_'+output_tag+'.txt' , sep=r"\s+", names=["fx", "fy", "fz"])
    #Convert forces from Hartree/Bohr to eV/A:
    forces.fx *= hartree_eV*bohr_angstrom
    forces.fy *= hartree_eV*bohr_angstrom
    forces.fz *= hartree_eV*bohr_angstrom

    #Combine the positions and forces into 1 dataframe:
    coords_forces = dp.join(forces)
    print(coords_forces)

    if len(coords_forces)/(num_atoms) == num_files:
        print('Number of configurations loaded correctly!', int(len(coords_forces)/(num_atoms)), 'configurations.')

    energy = np.loadtxt(filepath + 'energies_'+output_tag+'.txt', max_rows=num_files)*hartree_eV

    return coords_forces, energy

def load_dipoles_data(num_files, num_atoms, output_tag, filepath=str):
    hartree_eV = 27.211386
    #debye_eA = 0.208194 

    coords = pd.read_csv(filepath + 'positions_standard_'+output_tag+'.txt', sep=r"\s+", names=["id", "x", "y", "z"]) 
    #Convert the Atomic numbers into the atomic symbols:
    
    coords["id"] = coords["id"].map({6: 'C', 1: 'H', 8: 'O'})
    #dp["id"].replace(6,'C', inplace=True)
    #dp["id"].replace(1,'H', inplace=True)
    #dp["id"].replace(8,'O', inplace=True)

    if len(coords)/(num_atoms) == num_files:
        print('Number of configurations loaded correctly!', int(len(coords)/(num_atoms)), 'configurations.')

    dipoles = np.loadtxt(filepath + 'dipole_'+output_tag+'.txt', max_rows=num_files, usecols=(0,1,2))

    energy = np.loadtxt(filepath + 'energies_'+output_tag+'.txt', max_rows=num_files)*hartree_eV

    return coords, dipoles, energy

def write_energy_forces_xyz(coords_forces_data, energy_data, output_filename=str, num_configs=int, num_atoms=int):
    with open(output_filename, 'a') as file:    
        for i in range(num_configs):
            comment_line = [f"{num_atoms} \n", f"Lattice=\"18.0 0.0 0.0 0.0 18.0 0.0 0.0 0.0 18.0\" Properties=species:S:1:pos:R:3:dft_forces:R:3 dft_energy={energy_data[i]} pbc=\"F F F\" \n"] 
            file.writelines(comment_line)
            file.write(coords_forces_data.iloc[i*num_atoms:(i*num_atoms+num_atoms),:].to_string(header=False, index=False))
            file.write("\n")
    print(f'Wrote {num_configs} configurations to {output_filename}')
    return

def write_dipoles_xyz(coords_data, dipoles_data, energy_data, output_filename=str, num_configs=int, num_atoms=int):
    with open(output_filename, 'a') as file:    
        for i in range(num_configs):
            comment_line = [f"{num_atoms} \n", f"Lattice=\"18.0 0.0 0.0 0.0 18.0 0.0 0.0 0.0 18.0\" Properties=species:S:1:pos:R:3 dft_energy={energy_data[i]} pbc=\"F F F\" debye_dipole=\"{dipoles_data[i,0]} {dipoles_data[i,1]} {dipoles_data[i,2]}\" \n"] 
            file.writelines(comment_line)
            file.write(coords_data.iloc[i*num_atoms:(i*num_atoms+num_atoms),:].to_string(header=False, index=False))
            file.write("\n")
    print(f'Wrote {num_configs} configurations to {output_filename}')
    return