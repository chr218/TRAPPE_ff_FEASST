# -*- coding: utf-8 -*-
"""
Spyder Editor

The purpose of this code is to convert the parameters taken from Siepmann's
website (TRAPPE-UA) into FEASST readible input files. The input files are 
restrcited to rigid molecules.

Author: Christopher Rzepa
email: cvr5246@gmail.com
"""
import numpy as np
import time
import sys

from scipy import constants
from ase.data.pubchem import pubchem_atoms_search
from ase.geometry.analysis import Analysis

'''
Constants
'''

NA = constants.physical_constants['Avogadro constant']
kB = constants.physical_constants['Boltzmann constant']
#print(constants.calorie)

'''
Functions & Subroutines
'''

def get_psuedoatoms(atoms_obj,TRAP_atms,TRAP_bonds):
    
    ana_H = Analysis(atoms_obj)
    
    #Collect bond distances from ASE atoms.-Find easier way to do this!
    ASE_bond_dist = []
    ASE_bond_atms = []
    atm_bonds_list = ana_H.unique_bonds[0]
    for index, atm_bonds in enumerate(atm_bonds_list):
        for neighbor in atm_bonds:
            dist = np.linalg.norm(atoms_obj[index].position-atoms_obj[neighbor].position)
            ASE_bond_dist.append(dist)
            ASE_bond_atms.append([index,neighbor])
  
    H_index = []
    ASE_bond_atms_cp = ASE_bond_atms.copy()
    ASE_bond_dist_cp = ASE_bond_dist.copy()
    atoms_obj_cp = atoms_obj.copy()
    #If "H" is in psuedo potentials, we must include it and collect it's index.
    if 'H' in TRAP_atms[:,1]:    
        for bond in TRAP_bonds:
            if 'H' in bond[2].split('-'):
                H_bond_dist = float(bond[3])
                indx, value = find_nearest(ASE_bond_dist_cp,H_bond_dist)
                for i in (ASE_bond_atms_cp[indx]):
                    if atoms_obj_cp[i].symbol == 'H':
                        H_index.append(atoms_obj_cp[i].index)
                        
                        ASE_bond_atms_cp.pop(indx)
                        ASE_bond_dist_cp.pop(indx)

    #Determine hybridization. "atoms_C" does not strictly apply to carbon atoms.       
    atoms_H = atoms_obj.copy()
    del atoms_H[[atom.index for atom in atoms_H if not atom.symbol=='H']]
    
    atoms_C = atoms_obj.copy()
    del atoms_C[[atom.index for atom in atoms_C if atom.symbol=='H' and atom.index not in H_index]]   

    for i,C in enumerate(atoms_C):
        hybridization = 0
        if C.symbol == 'C':
            for j,H in enumerate(atoms_H):
                if np.linalg.norm(C.position-H.position) < 1.2:
                    hybridization += 1
            atoms_C[i].mass += hybridization*(atoms_H[0].mass)
            atoms_C[i].tag = hybridization

    return atoms_C

def get_TRAPPE_params(fname):
    main_full = []
    main = []
    f = open(fname,"r")
    for line in f:
        stripped_line = line.strip()
        if stripped_line == '': #Skip new lines
            continue
        else:
            main.append(stripped_line.split(','))
            main_full.append(stripped_line.split(' '))
    f.close()
    
    TRAPPE_name = main_full[0]
    
    TRAPPE=[[],[]] #List constaining TRAPPE atoms [0] and TRAPPE bonds [1] info.
    
    for i, line in enumerate(main):
        if 'stretch' in line:
            j = i
            while '#' not in main[j+1]:
                TRAPPE[1].append(main[j+1])
                j += 1
        if '(pseudo)atom' in line:
            j = i
            while '#' not in main[j+1]:
                TRAPPE[0].append(main[j+1])
                j += 1    
    
    #Convert TRAPPE list into array for psuedoatoms and array for bonds:
    TRAPPE_atoms_arr = (np.asarray(TRAPPE[0]).reshape((len(TRAPPE[0]),len(TRAPPE[0][0]))))
    TRAPPE_bonds_arr = (np.asarray(TRAPPE[1]).reshape((len(TRAPPE[1]),len(TRAPPE[1][0]))))
    
    return TRAPPE_atoms_arr, TRAPPE_bonds_arr, TRAPPE_name
  
#Used within "Order_atoms_wrt_TRAPPE" to find closest matching bond distance.
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

#Orders the ASE atom indieces w.r.t TRAPPE atom indices.
def Order_atoms_wrt_TRAPPE(ASE_atms, ASE_geom, TRAP_bonds,TRAP_atms):
    
    #Create list of psuedoatom strings.
    ASE_hybrid_list = []
    for i, atm in enumerate(ASE_atms):
        if atm.symbol =='C':              
            if atm.tag == 0:
                ASE_hybrid_list.append('C')
            elif atm.tag == 1:
                ASE_hybrid_list.append('CH')
            else:
                ASE_hybrid_list.append('CH' + str(atm.tag))
        elif atm.symbol != 'C':
            ASE_hybrid_list.append(atm.symbol)
    
    #Collect bond distances from ASE atoms.-Find easier way to do this!
    ASE_bond_dist_list = []
    ASE_bond_atms = []
    ASE_bond_hybrids = []
    atm_bonds_list = ana.unique_bonds[0]
    for index, atm_bonds in enumerate(atm_bonds_list):
        for neighbor in atm_bonds:
            dist = np.linalg.norm(psuedo_atoms[index].position-psuedo_atoms[neighbor].position)
            ASE_bond_dist_list.append(dist)
            ASE_bond_atms.append([index,neighbor])
            ASE_bond_hybrids.append([ASE_hybrid_list[index],ASE_hybrid_list[neighbor]])
    
    #Collect bond distances from TRAPPE params.
    #Strip quotation marks and convert to int.
    TRAP_bond_dist_list = []
    TRAP_bond_atms = []
    TRAP_bond_hybrids = []
    TRAP_bond_charges = []
    
    for i,bond in enumerate(TRAP_bonds):
        TRAP_bond_dist = float(bond[-1])
        TRAP_bond_dist_list.append(TRAP_bond_dist)
        TRAP_bond_desc = bond[1].split('-')
        TRAP_bond_atms.append(TRAP_bond_desc)
        TRAP_bond_hybrids.append([])
        TRAP_bond_charges.append([])

    #Clean TRAPPE psuedo atom params
    for i,bond2 in enumerate(TRAP_bond_atms):
        #Clean atom indices wihtin bond type.
        for j,el in enumerate(bond2):
            TRAP_bond_atms[i][j] = el.strip('\'\"')
            #Collect hybridization of psuedoatoms corresponding to bond.
            for k,line in enumerate(TRAP_atms):
                if TRAP_bond_atms[i][j].strip() == line[0]:
                    
                    TRAP_bond_hybrids[i].append(line[1])
                    TRAP_bond_charges[i].append(line[-1])
    
    skip_ASE_indices = []    
    skip_TRAPPE_indices = []
    ASE_new_positions = []
    #Match TRAPPE bonds/indices w/ ASE bonds/indices.
    for i, dist in enumerate(TRAP_bond_dist_list):
        
        #Find index and value of ASE bond length which is closest to TRAPPE bond length.
        ASE_bond_dist_list_tmp = ASE_bond_dist_list.copy()        
        indx, value = find_nearest(ASE_bond_dist_list,dist)
        
        #Ensure hybrids in ASE bond match TRAPPE bond.
        timeout_loop = time.time() + 60
        while not set(TRAP_bond_hybrids[i]) == set(ASE_bond_hybrids[indx]):
            ASE_bond_dist_list_tmp[indx] = 0
            indx, value = find_nearest(ASE_bond_dist_list_tmp,dist)            
            if time.time() > timeout_loop:
                sys.exit('Error: ASE Hybrids could not be matched with TRAPPE Hybrids. Check SMILES structure &/or adjust Hydrogen bond criterion.')
                  
        #Change ASE index to TRAPPE index by matching psuedoatom hybridization.
        for j,TRAPE_hybrid in enumerate(TRAP_bond_hybrids[i]):
            #print(TRAP_bond_atms[i][j])
            
            for k, ASE_hybrid in enumerate(ASE_bond_hybrids[indx]):
                
                #Ensure that closest bond distance assumption is valid.
                if TRAPE_hybrid not in ASE_bond_hybrids[indx]:
                    sys.exit('Error: Hybrid not found in ASE bond. Check distances.')
                
                elif TRAPE_hybrid == ASE_hybrid:
                    #Collect corresponding ASE info.
                    #print(ASE_bond_atms[indx],ASE_bond_dist_list[indx],ASE_bond_hybrids[indx])
                    #Collect corresponding TRAPPE info.
                    #print(TRAP_bond_atms[i],TRAP_bond_dist_list[i],TRAP_bond_hybrids[i])
                    #print(TRAP_bond_atms[i][j])
                    

                    for z,atm in enumerate(ASE_atms):
                        if atm.index == int(ASE_bond_atms[indx][k]) and atm.index not in skip_ASE_indices and int(TRAP_bond_atms[i][j]) not in skip_TRAPPE_indices:
                            ASE_new_positions.append([int(TRAP_bond_atms[i][j]),atm.mass,TRAP_bond_charges[i][j],atm.position])
                            skip_ASE_indices.append(atm.index)
                            skip_TRAPPE_indices.append(int(TRAP_bond_atms[i][j]))

                            
        ASE_bond_dist_list[indx] = False
        ASE_bond_atms[indx] = False
        ASE_bond_hybrids[indx] = False
                 
    return ASE_new_positions


'''
Main Program Begins Here
'''
#TRAPPE 
#Collect TRAPPE input params
TRAPPE_fname =  'trappe_parameters_12.txt'
TRAPPE_atoms, TRAPPE_bonds, TRAPPE_name = get_TRAPPE_params(TRAPPE_fname)

#ASE
#Convert Smiles string into Atoms object.
SMILES = 'COC'
atoms = pubchem_atoms_search(smiles=SMILES)

#Determine psuedo-atoms of molecule and remove hydrogens.
psuedo_atoms = get_psuedoatoms(atoms,TRAPPE_atoms,TRAPPE_bonds)

ana = Analysis(psuedo_atoms) #Create geometry analysis object from Atoms object.

organized_ASE_atoms = Order_atoms_wrt_TRAPPE(psuedo_atoms,ana,TRAPPE_bonds,TRAPPE_atoms)

#Center the psuedo_atoms. FEASST requires atom[0] to be located within origin!
for atm in organized_ASE_atoms:
    if atm[0] == 1:
        center_about = np.copy(atm[-1])
        break
        
for atm in organized_ASE_atoms:
    atm[-1] -= center_about

'''
File Output
'''
f = open('FEASST_'+TRAPPE_fname.split('_')[-1].split('.')[0]+'.in',"w")

#Heading
f.write('# FEASST data file \n# %s \n# %s \n# %s \n' % (TRAPPE_name[0],TRAPPE_fname,SMILES)) 

#Preface
f.write('\n%s atoms \n%s bonds \n\n%s atom types \n%s bond types \n' % 
        (len(TRAPPE_atoms[:,0]),len(TRAPPE_bonds[:,0]), 
         len(TRAPPE_atoms[:,1]),len(TRAPPE_bonds[:,2])))
f.write('\n%s %s xlo xhi \n%s %s ylo yhi \n%s %s zlo zhi \n' % 
        (np.amin(psuedo_atoms.positions[:,0]),np.amax(psuedo_atoms.positions[:,0]),
         np.amin(psuedo_atoms.positions[:,1]),np.amax(psuedo_atoms.positions[:,1]),
         np.amin(psuedo_atoms.positions[:,2]),np.amax(psuedo_atoms.positions[:,2])))

#Masses
f.write('\nMasses\n\n')
for i,atm in enumerate(organized_ASE_atoms):
    f.write('%s %s \n' % (atm[0],atm[1]))    

#Pair Coeffs
f.write('\nPair Coeffs\n\n')
for i,line in enumerate(TRAPPE_atoms):
    f.write('%s %s %s \n' % (TRAPPE_atoms[i,0],float(TRAPPE_atoms[i,3])*kB[0]*NA[0]/1000,TRAPPE_atoms[i,4]))

#Bond Coeffs
f.write('\nBond Coeffs\n\n')
for i,line in enumerate(TRAPPE_bonds):
    f.write('%s -1 %s \n' % (TRAPPE_bonds[i,0],TRAPPE_bonds[i,-1]))
    
#Atoms
f.write('\nAtoms\n\n')
for i,atm in enumerate(organized_ASE_atoms):
    f.write('%s %s %s %s %s %s %s 0 0 0\n' % (atm[0],1,atm[0],atm[2],atm[3][0],atm[3][1],atm[3][2]))    

#Bonds
f.write('\nBonds\n\n')
for i,bond in enumerate(TRAPPE_bonds):
    f.write('%s %s %s %s\n' % (bond[0],bond[0],bond[1].strip('"\'').split('-')[0].strip(),bond[1].strip('"\'').split('-')[1].strip()))
      
f.close()