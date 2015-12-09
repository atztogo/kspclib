#!/usr/bin/python
import sys
import logging
import numpy as np
import utils
import tetrahedron_method_interface

def main():
    # generate som simple bands based on the std. parabolic approximation

    # set testcell
    cell=np.array([[12.0,0.0,0.0],[0.0,12.0,0.0],[0.0,0.0,12.0]])

    # rec_basis (only cubic)
    rec_basis=2*math.pi/cell
                  
    # first, set kmesh (we can probably just use kgrid in this dir now, but need interface, fix tomorrow)
    kmesh = 

    # transform to cart, no need to save
    kmesh=dir_to_cart(kmesh,rec_basis)
    
    # effective mass
    effmass=np.array([1.0,1.0,1.0])
    
    # fetch energy
    energy=generate_energy(kmesh,effmass)

    # set dos energies
    num_samples=200
    energy_samples=np.linspace(-1,5,num_samples)
    num_kpoints_ibz=kmesh.shape[0]
    volume=np.linalg.det(rec_basis) 
    dos=np.zeros((1,num_samples),dtype='double')
    int_dos=np.zeros((1,num_samples),dtype='double')
    
    # call tetrahedron to calculate density of states
    tetrahedron_method_interface.calc_dos_lin_tetra(energies,spg_grid_address,spg_grid_mapping, kmesh, rec_basis, energy_samples, num_energy_samples, 1, num_kpoints_ibz, volume, dos, int_dos)

    # dump data to some file, fix tomorrow
    
def generate_energy(kmesh,effmass):
    # constant so that units are okeyish (eV energies)
    bandunit=3.81
    k2=kmesh*kmesh
    return bandunit*np.sum(k2/effmass,axis=1)

def dir_to_cart(kmesh,rec_basis):
    cart = np.dot(v,rec_basis)
    return cart

if __name__ == '__main__':
    main()
