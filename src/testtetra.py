#!/usr/bin/python
import sys
import logging
import numpy as np
import utils
import tetrahedron_method_interface
import spglib_interface
import math

def main():
    # generate som simple bands based on the std. parabolic approximation

    # set testcell
    lattice=np.array([[12.0,0.0,0.0],[0.0,12.0,0.0],[0.0,0.0,12.0]],dtype='double')
    positions=np.array([[0.0,0.0,0.0]],dtype='double')
    anumbers=np.array([1],dtype='intc')
    shift=np.array([0,0,0],dtype='intc')
    
    # rec_basis (only cubic)
    rec_basis=np.ascontiguousarray((2*math.pi*np.linalg.inv(lattice)).transpose())
                  
    # first, set kmesh
    ksampling=np.array([31,31,31],dtype='intc')
    k=np.zeros((np.prod(ksampling),3),dtype='intc')
    mapping_bz_to_ibz=np.zeros(np.prod(ksampling),dtype='intc')
    scaling=np.floor(ksampling/2.0)
    spglib_interface.get_reciprocal_mesh(ksampling,np.ascontiguousarray(lattice.T),positions,anumbers,shift,k,mapping_bz_to_ibz,is_time_reversal=True,symprec=1e-5)
    spg_k=k
    k=k/scaling/2
    k_ired=k[np.unique(mapping_bz_to_ibz,return_index=True)[1]]
    kmesh=k[np.lexsort((k[:,2],k[:,1],k[:,0]))]
    kmesh_ired=k_ired[np.lexsort((k_ired[:,2],k_ired[:,1],k_ired[:,0]))] 
    
    # transform to cart, no need to save
    kmesh_ired=dir_to_cart(kmesh_ired,rec_basis)

    np.savetxt('kmesh',np.column_stack((kmesh_ired[:,0].T,kmesh_ired[:,1].T, kmesh_ired[:,2].T)))
    
    # effective mass
    effmass=np.array([0.5,0.5,0.5])
    
    # fetch energy (single band, add artifical index)
    energies=np.ascontiguousarray(generate_energy(kmesh_ired,effmass)[None].T)

    # set dos energies
    num_samples=200
    energy_samples=np.linspace(-0.2,1.5,num_samples)
    num_kpoints_ibz=kmesh_ired.shape[0]
    volume=np.linalg.det(rec_basis) 
    dos=np.zeros((1,num_samples),dtype='double')
    int_dos=np.zeros((1,num_samples),dtype='double')

    # turn bloechl on (1) or off (0)
    bloechl=0
    
    # call tetrahedron to calculate density of states
    tetrahedron_method_interface.calc_density_of_states_interface(energies,spg_k,mapping_bz_to_ibz, ksampling, rec_basis, energy_samples, num_samples, 1, num_kpoints_ibz, volume, bloechl, dos, int_dos)

    # scale tetra dos and number of states (also multiply by two
    # since we do not have spin degeneracy in the tetrahedron per se)
    # also, to yield same units, multiply int_dos by 1e3
    dos=2.0*volume*dos/np.prod(ksampling)/np.power(2*math.pi,3.0)
    int_dos=1e3*2.0*volume*int_dos/np.prod(ksampling)/np.power(2*math.pi,3.0)

    # generate data for the parabolc band (exact analytic solution)
    # right now the effective mass terms are not done correctly,
    # be sure all the values are the same
    dos_analytic=calc_analytic_dos(energy_samples,effmass)

    # generate number of states
    number_of_states_analytic=calc_number_of_states(energy_samples,effmass)
    
    # save dos to file, only one band again
    diff_dos=(dos_analytic-dos[0])/dos_analytic
    np.savetxt('dos',np.column_stack((energy_samples.T,dos_analytic.T,dos[0].T,diff_dos.T)))

    # save number of states
    diff_ns=(number_of_states_analytic-int_dos[0])/number_of_states_analytic
    np.savetxt('ns',np.column_stack((energy_samples.T,number_of_states_analytic.T,int_dos[0].T,diff_ns.T)))

def calc_number_of_states(energy_samples,effmass):
    analytic_unit=1e9*4.0*np.power(0.5109989461/np.power(197.3269788,2.0),1.5)/np.power(math.pi,2.0)/3.0/np.sqrt(2.0)
    return analytic_unit*np.power(abs(energy_samples)*np.sum(effmass)/3,1.5)
    
def calc_analytic_dos(energy_samples,effmass):
    analytic_unit=1e6*np.sqrt(2.0)*np.power(0.5109989461/np.power(197.3269788,2.0),1.5)/np.power(math.pi,2.0)
    return analytic_unit*np.power(np.sum(effmass)/3,1.5)*np.sqrt(abs(energy_samples))

def generate_energy(kmesh,effmass):
    # constant so that units are okeyish (eV energies)
    bandunit=3.81
    k2=kmesh*kmesh
    return bandunit*np.sum(k2/effmass,axis=1)

def dir_to_cart(v,rec_basis):
    cart = np.dot(v,rec_basis)
    return cart

if __name__ == '__main__':
    main()
