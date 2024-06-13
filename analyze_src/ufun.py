import numpy as np
import matplotlib.pyplot as plt

def loadfile(filename,num_comp,num_beta,num_coord):
    data=np.loadtxt(filename,usecols=np.arange(2*num_beta+1+5))[0]
    Lbetas=data[:num_beta]
    Js=data[num_beta:2*num_beta]
    total_energy=data[2*num_beta]
    errors=data[2*num_beta+1:2*num_beta+1+5]
    data=np.loadtxt(filename,skiprows=1)
    # print(data.shape)
    phis=data[:num_comp*num_beta].reshape((num_comp,num_beta,num_coord))
    omegas=data[num_comp*num_beta:num_comp*num_beta*2].reshape((num_comp,num_beta,num_coord))
    psi=data[-num_beta:].reshape((num_beta,num_coord))


    return Lbetas,omegas,phis,psi,Js,total_energy,errors

def analyze_from_filename(filename,num_comp,num_beta,num_coord):
    Lbetas,omegas,phis,psi,Js,total_energy,errors=loadfile(filename,num_comp,num_beta,num_coord)

    properties=['r-','b-','r--','b--','k:']
    names=['p+','p-','e-','e+','S']

    for itr_beta in range(num_beta):
        xs=np.arange(num_coord)*Lbetas[itr_beta]/num_coord
        plt.figure()
        for itr_comp in range(num_comp):
            plt.plot(xs,phis[itr_comp,itr_beta],properties[itr_comp],label=str(names[itr_comp]))
        plt.legend()
        plt.title(rf'$\beta={itr_beta}$')
    plt.show()
    return Lbetas,omegas,phis,psi,Js,total_energy,errors

