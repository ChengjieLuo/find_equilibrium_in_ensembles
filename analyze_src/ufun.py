import numpy as np
import matplotlib.pyplot as plt
np.seterr(divide='ignore')

def loadfile(filename,num_comp,num_beta,num_coord):
    if (2*num_beta+1+5+1<num_coord):
        data=np.loadtxt(filename,usecols=np.arange(2*num_beta+1+5+1))[0]
    else:
        with open(filename) as f:
            # lines = (line for line in f if not line.startswith('#'))
            # FH = np.loadtxt(lines, delimiter=',', skiprows=1)
            for iline, line in enumerate(f):
                if iline==0:
                    data=np.array(line.split(),dtype=float)
                    break
            # print(f"data={data}")


    Lbetas=data[:num_beta]
    Js=data[num_beta:2*num_beta]
    total_energy=data[2*num_beta]
    errors=data[2*num_beta+1:2*num_beta+1+5]
    step=data[2*num_beta+1+5]
            

    data=np.loadtxt(filename,skiprows=1)
    # print(data.shape)
    phis=data[:num_comp*num_beta].reshape((num_comp,num_beta,num_coord))
    omegas=data[num_comp*num_beta:num_comp*num_beta*2].reshape((num_comp,num_beta,num_coord))
    psi=data[-num_beta:].reshape((num_beta,num_coord))


    return Lbetas,omegas,phis,psi,Js,total_energy,errors,step

def analyze_from_filename(filename,num_comp,num_beta,num_coord):
    Lbetas,omegas,phis,psi,Js,total_energy,errors,step=loadfile(filename,num_comp,num_beta,num_coord)

    properties=['r-','b-','r--','b--','k:']
    properties2=['ro','bo','r*','b*','ks']
    names=['p+','p-','e-','e+','S']

    for itr_beta in range(num_beta):
        xs=np.arange(num_coord)*Lbetas[itr_beta]/num_coord
        plt.figure()
        for itr_comp in range(num_comp):
            if(len(xs)>1):
                plt.plot(xs,phis[itr_comp,itr_beta],properties[itr_comp],label=str(names[itr_comp]))
            else:
                plt.plot(xs,phis[itr_comp,itr_beta],properties2[itr_comp],label=str(names[itr_comp]))
        plt.legend()
        plt.title(rf'$\beta={itr_beta}$')
    plt.show()
    return Lbetas,omegas,phis,psi,Js,total_energy,errors,step




### Fourier transform
def Fourier(fx):
    transformed=np.fft.fftn(fx)
    return transformed

### Inverse Fourier transform
def Inverse_Fourier(gw):
    ### note dw*N is needed for inverse case
    transformed=np.fft.ifftn(gw)
    return transformed

# def cal_laplace(phi_now,k2s):
#         phik=Fourier(phi_now)
#         phi_k2_k=-phik*k2s
#         laplace_phi=np.real(Inverse_Fourier(phi_k2_k))
#         return laplace_phi

### convolution with a kernal
def convol_x_k(X,K):
        return np.real(Inverse_Fourier(Fourier(X)*K))
# def solve_laplace(phi_now,k2s):
#         phik=Fourier(phi_now)
#         phi_k2_k=-phik/k2s
#         laplace_phi=np.real(Inverse_Fourier(phi_k2_k))
#         return laplace_phi

### check equilibrium
def check_eq(phi_means,
    chis,
    Ls,
    zs,
    kappas,
    v,
    omegas,
    phis,
    psi,
    Js,
    Lbetas):
    ### note omegas are not needed
    # print(f"{Lbetas=}")
    
    
    num_comps,  num_beta, num_coord= omegas.shape[:3]

    # print(num_comps,  num_beta, num_coord)

    #### k**2 which is beta dependent
    k2s=np.zeros((num_beta,num_coord))
    for i in range(num_beta):
        k2s[i]=(np.fft.fftfreq(num_coord, d=(Lbetas[i]/num_coord))*(2.0*np.pi))**2

    inverse_k2s=1./k2s
    for i in range(num_beta):
        inverse_k2s[i].flat[0]=0.0


    incomp=np.einsum('ijk->j',phis)/num_coord-1
    print(f'{incomp=}')
    # incomp=np.max(np.abs(1.0-np.einsum('ijk->jk',phis)))
    # print(f'{incomp=}')

    phimeancheck=np.mean(np.einsum('ijk->ij',phis)/128*Js,axis=-1)
    phimeandiff=phi_means-phimeancheck
    print(f'{phimeandiff=}')

    charge=np.zeros((num_beta))
    for itr_comp in range(num_comps):
        charge+=np.mean(phis[itr_comp]*zs[itr_comp],axis=-1)
    print(f"{charge=}")
    ### calculate total energy
    total_energy=0
    energy_beta=np.zeros((num_beta))
    
    
    chi_times_phis=np.tensordot(chis,phis,1)

    entropyterm=0
    chiterm=0
    psizsphiterm=0
    for itr_comp in range(num_comps):
        ### mean over x
        total_energy+=np.mean(np.mean(phis[itr_comp]*np.log(phis[itr_comp])/Ls[itr_comp],axis=-1)*Js)

        entropyterm+=np.mean(np.mean(phis[itr_comp]*np.log(phis[itr_comp])/Ls[itr_comp],axis=-1)*Js)

        total_energy+=np.mean(np.mean(0.5*phis[itr_comp]*chi_times_phis[itr_comp],axis=-1)*Js)

        chiterm+=np.mean(np.mean(0.5*phis[itr_comp]*chi_times_phis[itr_comp],axis=-1)*Js)

        total_energy+=np.mean(np.mean(psi*zs[itr_comp]*phis[itr_comp],axis=-1)*Js)
        psizsphiterm+=np.mean(np.mean(psi*zs[itr_comp]*phis[itr_comp],axis=-1)*Js)


        energy_beta+=np.mean(phis[itr_comp]*np.log(phis[itr_comp])/Ls[itr_comp]+ 0.5*phis[itr_comp]*chi_times_phis[itr_comp]+psi*zs[itr_comp]*phis[itr_comp],axis=-1)

    tmp=0
    for itr_beta in range(num_beta):
        for itr_comp in range(num_comps):
            tmp+=np.mean(-0.5*kappas[itr_comp]*phis[itr_comp,itr_beta]*convol_x_k(phis[itr_comp,itr_beta],-k2s[itr_beta]))*Js[itr_beta]

            energy_beta[itr_beta]+=np.mean(-0.5*kappas[itr_comp]*phis[itr_comp,itr_beta]*convol_x_k(phis[itr_comp,itr_beta],-k2s[itr_beta]))

    tmp/=num_beta
    kappaterm=tmp
    total_energy+=tmp

    tmp=0
    for itr_beta in range(num_beta):
        tmp+=(v/8./np.pi)* np.mean(psi[itr_beta]*convol_x_k(psi[itr_beta],-k2s[itr_beta]))*Js[itr_beta]

        energy_beta[itr_beta]+=(v/8./np.pi)* np.mean(psi[itr_beta]*convol_x_k(psi[itr_beta],-k2s[itr_beta]))

    tmp/=num_beta
    nablapsiterm=tmp
    total_energy+=tmp
    
    zeroterm=0
    for itr_comp in range(num_comps):
        total_energy-=phi_means[itr_comp]/Ls[itr_comp]*np.log(phi_means[itr_comp])
        zeroterm-=phi_means[itr_comp]/Ls[itr_comp]*np.log(phi_means[itr_comp])

    # print(f"{total_energy=}")
    # print(f"{energy_beta=}")
    # sum=entropyterm+chiterm+psizsphiterm+kappaterm+nablapsiterm+zeroterm
    # print(f"{sum=}")
    # sum1=np.mean(energy_beta*Js)+zeroterm
    # print(f"{sum1=}")
    # print(f"{np.mean(energy_beta*Js)-np.sum(phi_means/Ls*np.log(phi_means))}")

    # ### now calculate chemical potentials
    omega_temp = np.tensordot(chis,phis,1)
    for itr_comp in range(num_comps):
        for itr_beta in range(num_beta):
            # omega_temp[itr_comp,itr_beta]+=(-kappas[itr_comp]*convol_x_k(phis[itr_comp,itr_beta]-np.sum(phis[:,itr_beta],axis=0),-k2s[itr_beta])+psi[itr_beta]*zs[itr_comp])
            # +zs[itr_comp]*zetas[itr_beta]
            omega_temp[itr_comp,itr_beta]+=(-kappas[itr_comp]*convol_x_k(phis[itr_comp,itr_beta],-k2s[itr_beta])+psi[itr_beta]*zs[itr_comp])
            omega_temp[itr_comp,itr_beta]+=np.log(phis[itr_comp,itr_beta])/Ls[itr_comp]+1/Ls[itr_comp]
            # +zs[itr_comp]*zetas[itr_beta]
    # colors=['r','b','y','g','k']
    # markers=['o','-','v','s','--']
    # for itr_beta in range(num_beta):
    #     # plt.figure(itr_beta+1)
    #     for itr_comp in range(num_comps):
    #         plt.plot((omega_temp[itr_comp,itr_beta]-omega_temp[-1,itr_beta]),colors[itr_comp]+markers[itr_beta],mfc='w')
    # plt.ylabel(r'$\mu_i(\beta)-\mu_S(\beta)$')
    # plt.show()
    
    # for itr_beta in range(num_beta-1):
    #     plt.figure(itr_beta+1)
    #     for itr_comp in range(num_comps):
    #         plt.plot((omega_temp[itr_comp,itr_beta]-omega_temp[-1,itr_beta])-(omega_temp[itr_comp,-1]-omega_temp[-1,-1]),colors[itr_comp]+markers[itr_beta])
    #     plt.ylabel(r'$[\mu_i(\beta)-\mu_S(\beta)]-[\mu_i(\beta=-1)-\mu_S(\beta=-1)]$')
    #     plt.title(rf'$\beta={itr_beta+1}$')
    # plt.show()
    # for itr_beta in range(num_beta):
    #     plt.figure(itr_beta+1)
    #     for itr_comp in range(num_comps):
    #         plt.plot(omegas[itr_comp,itr_beta],'--')#*Js[itr_beta])
    # plt.show()

    return entropyterm, chiterm, psizsphiterm, kappaterm,nablapsiterm, zeroterm, total_energy, omega_temp


def compare_properties(p1,p2,threshold):
    for i in range(len(p1)):
        p1i=p1[i]
        p2i=p2[i]
        if np.max(np.abs(p1i-p2i))>threshold:
            return 0
    return 1

def get_phases(phi_means,
    chis,
    Ls,
    zs,
    kappas,
    v,
    omegas,
    phis,
    psi,
    Js,
    Lbetas):

    threshold=1e-4

    num_comps,  num_beta, num_coord= omegas.shape[:3]

    if(num_beta==0):
        newphis=phis.copy()
        newpsi=psi.copy()
        return newphis,newpsi,1
    else:
        properties=[]
        for ibeta in range(num_beta):
            maxphis=np.max(phis[:,ibeta,:],axis=-1)
            minphis=np.min(phis[:,ibeta,:],axis=-1)
            L=Lbetas[ibeta]
            properties.append([maxphis,minphis,L])

        unique_list=[0]
        newJs=[Js[0]]
        for i in range(1,num_beta):
            flag_new=1
            for ilist in unique_list:
                if compare_properties(properties[i],properties[ilist],threshold)==1: ## same phase
                    newJs[ilist]+=Js[i]
                    flag_new=0
                    break
            if(flag_new==1):## new phase
                unique_list.append(i)
                newJs.append(Js[i])
        print(unique_list)
        newphis=phis[:,unique_list,:]
        newpsi=psi[unique_list,:]
        return newphis,newpsi,newJs,unique_list,len(unique_list)

def plot_phis(phis,Lbetas):
    num_comps,  num_beta, num_coord= phis.shape[:3]
    properties=['r-','b-','r--','b--','k:']
    names=['p+','p-','e-','e+','S']
    for itr_beta in range(num_beta):
        xs=np.arange(num_coord)*Lbetas[itr_beta]/num_coord
        plt.figure()
        for itr_comp in range(num_comps):
            plt.plot(xs,phis[itr_comp,itr_beta],properties[itr_comp],label=str(names[itr_comp]))
        plt.legend()
        plt.title(rf'$\beta={itr_beta}$')
    plt.show()           
            


