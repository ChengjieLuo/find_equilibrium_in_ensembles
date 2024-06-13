import numpy as np
import matplotlib.pyplot as plt
from analyze_src.ufun import *

def Fourier(fx):
    transformed=np.fft.fftn(fx)
    return transformed

def Inverse_Fourier(gw):
    ### note dw*N is needed for inverse case
    transformed=np.fft.ifftn(gw)
    return transformed

# def cal_laplace(phi_now,k2s):
#         phik=Fourier(phi_now)
#         phi_k2_k=-phik*k2s
#         laplace_phi=np.real(Inverse_Fourier(phi_k2_k))
#         return laplace_phi

def convol_x_k(X,K):
        return np.real(Inverse_Fourier(Fourier(X)*K))

# def solve_laplace(phi_now,k2s):
#         phik=Fourier(phi_now)
#         phi_k2_k=-phik/k2s
#         laplace_phi=np.real(Inverse_Fourier(phi_k2_k))
#         return laplace_phi

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
    print(f"{Lbetas=}")
    
    
    num_comps,  num_beta, num_coord= omegas.shape[:3]

    print(num_comps,  num_beta, num_coord)

    #### k**2 which is beta dependent
    k2s=np.zeros((num_beta,num_coord))
    for i in range(num_beta):
        k2s[i]=(np.fft.fftfreq(num_coord, d=(Lbetas[i]/num_coord))*(2.0*np.pi))**2

    inverse_k2s=1./k2s
    for i in range(num_beta):
        inverse_k2s[i].flat[0]=0.0


    incomp=np.einsum('ijk->j',phis)/num_coord-1
    print(f'{incomp=}')
    incomp=np.max(np.abs(1.0-np.einsum('ijk->jk',phis)))
    print(f'{incomp=}')

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

    print(f"{total_energy=}")
    print(f"{energy_beta=}")
    sum=entropyterm+chiterm+psizsphiterm+kappaterm+nablapsiterm+zeroterm
    print(f"{sum=}")
    sum1=np.mean(energy_beta*Js)+zeroterm
    print(f"{sum1=}")
    # print(f"{np.mean(energy_beta*Js)-np.sum(phi_means/Ls*np.log(phi_means))}")

    # ### now calculate chemical reactions
    # omega_temp = np.tensordot(chis,phis,1)
    # for itr_comp in range(num_comps):
    #     for itr_beta in range(num_beta):
    #         # omega_temp[itr_comp,itr_beta]+=(-kappas[itr_comp]*convol_x_k(phis[itr_comp,itr_beta]-np.sum(phis[:,itr_beta],axis=0),-k2s[itr_beta])+psi[itr_beta]*zs[itr_comp])
    #         # +zs[itr_comp]*zetas[itr_beta]
    #         omega_temp[itr_comp,itr_beta]+=(-kappas[itr_comp]*convol_x_k(phis[itr_comp,itr_beta],-k2s[itr_beta])+psi[itr_beta]*zs[itr_comp])
    #         omega_temp[itr_comp,itr_beta]+=np.log(phis[itr_comp,itr_beta])/Ls[itr_comp]+1/Ls[itr_comp]
    #         # +zs[itr_comp]*zetas[itr_beta]
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

    return entropyterm, chiterm, psizsphiterm, kappaterm,nablapsiterm, zeroterm, total_energy





for chi in [-7.0, -6.0, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0]:
    zps=np.linspace(0.1,0.9,30)
    zms=np.linspace(-0.1,-0.9,30)
# cpp_code/data_from_dep/data_1_comp/zs_0.817241379310345_-0.210344827586207_-1.000000000000000_1.000000000000000_0.000000000000000__chi-5.0_alldata.txt
# for chi in [-7.0]:
#     zps=np.array([0.817241379310345])
#     zms=np.array([-0.210344827586207])
    # zps=np.array([0.1])
    # zms=np.array([-0.237931034482759])
    print(f"{zps=}")
    print(f"{zms=}")
    Nzp=len(zps)
    Nzm=len(zms)
    Lall=np.zeros((Nzp,Nzm))
    errorsall=np.zeros((Nzp,Nzm))
    errorsLbetaall=np.zeros((Nzp,Nzm))
    energyall=np.zeros((Nzp,Nzm))
    finished=np.zeros((Nzp,Nzm))
    abscharges=np.zeros((Nzp,Nzm))

    entropy_all=np.zeros((Nzp,Nzm))
    chi_all=np.zeros((Nzp,Nzm))
    psizsphi_all=np.zeros((Nzp,Nzm))
    kappa_all=np.zeros((Nzp,Nzm))
    nablapsi_all=np.zeros((Nzp,Nzm))
    zero_all=np.zeros((Nzp,Nzm))
    totale_all=np.zeros((Nzp,Nzm))

    for izp,zp in enumerate(zps):
        for izm,zm in enumerate(zms):
            num_comps=5
            num_beta=1
            num_coord=128

            # chi=-5
            Ls=np.array([10.,10.,1.,1.,1.])

            zs=np.array([zp ,zm, -1.0 ,1.0, 0.0])

            #### phi_means for ions is controlled by neutral charge 
            phi_means=np.array([0.1,0.1,0,0,0])
            phi_means[2]=-zs[0]*phi_means[0]/zs[2]
            phi_means[3]=-zs[1]*phi_means[1]/zs[3]
            phi_means[4]=1-phi_means[:4].sum()
            # print(f"total charge={np.sum(zs*phi_means)}")


            kappas=np.array([1.,1.,1.,1.,0.])*10
            v=100
            # [file filenameomega || random amplitude] 
            omega_type="random"
            omega_value=5.0
            steps_inner1= 100001
            steps_inner2= 300001
            acceptance_omega= 0.001
            acceptance_J= 0.001
            acceptance_Lbeta= 10
            acceptance_zeta= 10
            # [file filenameJs||random amplitude] 
            Js_type="random"
            Js_value=0.0

            # [file filenameLbetas||equal valueL] 
            Lbeta_type="equal"
            Lbeta_value=30

            # [file filenamezetas||equal valuezeta] 
            zetas_type="equal"
            zetas_value=0.0

            flag_C= "true"
            C= 100.0

            # [file filenameps||equal valuep]
            ps_type="equal"
            ps_value=1.0

            flag_zetas= "true"
            flag_ps= "false"

            flag_save_separate="true"

            ###convert vectors to strings
            phi_means_str=" "
            Ls_str=" "
            zs_str=" "
            kappas_str=" "
            for i in range(num_comps):
                phi_means_str=phi_means_str+f"{phi_means[i]:.15f} "
                Ls_str=Ls_str+f"{Ls[i]:.15f} "
                zs_str=zs_str+f"{zs[i]:.15f} "
                kappas_str=kappas_str+f"{kappas[i]:.15f} "


            savedata_pre= f"{num_comps}_{num_beta}_{num_coord}_{phi_means_str}_{chi}_{Ls_str}_{zs_str}_{kappas_str}_{v}_{omega_type}_{omega_value}_{steps_inner1}_{steps_inner2}_{acceptance_omega}_{acceptance_J}_{acceptance_Lbeta}_{acceptance_zeta}_{Js_type}_{Js_value}_{Lbeta_type}_{Lbeta_value}_{zetas_type}_{zetas_value}_{flag_C}_{C}_{ps_type}_{ps_value}_{flag_zetas}_{flag_ps}_{flag_save_separate}"

            # import re
            # savedata_pre=re.sub(' +','_',savedata_pre)

            savedata_pre= f"zs{zs_str}_chi{chi}"
            savedata_pre='_'.join(savedata_pre.split(' '))

            # print(f"{savedata_pre=}")

            savedata_folder1='./data_1_comp/'
            filename1=f'{savedata_folder1+savedata_pre}_alldata.txt'

            savedata_folder2='./data_2_comp/'
            filename2=f'{savedata_folder2+savedata_pre}_alldata.txt'

            savedata_folder2constantomega='./data_2_comp_constantomega/'
            filename2constantomega=f'{savedata_folder2constantomega+savedata_pre}_alldata.txt'

            try:
                # Lbetas,omegas,phis,psi,Js,total_energy,errors=analyze_from_filename(filename,5,1,128)
                ####errors=max_abs_incomp,max_omega_diff,max_J_diff,max_Lbeta_error,max_zeta_error
                


                chis=np.zeros((num_comps,num_comps))
                chis[0,1]=chis[1,0]=chi
                chis-=np.min(chis)
                print(chis)

                Lbetas1,omegas1,phis1,psi1,Js1,total_energy1,errors1=loadfile(filename1,5,1,128)
                psi1=np.loadtxt(f'{savedata_folder1+savedata_pre}_psi.txt').reshape((1,128))
                entropyterm1, chiterm1, psizsphiterm1, kappaterm1,nablapsiterm1, zeroterm1, totale1=check_eq(phi_means,
                    chis,
                    Ls,
                    zs,
                    kappas,
                    v,
                    omegas1,
                    phis1,
                    psi1,
                    Js1,
                    Lbetas1)
                print(f"{total_energy1=}")
                print(entropyterm1, chiterm1, psizsphiterm1, kappaterm1,nablapsiterm1, zeroterm1, totale1)
                # print(f"{psi1=}")
            


                Lbetas2,omegas2,phis2,psi2,Js2,total_energy2,errors2=loadfile(filename2,5,2,128)
                print(errors2)
                psi2=np.loadtxt(f'{savedata_folder2+savedata_pre}_psi.txt').reshape((2,128))
                entropyterm2, chiterm2, psizsphiterm2, kappaterm2,nablapsiterm2, zeroterm2, totale2=check_eq(phi_means,
                    chis,
                    Ls,
                    zs,
                    kappas,
                    v,
                    omegas2,
                    phis2,
                    psi2,
                    Js2,
                    Lbetas2)
                # print(phis2)
                # print(Js2)
                print(f"{total_energy2=}")
                print(entropyterm2, chiterm2, psizsphiterm2, kappaterm2,nablapsiterm2, zeroterm2, totale2)
                # print(f"{psi2=}")


                Lbetas2constantomega,omegas2constantomega,phis2constantomega,psi2constantomega,Js2constantomega,total_energy2constantomega,errors2constantomega=loadfile(filename2constantomega,5,2,128)
                psi2constantomega=np.loadtxt(f'{savedata_folder2constantomega+savedata_pre}_psi.txt').reshape((2,128))
                entropyterm2constantomega, chiterm2constantomega, psizsphiterm2constantomega, kappaterm2constantomega,nablapsiterm2constantomega, zeroterm2constantomega, totale2constantomega=check_eq(phi_means,
                    chis,
                    Ls,
                    zs,
                    kappas,
                    v,
                    omegas2constantomega,
                    phis2constantomega,
                    psi2constantomega,
                    Js2constantomega,
                    Lbetas2constantomega)
                
                print(f"{total_energy2constantomega=}")
            

                title=f'zs_{zs[0]:.3f}{zs[1]:.3f}{zs[2]:.3f}{zs[3]:.3f}_chi{chi:.3f}'
                print(title)
                pattern_amplitude_value=1e-3
                different_phase_value=1e-3
                if(total_energy1<=total_energy2constantomega and total_energy1<=total_energy2):

                    compart1amplitude=np.max(phis1[-1,0,:])-np.min(phis1[-1,0,:])

                    if (compart1amplitude<pattern_amplitude_value):

                        properties=['r-','b-','r--','b--','k:']
                        names=['p+','p-','e-','e+','S']
                        # plt.figure()
                        # for itr_beta in range(1):
                        #     xs=np.arange(128)*Lbetas1[itr_beta]/128
                        #     for itr_comp in range(5):
                        #         plt.plot(xs,phis1[itr_comp,itr_beta],properties[itr_comp],label=str(names[itr_comp]))
                        #     plt.legend()
                        # plt.title(f'{title}_Homo_error{errors1[1]}')
                        # plt.xlabel('x')
                        # plt.ylabel('phi')
                        # plt.savefig(f'figures/phis/phis_{savedata_pre}.pdf')
                        # plt.close()
                        Lall[izp,izm]=Lbetas1[0]
                        errorsall[izp,izm]=errors1[1]
                        errorsLbetaall[izp,izm]=errors1[3]
                        energyall[izp,izm]=total_energy1
                        finished[izp,izm]=3               
                        abscharges[izp,izm]=np.mean(np.abs(np.einsum('i,ij->j',zs,phis1[:,0,:])))
                        entropy_all[izp,izm]=entropyterm1
                        chi_all[izp,izm]=chiterm1
                        psizsphi_all[izp,izm]=psizsphiterm1
                        kappa_all[izp,izm]=kappaterm1
                        nablapsi_all[izp,izm]=nablapsiterm1
                        zero_all[izp,izm]=zeroterm1
                        totale_all[izp,izm]=totale1
                    else:
                        properties=['r-','b-','r--','b--','k:']
                        names=['p+','p-','e-','e+','S']
                        # plt.figure()
                        # for itr_beta in range(1):
                        #     xs=np.arange(128)*Lbetas1[itr_beta]/128
                        #     for itr_comp in range(5):
                        #         plt.plot(xs,phis1[itr_comp,itr_beta],properties[itr_comp],label=str(names[itr_comp]))
                        #     plt.legend()
                        # plt.title(f'{title}_Pattern_error{errors1[1]}')
                        # plt.xlabel('x')
                        # plt.ylabel('phi')
                        # plt.savefig(f'figures/phis/phis_{savedata_pre}.pdf')
                        # plt.close()
                        Lall[izp,izm]=Lbetas1[0]
                        errorsall[izp,izm]=errors1[1]
                        errorsLbetaall[izp,izm]=errors1[3]
                        energyall[izp,izm]=total_energy1
                        finished[izp,izm]=0
                        abscharges[izp,izm]=np.mean(np.abs(np.einsum('i,ij->j',zs,phis1[:,0,:])))
                        entropy_all[izp,izm]=entropyterm1
                        chi_all[izp,izm]=chiterm1
                        psizsphi_all[izp,izm]=psizsphiterm1
                        kappa_all[izp,izm]=kappaterm1
                        nablapsi_all[izp,izm]=nablapsiterm1
                        zero_all[izp,izm]=zeroterm1
                        totale_all[izp,izm]=totale1


                elif(total_energy2<=total_energy2constantomega and total_energy2<=total_energy1):

                    ##### this case can be either (two identical patterned compartments and hence in fact only one patterned phase) or (one flat and one patterned phase)

                    compart1amplitude=np.max(phis2[-1,0,:])-np.min(phis2[-1,0,:])
                    compart2amplitude=np.max(phis2[-1,1,:])-np.min(phis2[-1,1,:])

                    compart1max=np.max(phis2[-1,0,:])
                    compart2max=np.max(phis2[-1,1,:])

                    compart1min=np.min(phis2[-1,0,:])
                    compart2min=np.min(phis2[-1,1,:])


                    if(compart1amplitude>1e-2):
                        flag1=0 ##pattern
                    else:
                        flag1=1 ##flat
                    if(compart2amplitude>1e-2):
                        flag2=0 ##pattern
                    else:
                        flag2=1 ##flat

                    if(flag1!=flag2):#### real mixed case
                        properties=['r-','b-','r--','b--','k:']
                        names=['p+','p-','e-','e+','S']
                        # fig,axs=plt.subplots(2,figsize=(5,10))
                        # for itr_beta in range(2):
                        #     xs=np.arange(128)*Lbetas2[itr_beta]/128
                        #     for itr_comp in range(5):
                        #         axs[itr_beta].plot(xs,phis2[itr_comp,itr_beta],properties[itr_comp],label=str(names[itr_comp]))
                        #     axs[itr_beta].legend()
                        #     axs[itr_beta].set_xlabel('x')
                        #     axs[itr_beta].set_ylabel('phi')
                        # fig.suptitle(f'{title}_PatternFlat_error{errors2[1]}')
                        # fig.savefig(f'figures/phis/phis_{savedata_pre}.pdf')
                        # plt.close()

                        if(flag1==0 and flag2==1):
                            Lall[izp,izm]=Lbetas2[0]##length is of the patterned phase
                            abscharges[izp,izm]=np.mean(np.abs(np.einsum('i,ij->j',zs,phis2[:,0,:])))
                        if(flag1==1 and flag2==0):
                            Lall[izp,izm]=Lbetas2[1]
                            abscharges[izp,izm]=np.mean(np.abs(np.einsum('i,ij->j',zs,phis2[:,1,:])))

                        errorsall[izp,izm]=errors2[1]
                        errorsLbetaall[izp,izm]=errors2[3]
                        energyall[izp,izm]=total_energy2
                        finished[izp,izm]=1
                        
                        entropy_all[izp,izm]=entropyterm2
                        chi_all[izp,izm]=chiterm2
                        psizsphi_all[izp,izm]=psizsphiterm2
                        kappa_all[izp,izm]=kappaterm2
                        nablapsi_all[izp,izm]=nablapsiterm2
                        zero_all[izp,izm]=zeroterm2
                        totale_all[izp,izm]=totale2
                    elif(flag1==flag2 and flag1==1):###flat phases

                        if (np.abs(compart1max-compart2max)<different_phase_value):##Homogeneous
                            properties=['r-','b-','r--','b--','k:']
                            names=['p+','p-','e-','e+','S']
                            # plt.figure()
                            # for itr_beta in range(1):
                            #     xs=np.arange(128)*Lbetas2[itr_beta]/128
                            #     for itr_comp in range(5):
                            #         plt.plot(xs,phis2[itr_comp,itr_beta],properties[itr_comp],label=str(names[itr_comp]))
                            #     plt.legend()
                            #     plt.xlabel('x')
                            #     plt.ylabel('phi')
                            # plt.title(f'{title}_Homo_error{errors2[1]}')
                            # plt.savefig(f'figures/phis/phis_{savedata_pre}.pdf')
                            # plt.close()
                            abscharges[izp,izm]=np.mean(np.abs(np.einsum('i,ij->j',zs,phis2[:,0,:])))
                            Lall[izp,izm]=np.min(Lbetas2)
                            errorsall[izp,izm]=errors2[1]
                            errorsLbetaall[izp,izm]=errors2[3]
                            energyall[izp,izm]=total_energy2
                            finished[izp,izm]=3
                            entropy_all[izp,izm]=entropyterm2
                            chi_all[izp,izm]=chiterm2
                            psizsphi_all[izp,izm]=psizsphiterm2
                            kappa_all[izp,izm]=kappaterm2
                            nablapsi_all[izp,izm]=nablapsiterm2
                            zero_all[izp,izm]=zeroterm2
                            totale_all[izp,izm]=totale2

                        else: ### Flat+Flat
                            properties=['r-','b-','r--','b--','k:']
                            names=['p+','p-','e-','e+','S']
                            fig,axs=plt.subplots(2,figsize=(5,10))
                            for itr_beta in range(2):
                                xs=np.arange(128)*Lbetas2[itr_beta]/128
                                for itr_comp in range(5):
                                    axs[itr_beta].plot(xs,phis2[itr_comp,itr_beta],properties[itr_comp],label=str(names[itr_comp]))
                                axs[itr_beta].legend()
                                axs[itr_beta].set_xlabel('x')
                                axs[itr_beta].set_ylabel('phi')
                            fig.suptitle(f'{title}_FlatFlat_error{errors2[1]}')
                            fig.savefig(f'figures/phis/phis_{savedata_pre}.pdf')
                            plt.close()
                            abscharges[izp,izm]=np.mean(np.abs(np.einsum('i,ij->j',zs,phis2[:,np.argmin(Lbetas2),:])))
                            
                            Lall[izp,izm]=np.min(Lbetas2)
                            errorsall[izp,izm]=errors2[1]
                            errorsLbetaall[izp,izm]=errors2[3]
                            energyall[izp,izm]=total_energy2
                            finished[izp,izm]=2
                            entropy_all[izp,izm]=entropyterm2
                            chi_all[izp,izm]=chiterm2
                            psizsphi_all[izp,izm]=psizsphiterm2
                            kappa_all[izp,izm]=kappaterm2
                            nablapsi_all[izp,izm]=nablapsiterm2
                            zero_all[izp,izm]=zeroterm2
                            totale_all[izp,izm]=totale2
                    elif(flag1==flag2 and flag1==0):###two patterns
                        if(np.abs(compart1max-compart2max)<different_phase_value and np.abs(compart1min-compart2min)<different_phase_value and np.abs(compart1amplitude-compart2amplitude)<different_phase_value):### same patterns in the two phases, so plot only one phase
                            index_beta=np.argmin(Lbetas2)
                            properties=['r-','b-','r--','b--','k:']
                            names=['p+','p-','e-','e+','S']
                            # plt.figure()
                            # for itr_beta in range(1):
                            #     xs=np.arange(128)*Lbetas2[index_beta]/128
                            #     for itr_comp in range(5):
                            #         plt.plot(xs,phis2[itr_comp,index_beta],properties[itr_comp],label=str(names[itr_comp]))
                            #     plt.legend()
                            # plt.title(f'{title}_Pattern_error{errors2[1]}')
                            # plt.xlabel('x')
                            # plt.ylabel('phi')
                            # plt.savefig(f'figures/phis/phis_{savedata_pre}.pdf')
                            # plt.close()

                            abscharges[izp,izm]=np.mean(np.abs(np.einsum('i,ij->j',zs,phis2[:,index_beta,:])))
                            Lall[izp,izm]=Lbetas2[index_beta]
                            errorsall[izp,izm]=errors2[1]
                            errorsLbetaall[izp,izm]=errors2[3]
                            energyall[izp,izm]=total_energy2
                            finished[izp,izm]=0
                            entropy_all[izp,izm]=entropyterm2
                            chi_all[izp,izm]=chiterm2
                            psizsphi_all[izp,izm]=psizsphiterm2
                            kappa_all[izp,izm]=kappaterm2
                            nablapsi_all[izp,izm]=nablapsiterm2
                            zero_all[izp,izm]=zeroterm2
                            totale_all[izp,izm]=totale2
                        else:### pattern+pattern
                            properties=['r-','b-','r--','b--','k:']
                            names=['p+','p-','e-','e+','S']
                            fig,axs=plt.subplots(2,figsize=(5,10))
                            for itr_beta in range(2):
                                xs=np.arange(128)*Lbetas2[itr_beta]/128
                                for itr_comp in range(5):
                                    axs[itr_beta].plot(xs,phis2[itr_comp,itr_beta],properties[itr_comp],label=str(names[itr_comp]))
                                axs[itr_beta].legend()
                                axs[itr_beta].set_xlabel('x')
                                axs[itr_beta].set_ylabel('phi')
                            fig.suptitle(f'{title}_PatternPattern_error{errors2[1]}')
                            fig.savefig(f'figures/phis/phis_{savedata_pre}.pdf')
                            plt.close()

                            if(flag1==0 and flag2==1):
                                Lall[izp,izm]=Lbetas2[0]
                            if(flag1==1 and flag2==0):
                                Lall[izp,izm]=Lbetas2[1]

                            errorsall[izp,izm]=errors2[1]
                            errorsLbetaall[izp,izm]=errors2[3]
                            energyall[izp,izm]=total_energy2
                            finished[izp,izm]=4
                            entropy_all[izp,izm]=entropyterm2
                            chi_all[izp,izm]=chiterm2
                            psizsphi_all[izp,izm]=psizsphiterm2
                            kappa_all[izp,izm]=kappaterm2
                            nablapsi_all[izp,izm]=nablapsiterm2
                            zero_all[izp,izm]=zeroterm2
                            totale_all[izp,izm]=totale2

                elif(total_energy2constantomega<=total_energy2 and total_energy2constantomega<=total_energy1):

                    ### rename to reuse the code above
                    Lbeats2=Lbetas2constantomega
                    omegas2=omegas2constantomega
                    phis2=phis2constantomega
                    psi2=psi2constantomega
                    Js2=Js2constantomega
                    total_energy2=total_energy2constantomega
                    errors2=errors2constantomega

                    entropyterm2=entropyterm2constantomega
                    chiterm2=chiterm2constantomega
                    psizsphiterm2=psizsphiterm2constantomega
                    kappaterm2=kappaterm2constantomega
                    nablapsiterm2=nablapsiterm2constantomega
                    zeroterm2=zeroterm2constantomega
                    totale2=totale2constantomega

                    ##### this case can be either (two identical patterned compartments and hence in fact only one patterned phase) or (one flat and one patterned phase)

                    compart1amplitude=np.max(phis2[-1,0,:])-np.min(phis2[-1,0,:])
                    compart2amplitude=np.max(phis2[-1,1,:])-np.min(phis2[-1,1,:])

                    compart1max=np.max(phis2[-1,0,:])
                    compart2max=np.max(phis2[-1,1,:])

                    compart1min=np.min(phis2[-1,0,:])
                    compart2min=np.min(phis2[-1,1,:])


                    if(compart1amplitude>1e-2):
                        flag1=0 ##pattern
                    else:
                        flag1=1 ##flat
                    if(compart2amplitude>1e-2):
                        flag2=0 ##pattern
                    else:
                        flag2=1 ##flat

                    if(flag1!=flag2):#### real mixed case
                        properties=['r-','b-','r--','b--','k:']
                        names=['p+','p-','e-','e+','S']
                        # fig,axs=plt.subplots(2,figsize=(5,10))
                        # for itr_beta in range(2):
                        #     xs=np.arange(128)*Lbetas2[itr_beta]/128
                        #     for itr_comp in range(5):
                        #         axs[itr_beta].plot(xs,phis2[itr_comp,itr_beta],properties[itr_comp],label=str(names[itr_comp]))
                        #     axs[itr_beta].legend()
                        #     axs[itr_beta].set_xlabel('x')
                        #     axs[itr_beta].set_ylabel('phi')
                        # fig.suptitle(f'{title}_PatternFlat_error{errors2[1]}')
                        # fig.savefig(f'figures/phis/phis_{savedata_pre}.pdf')
                        # plt.close()

                        if(flag1==0 and flag2==1):
                            Lall[izp,izm]=Lbetas2[0]##length is of the patterned phase
                            abscharges[izp,izm]=np.mean(np.abs(np.einsum('i,ij->j',zs,phis2[:,0,:])))
                        if(flag1==1 and flag2==0):
                            Lall[izp,izm]=Lbetas2[1]
                            abscharges[izp,izm]=np.mean(np.abs(np.einsum('i,ij->j',zs,phis2[:,1,:])))

                        errorsall[izp,izm]=errors2[1]
                        errorsLbetaall[izp,izm]=errors2[3]
                        energyall[izp,izm]=total_energy2
                        finished[izp,izm]=1
                        entropy_all[izp,izm]=entropyterm2
                        chi_all[izp,izm]=chiterm2
                        psizsphi_all[izp,izm]=psizsphiterm2
                        kappa_all[izp,izm]=kappaterm2
                        nablapsi_all[izp,izm]=nablapsiterm2
                        zero_all[izp,izm]=zeroterm2
                        totale_all[izp,izm]=totale2
                    elif(flag1==flag2 and flag1==1):###flat phases

                        if (np.abs(compart1max-compart2max)<different_phase_value):##Homogeneous
                            properties=['r-','b-','r--','b--','k:']
                            names=['p+','p-','e-','e+','S']
                            # plt.figure()
                            # for itr_beta in range(1):
                            #     xs=np.arange(128)*Lbetas2[itr_beta]/128
                            #     for itr_comp in range(5):
                            #         plt.plot(xs,phis2[itr_comp,itr_beta],properties[itr_comp],label=str(names[itr_comp]))
                            #     plt.legend()
                            #     plt.xlabel('x')
                            #     plt.ylabel('phi')
                            # plt.title(f'{title}_Homo_error{errors2[1]}')
                            # plt.savefig(f'figures/phis/phis_{savedata_pre}.pdf')
                            # plt.close()
                            abscharges[izp,izm]=np.mean(np.abs(np.einsum('i,ij->j',zs,phis2[:,0,:])))
                            Lall[izp,izm]=np.min(Lbetas2)
                            errorsall[izp,izm]=errors2[1]
                            errorsLbetaall[izp,izm]=errors2[3]
                            energyall[izp,izm]=total_energy2
                            finished[izp,izm]=3
                            entropy_all[izp,izm]=entropyterm2
                            chi_all[izp,izm]=chiterm2
                            psizsphi_all[izp,izm]=psizsphiterm2
                            kappa_all[izp,izm]=kappaterm2
                            nablapsi_all[izp,izm]=nablapsiterm2
                            zero_all[izp,izm]=zeroterm2
                            totale_all[izp,izm]=totale2

                        else: ### Flat+Flat
                            properties=['r-','b-','r--','b--','k:']
                            names=['p+','p-','e-','e+','S']
                            fig,axs=plt.subplots(2,figsize=(5,10))
                            for itr_beta in range(2):
                                xs=np.arange(128)*Lbetas2[itr_beta]/128
                                for itr_comp in range(5):
                                    axs[itr_beta].plot(xs,phis2[itr_comp,itr_beta],properties[itr_comp],label=str(names[itr_comp]))
                                axs[itr_beta].legend()
                                axs[itr_beta].set_xlabel('x')
                                axs[itr_beta].set_ylabel('phi')
                            fig.suptitle(f'{title}_FlatFlat_error{errors2[1]}')
                            fig.savefig(f'figures/phis/phis_{savedata_pre}.pdf')
                            plt.close()
                            abscharges[izp,izm]=np.mean(np.abs(np.einsum('i,ij->j',zs,phis2[:,np.argmin(Lbetas2),:])))
                            
                            Lall[izp,izm]=np.min(Lbetas2)
                            errorsall[izp,izm]=errors2[1]
                            errorsLbetaall[izp,izm]=errors2[3]
                            energyall[izp,izm]=total_energy2
                            finished[izp,izm]=2
                            entropy_all[izp,izm]=entropyterm2
                            chi_all[izp,izm]=chiterm2
                            psizsphi_all[izp,izm]=psizsphiterm2
                            kappa_all[izp,izm]=kappaterm2
                            nablapsi_all[izp,izm]=nablapsiterm2
                            zero_all[izp,izm]=zeroterm2
                            totale_all[izp,izm]=totale2
                    elif(flag1==flag2 and flag1==0):###two patterns
                        if(np.abs(compart1max-compart2max)<different_phase_value and np.abs(compart1min-compart2min)<different_phase_value and np.abs(compart1amplitude-compart2amplitude)<different_phase_value):### same patterns in the two phases, so plot only one phase
                            index_beta=np.argmin(Lbetas2)
                            properties=['r-','b-','r--','b--','k:']
                            names=['p+','p-','e-','e+','S']
                            # plt.figure()
                            # for itr_beta in range(1):
                            #     xs=np.arange(128)*Lbetas2[index_beta]/128
                            #     for itr_comp in range(5):
                            #         plt.plot(xs,phis2[itr_comp,index_beta],properties[itr_comp],label=str(names[itr_comp]))
                            #     plt.legend()
                            # plt.title(f'{title}_Pattern_error{errors2[1]}')
                            # plt.xlabel('x')
                            # plt.ylabel('phi')
                            # plt.savefig(f'figures/phis/phis_{savedata_pre}.pdf')
                            # plt.close()

                            abscharges[izp,izm]=np.mean(np.abs(np.einsum('i,ij->j',zs,phis2[:,index_beta,:])))
                            Lall[izp,izm]=Lbetas2[index_beta]
                            errorsall[izp,izm]=errors2[1]
                            errorsLbetaall[izp,izm]=errors2[3]
                            energyall[izp,izm]=total_energy2
                            finished[izp,izm]=0
                            entropy_all[izp,izm]=entropyterm2
                            chi_all[izp,izm]=chiterm2
                            psizsphi_all[izp,izm]=psizsphiterm2
                            kappa_all[izp,izm]=kappaterm2
                            nablapsi_all[izp,izm]=nablapsiterm2
                            zero_all[izp,izm]=zeroterm2
                            totale_all[izp,izm]=totale2
                        else:### pattern+pattern
                            properties=['r-','b-','r--','b--','k:']
                            names=['p+','p-','e-','e+','S']
                            fig,axs=plt.subplots(2,figsize=(5,10))
                            for itr_beta in range(2):
                                xs=np.arange(128)*Lbetas2[itr_beta]/128
                                for itr_comp in range(5):
                                    axs[itr_beta].plot(xs,phis2[itr_comp,itr_beta],properties[itr_comp],label=str(names[itr_comp]))
                                axs[itr_beta].legend()
                                axs[itr_beta].set_xlabel('x')
                                axs[itr_beta].set_ylabel('phi')
                            fig.suptitle(f'{title}_PatternPattern_error{errors2[1]}')
                            fig.savefig(f'figures/phis/phis_{savedata_pre}.pdf')
                            plt.close()

                            if(flag1==0 and flag2==1):
                                Lall[izp,izm]=Lbetas2[0]
                            if(flag1==1 and flag2==0):
                                Lall[izp,izm]=Lbetas2[1]

                            errorsall[izp,izm]=errors2[1]
                            errorsLbetaall[izp,izm]=errors2[3]
                            energyall[izp,izm]=total_energy2
                            finished[izp,izm]=4
                            entropy_all[izp,izm]=entropyterm2
                            chi_all[izp,izm]=chiterm2
                            psizsphi_all[izp,izm]=psizsphiterm2
                            kappa_all[izp,izm]=kappaterm2
                            nablapsi_all[izp,izm]=nablapsiterm2
                            zero_all[izp,izm]=zeroterm2
                            totale_all[izp,izm]=totale2

            except:
                print(f"error in {savedata_pre}")
    zpplot,zmplot=np.meshgrid(zps,-zms,indexing='ij')
    # plt.figure()
    # plt.pcolor(zpplot,zmplot,Lall)
    # plt.xlabel('zp')
    # plt.ylabel('-zm')
    # plt.colorbar()
    # plt.title('Length')
    # plt.savefig(f'figures/lowestenergy_Lall_chi{chi}.pdf',dpi=300)
    # plt.close()

    # plt.figure()
    # plt.pcolor(zpplot,zmplot,np.log10(errorsall))
    # plt.xlabel('zp')
    # plt.ylabel('-zm')
    # plt.colorbar()
    # plt.title('log10(error omega)')
    # plt.savefig(f'figures/lowestenergy_erroromega_chi{chi}.pdf',dpi=300)
    # plt.close()

    # plt.figure()
    # plt.pcolor(zpplot,zmplot,np.log10(errorsLbetaall))
    # plt.xlabel('zp')
    # plt.ylabel('-zm')
    # plt.colorbar()
    # plt.title('log10(error Lbeta)')
    # plt.savefig(f'figures/lowestenergy_errorsLbeta_chi{chi}.pdf',dpi=300)
    # plt.close()

    # plt.figure()
    # plt.pcolor(zpplot,zmplot,energyall)
    # plt.xlabel('zp')
    # plt.ylabel('-zm')
    # plt.title('energyall')
    # plt.colorbar()
    # plt.savefig(f'figures/lowestenergy_energy_chi{chi}.pdf',dpi=300)
    # plt.close()

    # plt.figure()
    # plt.pcolor(zpplot,zmplot,finished)
    # plt.xlabel('zp')
    # plt.ylabel('-zm')
    # plt.title('types')
    # plt.colorbar()
    # plt.savefig(f'figures/lowestenergy_types_chi{chi}.pdf',dpi=300)
    # plt.close()

    # plt.figure()
    # plt.pcolor(zpplot,zmplot,abscharges)
    # plt.xlabel('zp')
    # plt.ylabel('-zm')
    # plt.title('charges')
    # plt.colorbar()
    # plt.savefig(f'figures/lowestenergy_abscharges_chi{chi}.pdf',dpi=300)
    # plt.close()

    plt.figure()
    plt.pcolor(zpplot,zmplot,entropy_all)
    plt.xlabel('zp')
    plt.ylabel('-zm')
    plt.title('entropy')
    plt.colorbar()
    plt.savefig(f'figures/lowestenergy_entropy_chi{chi}.pdf',dpi=300)
    plt.close()

    plt.figure()
    plt.pcolor(zpplot,zmplot,chi_all)
    plt.xlabel('zp')
    plt.ylabel('-zm')
    plt.title('chiterm')
    plt.colorbar()
    plt.savefig(f'figures/lowestenergy_chiterm_chi{chi}.pdf',dpi=300)
    plt.close()

    plt.figure()
    plt.pcolor(zpplot,zmplot,psizsphi_all)
    plt.xlabel('zp')
    plt.ylabel('-zm')
    plt.title('psizsphi')
    plt.colorbar()
    plt.savefig(f'figures/lowestenergy_psizsphiterm_chi{chi}.pdf',dpi=300)
    plt.close()

    plt.figure()
    plt.pcolor(zpplot,zmplot,kappa_all)
    plt.xlabel('zp')
    plt.ylabel('-zm')
    plt.title('kappa')
    plt.colorbar()
    plt.savefig(f'figures/lowestenergy_kappaterm_chi{chi}.pdf',dpi=300)
    plt.close()

    plt.figure()
    plt.pcolor(zpplot,zmplot,nablapsi_all)
    plt.xlabel('zp')
    plt.ylabel('-zm')
    plt.title('nablapsi')
    plt.colorbar()
    plt.savefig(f'figures/lowestenergy_nablapsiterm_chi{chi}.pdf',dpi=300)
    plt.close()

    
    plt.figure()
    plt.pcolor(zpplot,zmplot,zero_all)
    plt.xlabel('zp')
    plt.ylabel('-zm')
    plt.title('zero')
    plt.colorbar()
    plt.savefig(f'figures/lowestenergy_zeroterm_chi{chi}.pdf',dpi=300)
    plt.close()

    plt.figure()
    plt.pcolor(zpplot,zmplot,totale_all)
    plt.xlabel('zp')
    plt.ylabel('-zm')
    plt.title('totale')
    plt.colorbar()
    plt.savefig(f'figures/lowestenergy_totaleterm_chi{chi}.pdf',dpi=300)
    plt.close()

    # np.savetxt(f"Lall_chi{chi}.txt",Lall)
    # np.savetxt(f"energyall_chi{chi}.txt",energyall)
    # np.savetxt(f"types_chi{chi}.txt",finished)
    # np.savetxt(f"abscharges_chi{chi}.txt",abscharges)
    np.savetxt(f"entropy_chi{chi}.txt",entropy_all)
    np.savetxt(f"chi_chi{chi}.txt",chi_all)
    np.savetxt(f"psizsphi_chi{chi}.txt",psizsphi_all)
    np.savetxt(f"kappa_chi{chi}.txt",kappa_all)
    np.savetxt(f"nablapsi_chi{chi}.txt",nablapsi_all)
    np.savetxt(f"zero_chi{chi}.txt",zero_all)
    np.savetxt(f"totale_chi{chi}.txt",totale_all)

