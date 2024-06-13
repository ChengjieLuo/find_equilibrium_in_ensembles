import numpy as np
import matplotlib.pyplot as plt
from analyze_src.ufun import *
for chi in [-7.0, -6.0, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0]:
    zps=np.linspace(0.1,0.9,30)
    zms=np.linspace(-0.1,-0.9,30)

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
                Lbetas1,omegas1,phis1,psi1,Js1,total_energy1,errors1=loadfile(filename1,5,1,128)
                

                Lbetas2,omegas2,phis2,psi2,Js2,total_energy2,errors2=loadfile(filename2,5,2,128)

                Lbetas2constantomega,omegas2constantomega,phis2constantomega,psi2constantomega,Js2constantomega,total_energy2constantomega,errors2constantomega=loadfile(filename2constantomega,5,2,128)

                title=f'zs_{zs[0]:.3f}{zs[1]:.3f}{zs[2]:.3f}{zs[3]:.3f}_chi{chi:.3f}'
                print(title)
                pattern_amplitude_value=1e-3
                different_phase_value=1e-3
                if(total_energy1<=total_energy2constantomega and total_energy1<=total_energy2):

                    compart1amplitude=np.max(phis1[-1,0,:])-np.min(phis1[-1,0,:])

                    if (compart1amplitude<pattern_amplitude_value):

                        properties=['r-','b-','r--','b--','k:']
                        names=['p+','p-','e-','e+','S']
                        plt.figure()
                        for itr_beta in range(1):
                            xs=np.arange(128)*Lbetas1[itr_beta]/128
                            for itr_comp in range(5):
                                plt.plot(xs,phis1[itr_comp,itr_beta],properties[itr_comp],label=str(names[itr_comp]))
                            plt.legend()
                        plt.title(f'{title}_Homo_error{errors1[1]}')
                        plt.xlabel('x')
                        plt.ylabel('phi')
                        plt.savefig(f'figures/phis/phis_{savedata_pre}.pdf')
                        plt.close()
                        Lall[izp,izm]=Lbetas1[0]
                        errorsall[izp,izm]=errors1[1]
                        errorsLbetaall[izp,izm]=errors1[3]
                        energyall[izp,izm]=total_energy1
                        finished[izp,izm]=3
                    else:
                        properties=['r-','b-','r--','b--','k:']
                        names=['p+','p-','e-','e+','S']
                        plt.figure()
                        for itr_beta in range(1):
                            xs=np.arange(128)*Lbetas1[itr_beta]/128
                            for itr_comp in range(5):
                                plt.plot(xs,phis1[itr_comp,itr_beta],properties[itr_comp],label=str(names[itr_comp]))
                            plt.legend()
                        plt.title(f'{title}_Pattern_error{errors1[1]}')
                        plt.xlabel('x')
                        plt.ylabel('phi')
                        plt.savefig(f'figures/phis/phis_{savedata_pre}.pdf')
                        plt.close()
                        Lall[izp,izm]=Lbetas1[0]
                        errorsall[izp,izm]=errors1[1]
                        errorsLbetaall[izp,izm]=errors1[3]
                        energyall[izp,izm]=total_energy1
                        finished[izp,izm]=0


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
                        fig,axs=plt.subplots(2,figsize=(5,10))
                        for itr_beta in range(2):
                            xs=np.arange(128)*Lbetas2[itr_beta]/128
                            for itr_comp in range(5):
                                axs[itr_beta].plot(xs,phis2[itr_comp,itr_beta],properties[itr_comp],label=str(names[itr_comp]))
                            axs[itr_beta].legend()
                            axs[itr_beta].set_xlabel('x')
                            axs[itr_beta].set_ylabel('phi')
                        fig.suptitle(f'{title}_PatternFlat_error{errors2[1]}')
                        fig.savefig(f'figures/phis/phis_{savedata_pre}.pdf')
                        plt.close()

                        if(flag1==0 and flag2==1):
                            Lall[izp,izm]=Lbetas2[0]##length is of the patterned phase
                        if(flag1==1 and flag2==0):
                            Lall[izp,izm]=Lbetas2[1]

                        errorsall[izp,izm]=errors2[1]
                        errorsLbetaall[izp,izm]=errors2[3]
                        energyall[izp,izm]=total_energy2
                        finished[izp,izm]=1
                    elif(flag1==flag2 and flag1==1):###flat phases

                        if (np.abs(compart1max-compart2max)<different_phase_value):##Homogeneous
                            properties=['r-','b-','r--','b--','k:']
                            names=['p+','p-','e-','e+','S']
                            plt.figure()
                            for itr_beta in range(1):
                                xs=np.arange(128)*Lbetas2[itr_beta]/128
                                for itr_comp in range(5):
                                    plt.plot(xs,phis2[itr_comp,itr_beta],properties[itr_comp],label=str(names[itr_comp]))
                                plt.legend()
                                plt.xlabel('x')
                                plt.ylabel('phi')
                            plt.title(f'{title}_Homo_error{errors2[1]}')
                            plt.savefig(f'figures/phis/phis_{savedata_pre}.pdf')
                            plt.close()
                            Lall[izp,izm]=np.min(Lbetas2)
                            errorsall[izp,izm]=errors2[1]
                            errorsLbetaall[izp,izm]=errors2[3]
                            energyall[izp,izm]=total_energy2
                            finished[izp,izm]=3

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
                            Lall[izp,izm]=np.min(Lbetas2)
                            errorsall[izp,izm]=errors2[1]
                            errorsLbetaall[izp,izm]=errors2[3]
                            energyall[izp,izm]=total_energy2
                            finished[izp,izm]=2
                    elif(flag1==flag2 and flag1==0):###two patterns
                        if(np.abs(compart1max-compart2max)<different_phase_value and np.abs(compart1min-compart2min)<different_phase_value and np.abs(compart1amplitude-compart2amplitude)<different_phase_value):### same patterns in the two phases, so plot only one phase
                            index_beta=np.argmin(Lbetas2)
                            properties=['r-','b-','r--','b--','k:']
                            names=['p+','p-','e-','e+','S']
                            plt.figure()
                            for itr_beta in range(1):
                                xs=np.arange(128)*Lbetas2[index_beta]/128
                                for itr_comp in range(5):
                                    plt.plot(xs,phis2[itr_comp,index_beta],properties[itr_comp],label=str(names[itr_comp]))
                                plt.legend()
                            plt.title(f'{title}_Pattern_error{errors2[1]}')
                            plt.xlabel('x')
                            plt.ylabel('phi')
                            plt.savefig(f'figures/phis/phis_{savedata_pre}.pdf')
                            plt.close()
                            Lall[izp,izm]=Lbetas2[index_beta]
                            errorsall[izp,izm]=errors2[1]
                            errorsLbetaall[izp,izm]=errors2[3]
                            energyall[izp,izm]=total_energy2
                            finished[izp,izm]=0
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

                elif(total_energy2constantomega<=total_energy2 and total_energy2constantomega<=total_energy1):

                    ### rename to reuse the code above
                    Lbeats2=Lbetas2constantomega
                    omegas2=omegas2constantomega
                    phis2=phis2constantomega
                    psi2=psi2constantomega
                    Js2=Js2constantomega
                    total_energy2=total_energy2constantomega
                    errors2=errors2constantomega

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

                    if(flag1!=flag2):#### real mixed case, pattern+flat
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
                        fig.suptitle(f'{title}_PatternFlat_error{errors2[1]}')
                        fig.savefig(f'figures/phis/phis_{savedata_pre}.pdf')
                        plt.close()

                        if(flag1==0 and flag2==1):
                            Lall[izp,izm]=Lbetas2[0]##length is of the patterned phase
                        if(flag1==1 and flag2==0):
                            Lall[izp,izm]=Lbetas2[1]

                        errorsall[izp,izm]=errors2[1]
                        errorsLbetaall[izp,izm]=errors2[3]
                        energyall[izp,izm]=total_energy2
                        finished[izp,izm]=1
                    elif(flag1==flag2 and flag1==1):###flat phases

                        if (np.abs(compart1max-compart2max)<different_phase_value):##Homogeneous
                            properties=['r-','b-','r--','b--','k:']
                            names=['p+','p-','e-','e+','S']
                            plt.figure()
                            for itr_beta in range(1):
                                xs=np.arange(128)*Lbetas2[itr_beta]/128
                                for itr_comp in range(5):
                                    plt.plot(xs,phis2[itr_comp,itr_beta],properties[itr_comp],label=str(names[itr_comp]))
                                plt.legend()
                                plt.xlabel('x')
                                plt.ylabel('phi')
                            plt.title(f'{title}_Homo_error{errors2[1]}')
                            plt.savefig(f'figures/phis/phis_{savedata_pre}.pdf')
                            plt.close()
                            Lall[izp,izm]=np.min(Lbetas2)
                            errorsall[izp,izm]=errors2[1]
                            errorsLbetaall[izp,izm]=errors2[3]
                            energyall[izp,izm]=total_energy2
                            finished[izp,izm]=3

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
                            Lall[izp,izm]=np.min(Lbetas2)
                            errorsall[izp,izm]=errors2[1]
                            errorsLbetaall[izp,izm]=errors2[3]
                            energyall[izp,izm]=total_energy2
                            finished[izp,izm]=2
                    elif(flag1==flag2 and flag1==0):###two patterns
                        if(np.abs(compart1max-compart2max)<different_phase_value and np.abs(compart1min-compart2min)<different_phase_value and np.abs(compart1amplitude-compart2amplitude)<different_phase_value):### same patterns in the two phases, so plot only one phase
                            index_beta=np.argmin(Lbetas2)
                            properties=['r-','b-','r--','b--','k:']
                            names=['p+','p-','e-','e+','S']
                            plt.figure()
                            for itr_beta in range(1):
                                xs=np.arange(128)*Lbetas2[index_beta]/128
                                for itr_comp in range(5):
                                    plt.plot(xs,phis2[itr_comp,index_beta],properties[itr_comp],label=str(names[itr_comp]))
                                plt.legend()
                            plt.title(f'{title}_Pattern_error{errors2[1]}')
                            plt.xlabel('x')
                            plt.ylabel('phi')
                            plt.savefig(f'figures/phis/phis_{savedata_pre}.pdf')
                            plt.close()
                            Lall[izp,izm]=Lbetas2[index_beta]
                            errorsall[izp,izm]=errors2[1]
                            errorsLbetaall[izp,izm]=errors2[3]
                            energyall[izp,izm]=total_energy2
                            finished[izp,izm]=0
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

            except:
                print(f"error in {savedata_pre}")
    zpplot,zmplot=np.meshgrid(zps,-zms,indexing='ij')
    plt.figure()
    plt.pcolor(zpplot,zmplot,Lall)
    plt.xlabel('zp')
    plt.ylabel('-zm')
    plt.colorbar()
    plt.title('Length')
    plt.savefig(f'figures/lowestenergy_Lall_chi{chi}.pdf',dpi=300)
    plt.close()

    plt.figure()
    plt.pcolor(zpplot,zmplot,np.log10(errorsall))
    plt.xlabel('zp')
    plt.ylabel('-zm')
    plt.colorbar()
    plt.title('log10(error omega)')
    plt.savefig(f'figures/lowestenergy_erroromega_chi{chi}.pdf',dpi=300)
    plt.close()

    plt.figure()
    plt.pcolor(zpplot,zmplot,np.log10(errorsLbetaall))
    plt.xlabel('zp')
    plt.ylabel('-zm')
    plt.colorbar()
    plt.title('log10(error Lbeta)')
    plt.savefig(f'figures/lowestenergy_errorsLbeta_chi{chi}.pdf',dpi=300)
    plt.close()

    plt.figure()
    plt.pcolor(zpplot,zmplot,energyall)
    plt.xlabel('zp')
    plt.ylabel('-zm')
    plt.title('energyall')
    plt.colorbar()
    plt.savefig(f'figures/lowestenergy_energy_chi{chi}.pdf',dpi=300)
    plt.close()

    plt.figure()
    plt.pcolor(zpplot,zmplot,finished)
    plt.xlabel('zp')
    plt.ylabel('-zm')
    plt.title('types')
    plt.colorbar()
    plt.savefig(f'figures/lowestenergy_types_chi{chi}.pdf',dpi=300)
    plt.close()

    np.savetxt(f"Lall_chi{chi}.txt",Lall)
    np.savetxt(f"energyall_chi{chi}.txt",energyall)
    np.savetxt(f"types_chi{chi}.txt",finished)
