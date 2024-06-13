import os
import numpy as np


flag_one_compartment=False
flag_two_compartment=False
flag_four_compartment=False
flag_eight_compartment=True

# zps=np.linspace(0.1,0.9,30)
# zms=np.linspace(-0.1,-0.9,30)
# print(f"{zps=}")
# print(f"{zms=}")

np.set_printoptions(precision=8)
for chi in [-5.0]:
    for zp in [0.8]:
        # for zm in [-0.2]:
        for zm in [-0.6]:
            for num_beta in [1,2,4,8]:

                if(num_beta==1 and flag_one_compartment==False):
                    continue

                if(num_beta==2 and flag_two_compartment==False):
                    continue

                if(num_beta==4 and flag_four_compartment==False):
                    continue

                if(num_beta==8 and flag_eight_compartment==False):
                    continue


                num_comps=5
                # num_beta=2
                num_coord=128

                # chi=-5
                chis=np.zeros((num_comps,num_comps))
                chis[0,1]=chis[1,0]=chi
                chis-=np.min(chis)

                

                Ls=np.array([10.,10.,1.,1.,1.])

                zs=np.array([zp ,zm, -1.0 ,1.0, 0.0])

                #### phi_means for ions is controlled by neutral charge 
                phi_means=np.array([0.1,0.1,0,0,0])
                phi_means[2]=-zs[0]*phi_means[0]/zs[2]
                phi_means[3]=-zs[1]*phi_means[1]/zs[3]
                phi_means[4]=1-phi_means[:4].sum()
                print(f"total charge={np.sum(zs*phi_means)}")


                kappas=np.array([1.,1.,1.,1.,0.])*10
                v=100
                # [file filenameomega || random amplitude] 
                omega_type="random"
                omega_value=5.0
                steps_inner1= 100001
                steps_inner2= 1000001
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

                threshold_incomp=1e-9
                threshold_omega=1e-7
                threshold_J=1e-6
                threshold_Lbeta=1e-8
                threshold_zeta=1e-7

                ###convert vectors to strings
                chis_str=" "
                phi_means_str=" "
                Ls_str=" "
                zs_str=" "
                kappas_str=" "
                for i in range(num_comps):
                    phi_means_str=phi_means_str+f"{phi_means[i]:.15f} "
                    Ls_str=Ls_str+f"{Ls[i]:.15f} "
                    zs_str=zs_str+f"{zs[i]:.15f} "
                    kappas_str=kappas_str+f"{kappas[i]:.15f} "
                for i in range(num_comps):
                    for j in range(num_comps):
                        chis_str=chis_str+f"{chis[i,j]:.15f} "


                savedata_pre= f"{num_comps}_{num_beta}_{num_coord}_{phi_means_str}_{chis_str}_{Ls_str}_{zs_str}_{kappas_str}_{v}_{omega_type}_{omega_value}_{steps_inner1}_{steps_inner2}_{acceptance_omega}_{acceptance_J}_{acceptance_Lbeta}_{acceptance_zeta}_{Js_type}_{Js_value}_{Lbeta_type}_{Lbeta_value}_{zetas_type}_{zetas_value}_{flag_C}_{C}_{ps_type}_{ps_value}_{flag_zetas}_{flag_ps}_{flag_save_separate}_{threshold_incomp}_{threshold_omega}_{threshold_J}_{threshold_Lbeta}_{threshold_zeta}"

                # import re
                # savedata_pre=re.sub(' +','_',savedata_pre)

                savedata_pre= f"zs{zs_str}_chi{chi}_numbeta{num_beta}"
                savedata_pre='_'.join(savedata_pre.split(' '))

                print(f"{savedata_pre=}")

                savedata_folder='./data_example/'


                # command=f"OMP_NUM_THREADS=1; qsub -S /bin/bash -b y -o out/out{savedata_pre}.txt -e err/err{savedata_pre}.txt -q teutates.q -N recon0 -cwd {exe} {num_comps} {num_beta} {num_coord} {phi_means_str} {chis_str} {Ls_str} {zs_str} {kappas_str} {v} {omega_type} {omega_value} {steps_inner1} {steps_inner2} {acceptance_omega} {acceptance_J} {acceptance_Lbeta} {acceptance_zeta} {Js_type} {Js_value} {Lbeta_type} {Lbeta_value} {zetas_type} {zetas_value} {flag_C} {C} {ps_type} {ps_value} {flag_zetas} {flag_ps} {savedata_folder+savedata_pre} {flag_save_separate} {threshold_incomp} {threshold_omega} {threshold_J} {threshold_Lbeta} {threshold_zeta}"

                exe='../src/solve_gibbs_with_input.out'
                command=f"{exe} {num_comps} {num_beta} {num_coord} {phi_means_str} {chis_str} {Ls_str} {zs_str} {kappas_str} {v} {omega_type} {omega_value} {steps_inner1} {steps_inner2} {acceptance_omega} {acceptance_J} {acceptance_Lbeta} {acceptance_zeta} {Js_type} {Js_value} {Lbeta_type} {Lbeta_value} {zetas_type} {zetas_value} {flag_C} {C} {ps_type} {ps_value} {flag_zetas} {flag_ps} {savedata_folder+savedata_pre} {flag_save_separate} {threshold_incomp} {threshold_omega} {threshold_J} {threshold_Lbeta} {threshold_zeta}"

                print("command=")
                print(command)


                # command=f"./solve_gibbs_with_input.out {num_comps} {num_beta} {num_coord} {phi_means_str} {chi} {Ls_str} {zs_str} {kappas_str} {v} {omega_type} {omega_value} {steps_inner1} {steps_inner2} {acceptance_omega} {acceptance_J} {acceptance_Lbeta} {acceptance_zeta} {Js_type} {Js_value} {Lbeta_type} {Lbeta_value} {zetas_type} {zetas_value} {flag_C} {C} {ps_type} {ps_value} {flag_zetas} {flag_ps} {savedata_folder+savedata_pre} {flag_save_separate}"# > {savedata_folder+savedata_pre}.log"

                os.system(command)


