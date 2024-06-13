import os
import numpy as np

# num_comps=5
# num_beta=1
# num_coord=128
# ### num_comps num_beta num_coord phi_means chi Ls zs kappas v [file filenameomega || random amplitude] steps_inner1 steps_inner2 acceptance_omega acceptance_J acceptance_Lbeta acceptance_zeta [file filenameJs||random amplitude] [file filenameLbetas||equal valueL] [file filenamezetas||equal valuezeta] flag_C C [file filenameps||equal valuep] flag_zetas flag_ps savedata_pre flag_save_separate
# ./solve_gibbs_with_input.out $num_comps $num_beta $num_coord 0.1 0.1 0.06 0.015 0.725 -5.0 10.0 10.0 1.0 1.0 1.0 0.6 -0.15 -1.0 1.0 0.0 10.0 10.0 10.0 10.0 0.0 100 random 5.0 100001 300001 0.001 0.001 10.0 10.0 random 0.0 equal 30 equal 0.0 true 100 equal 1.0 true false random1_0.6_-0.15 true

zps=np.linspace(0.1,0.9,30)
zms=np.linspace(-0.1,-0.9,30)
print(f"{zps=}")
print(f"{zms=}")

np.set_printoptions(precision=8)
for chi in [-7.0, -6.0, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0]:
    for zp in zps:
        for zm in zms:
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

            print(f"{savedata_pre=}")

            savedata_folder='./data_1_comp/'


            command=f"OMP_NUM_THREADS=1; qsub -S /bin/bash -b y -o out/out{savedata_pre}.txt -e err/err{savedata_pre}.txt -q teutates.q -N recon0 -cwd ./solve_gibbs_with_input.out {num_comps} {num_beta} {num_coord} {phi_means_str} {chi} {Ls_str} {zs_str} {kappas_str} {v} {omega_type} {omega_value} {steps_inner1} {steps_inner2} {acceptance_omega} {acceptance_J} {acceptance_Lbeta} {acceptance_zeta} {Js_type} {Js_value} {Lbeta_type} {Lbeta_value} {zetas_type} {zetas_value} {flag_C} {C} {ps_type} {ps_value} {flag_zetas} {flag_ps} {savedata_folder+savedata_pre} {flag_save_separate}"

            # command=f"./solve_gibbs_with_input.out {num_comps} {num_beta} {num_coord} {phi_means_str} {chi} {Ls_str} {zs_str} {kappas_str} {v} {omega_type} {omega_value} {steps_inner1} {steps_inner2} {acceptance_omega} {acceptance_J} {acceptance_Lbeta} {acceptance_zeta} {Js_type} {Js_value} {Lbeta_type} {Lbeta_value} {zetas_type} {zetas_value} {flag_C} {C} {ps_type} {ps_value} {flag_zetas} {flag_ps} {savedata_folder+savedata_pre} {flag_save_separate}"

            print("command=")
            print(command)


            # command=f"./solve_gibbs_with_input.out {num_comps} {num_beta} {num_coord} {phi_means_str} {chi} {Ls_str} {zs_str} {kappas_str} {v} {omega_type} {omega_value} {steps_inner1} {steps_inner2} {acceptance_omega} {acceptance_J} {acceptance_Lbeta} {acceptance_zeta} {Js_type} {Js_value} {Lbeta_type} {Lbeta_value} {zetas_type} {zetas_value} {flag_C} {C} {ps_type} {ps_value} {flag_zetas} {flag_ps} {savedata_folder+savedata_pre} {flag_save_separate}"# > {savedata_folder+savedata_pre}.log"

            

            os.system(command)