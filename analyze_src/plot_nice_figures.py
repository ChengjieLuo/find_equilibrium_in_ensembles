import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcol
from analyze_src.ufun import *




fig,axs=plt.subplots(4,4,figsize=(16,12),sharey=True,sharex=True)
fig2,axs2=plt.subplots()
for ichi, chi in enumerate([-7.0, -5.0, -4.0, -3.0]):
    zps=np.linspace(0.1,0.9,30)
    zms=np.linspace(-0.1,-0.9,30)
    Lall=np.loadtxt(f"Lall_chi{chi}.txt")
    energyall=np.loadtxt(f"energyall_chi{chi}.txt")
    types=np.loadtxt(f"types_chi{chi}.txt")
    abscharges=np.loadtxt(f"abscharges_chi{chi}.txt")



    
    zpplot,zmplot=np.meshgrid(zps,-zms,indexing='ij')


    bounds=[0,1,2,3]
    cmap=cm.get_cmap("Paired",lut=len(bounds)+1)
    cmap_bounds = np.arange(len(bounds)+1) - 0.5
    norm = mcol.BoundaryNorm(cmap_bounds,cmap.N)
    pc0=axs[0,ichi].pcolormesh(zpplot,zmplot,types,cmap=cmap,norm=norm)
    # axs[0,ichi].set_xlabel('zp')
    axs[0,ichi].set_ylabel(r'$-Z_a$')
    cbar=plt.colorbar(pc0,ax=axs[0,ichi],ticks=[0,1,2,3])
    cbar.set_ticklabels(["P","P+H","H+H","H"])
    axs[0,ichi].set_title(f'chi={chi}',fontsize=15)

    if(ichi==0):
        axs[0,ichi].text(0.2,0.8,'P',fontsize=15)
        axs[0,ichi].text(0.8,0.2,'P',fontsize=15)
        axs[0,ichi].text(0.6,0.6,'H+H',fontsize=15)
    if(ichi==1):
        axs[0,ichi].text(0.2,0.8,'P',fontsize=15)
        axs[0,ichi].text(0.8,0.2,'P',fontsize=15)
        axs[0,ichi].text(0.6,0.6,'H+H',fontsize=15)
    if(ichi==2):
        axs[0,ichi].text(0.2,0.8,'H',fontsize=15)
        axs[0,ichi].text(0.8,0.2,'H',fontsize=15)
        axs[0,ichi].text(0.6,0.6,'H+H',fontsize=15)
        axs[0,ichi].text(0.1,0.25,'P',fontsize=15)
        axs[0,ichi].text(0.25,0.08,'P',fontsize=15)
    if(ichi==3):
        axs[0,ichi].text(0.5,0.5,'H',fontsize=15)

    # plt.title('types')
    # plt.colorbar()
    # plt.savefig(f'figures/lowestenergy_types.pdf',dpi=300)

    Lallplot=Lall.copy()
    Lallplot[types==2]=np.nan
    Lallplot[types==3]=np.nan
    cmap=cm.get_cmap('plasma').copy()
    cmap.set_bad(color = 'gray')
    pc1=axs[1,ichi].pcolormesh(zpplot,zmplot,Lallplot,vmin=20,vmax=200,cmap=cmap)

    if(ichi==0 or ichi==1):
        fig1,axs1=plt.subplots()
        contour=axs1.contour(zpplot,zmplot,Lall,20,vmin=30,vmax=210)
        plt.colorbar(contour,ax=axs1)
        fig1.savefig(f'contour_L_chi{chi}.pdf')
    
    if(ichi==0):
        contour=axs2.contour(zpplot,zmplot,Lall,20,vmin=30,vmax=210,colors='r')
        plt.colorbar(contour,ax=axs2)

        Lall7=Lall.copy()
        
    if(ichi==1):
        contour=axs2.contour(zpplot,zmplot,Lall,20,vmin=30,vmax=210,colors='g')
        plt.colorbar(contour,ax=axs2)
        fig2.savefig(f'contour_L.pdf')

        ff,aa=plt.subplots(1,1)
        cmap=cm.get_cmap('bwr').copy()
        cmap.set_bad('gray')
        diffL=Lall7-Lall
        diffL[types==2]=np.nan
        diffL[types==3]=np.nan
        paa=aa.pcolormesh(zpplot,zmplot,diffL,vmin=-20,vmax=20,cmap=cmap)
        cb=plt.colorbar(paa,ax=aa)
        cb.set_label(r'$L(\chi=-7)-L(\chi=-5)$')
        aa.set_xlabel(r'$Z_c$')
        aa.set_ylabel(r'$-Z_a$')
        ff.savefig('diff_L.pdf')
    

    # axs[1,ichi].set_xlabel('zp')
    axs[1,ichi].set_ylabel(r'$-Z_a$')
    plt.colorbar(pc1,ax=axs[1,ichi])
    # plt.colorbar()
    # plt.title('Length')
    # plt.savefig(f'figures/lowestenergy_Lall.pdf',dpi=300)

    
    pc2=axs[2,ichi].pcolormesh(zpplot,zmplot,energyall,cmap='cool')
    axs[2,ichi].set_xlabel(r'$Z_c$')
    axs[2,ichi].set_ylabel(r'$-Z_a$')
    plt.colorbar(pc2,ax=axs[2,ichi])
    # plt.title('energyall')
    # plt.colorbar()
    # plt.savefig(f'figures/lowestenergy_energy.pdf',dpi=300)

    pc3=axs[3,ichi].pcolormesh(zpplot,zmplot,abscharges,cmap='viridis',vmin=abscharges.min()-1e-3,vmax=abscharges.max()+1e-3)
    # print(abscharges)
    axs[3,ichi].set_xlabel(r'$Z_c$')
    axs[3,ichi].set_ylabel(r'$-Z_a$')
    plt.colorbar(pc3,ax=axs[3,ichi])
    # plt.title('energyall')
    # plt.colorbar()
    # plt.savefig(f'figures/lowestenergy_energy.pdf',dpi=300)

    if(ichi==0):
        axs[0,0].text(-0.3,0.4,'Type',rotation='vertical',fontsize=15)
        axs[1,0].text(-0.3,0.4,'Length',rotation='vertical',fontsize=15)
        axs[2,0].text(-0.3,0.2,'Energy density',rotation='vertical',fontsize=15)
        axs[3,0].text(-0.3,0.2,'Absolute charge',rotation='vertical',fontsize=15)
fig.savefig('all.pdf',dpi=300,bbox_inches='tight')
plt.close()




######### plot contributions

fig,axs=plt.subplots(8,4,figsize=(16,20),sharey=True,sharex=True)
fig2,axs2=plt.subplots()
for ichi, chi in enumerate([-7.0, -5.0, -4.0, -3.0]):
    zps=np.linspace(0.1,0.9,30)
    zms=np.linspace(-0.1,-0.9,30)
    Lall=np.loadtxt(f"Lall_chi{chi}.txt")
    energyall=np.loadtxt(f"energyall_chi{chi}.txt")
    types=np.loadtxt(f"types_chi{chi}.txt")
    abscharges=np.loadtxt(f"abscharges_chi{chi}.txt")
    entropy_all=np.loadtxt(f"entropy_chi{chi}.txt")
    chi_all=np.loadtxt(f"chi_chi{chi}.txt")
    psizsphi_all=np.loadtxt(f"psizsphi_chi{chi}.txt")
    kappa_all=np.loadtxt(f"kappa_chi{chi}.txt")
    nablapsi_all=np.loadtxt(f"nablapsi_chi{chi}.txt")
    zero_all=np.loadtxt(f"zero_chi{chi}.txt")
    totale_all=np.loadtxt(f"totale_chi{chi}.txt")



    
    zpplot,zmplot=np.meshgrid(zps,-zms,indexing='ij')


    # bounds=[0,1,2,3]
    # cmap=cm.get_cmap("Paired",lut=len(bounds)+1)
    # cmap_bounds = np.arange(len(bounds)+1) - 0.5
    # norm = mcol.BoundaryNorm(cmap_bounds,cmap.N)
    # pc0=axs[0,ichi].pcolormesh(zpplot,zmplot,types,cmap=cmap,norm=norm)
    # # axs[0,ichi].set_xlabel('zp')
    # axs[0,ichi].set_ylabel(r'$-Z_a$')
    # cbar=plt.colorbar(pc0,ax=axs[0,ichi],ticks=[0,1,2,3])
    # cbar.set_ticklabels(["P","P+H","H+H","H"])
    # axs[0,ichi].set_title(f'chi={chi}',fontsize=15)

    # if(ichi==0):
    #     axs[0,ichi].text(0.2,0.8,'P',fontsize=15)
    #     axs[0,ichi].text(0.8,0.2,'P',fontsize=15)
    #     axs[0,ichi].text(0.6,0.6,'H+H',fontsize=15)
    # if(ichi==1):
    #     axs[0,ichi].text(0.2,0.8,'P',fontsize=15)
    #     axs[0,ichi].text(0.8,0.2,'P',fontsize=15)
    #     axs[0,ichi].text(0.6,0.6,'H+H',fontsize=15)
    # if(ichi==2):
    #     axs[0,ichi].text(0.2,0.8,'H',fontsize=15)
    #     axs[0,ichi].text(0.8,0.2,'H',fontsize=15)
    #     axs[0,ichi].text(0.6,0.6,'H+H',fontsize=15)
    #     axs[0,ichi].text(0.1,0.25,'P',fontsize=15)
    #     axs[0,ichi].text(0.25,0.08,'P',fontsize=15)
    # if(ichi==3):
    #     axs[0,ichi].text(0.5,0.5,'H',fontsize=15)

    # plt.title('types')
    # plt.colorbar()
    # plt.savefig(f'figures/lowestenergy_types.pdf',dpi=300)

    # Lallplot=Lall.copy()
    # Lallplot[types==2]=np.nan
    # Lallplot[types==3]=np.nan
    # cmap=cm.get_cmap('plasma').copy()
    # cmap.set_bad(color = 'gray')
    # pc1=axs[1,ichi].pcolormesh(zpplot,zmplot,Lallplot,vmin=20,vmax=200,cmap=cmap)

    # if(ichi==0 or ichi==1):
    #     fig1,axs1=plt.subplots()
    #     contour=axs1.contour(zpplot,zmplot,Lall,20,vmin=30,vmax=210)
    #     plt.colorbar(contour,ax=axs1)
    #     fig1.savefig(f'contour_L_chi{chi}.pdf')
    
    # if(ichi==0):
    #     contour=axs2.contour(zpplot,zmplot,Lall,20,vmin=30,vmax=210,colors='r')
    #     plt.colorbar(contour,ax=axs2)

    #     Lall7=Lall.copy()
        
    # if(ichi==1):
    #     contour=axs2.contour(zpplot,zmplot,Lall,20,vmin=30,vmax=210,colors='g')
    #     plt.colorbar(contour,ax=axs2)
    #     fig2.savefig(f'contour_L.pdf')

    #     ff,aa=plt.subplots(1,1)
    #     cmap=cm.get_cmap('bwr').copy()
    #     cmap.set_bad('gray')
    #     diffL=Lall7-Lall
    #     diffL[types==2]=np.nan
    #     diffL[types==3]=np.nan
    #     paa=aa.pcolormesh(zpplot,zmplot,diffL,vmin=-20,vmax=20,cmap=cmap)
    #     cb=plt.colorbar(paa,ax=aa)
    #     cb.set_label(r'$L(\chi=-7)-L(\chi=-5)$')
    #     aa.set_xlabel(r'$Z_c$')
    #     aa.set_ylabel(r'$-Z_a$')
    #     ff.savefig('diff_L.pdf')
    

    # # axs[1,ichi].set_xlabel('zp')
    # axs[1,ichi].set_ylabel(r'$-Z_a$')
    # plt.colorbar(pc1,ax=axs[1,ichi])
    # # plt.colorbar()
    # # plt.title('Length')
    # # plt.savefig(f'figures/lowestenergy_Lall.pdf',dpi=300)

    
    # pc2=axs[2,ichi].pcolormesh(zpplot,zmplot,energyall,cmap='cool')
    # axs[2,ichi].set_xlabel(r'$Z_c$')
    # axs[2,ichi].set_ylabel(r'$-Z_a$')
    # plt.colorbar(pc2,ax=axs[2,ichi])
    # # plt.title('energyall')
    # # plt.colorbar()
    # # plt.savefig(f'figures/lowestenergy_energy.pdf',dpi=300)

    # pc3=axs[3,ichi].pcolormesh(zpplot,zmplot,abscharges,cmap='viridis',vmin=abscharges.min()-1e-3,vmax=abscharges.max()+1e-3)
    # # print(abscharges)
    # axs[3,ichi].set_xlabel(r'$Z_c$')
    # axs[3,ichi].set_ylabel(r'$-Z_a$')
    # plt.colorbar(pc3,ax=axs[3,ichi])
    # # plt.title('energyall')
    # # plt.colorbar()
    # # plt.savefig(f'figures/lowestenergy_energy.pdf',dpi=300)

    # if(ichi==0):
    #     axs[0,0].text(-0.3,0.4,'Type',rotation='vertical',fontsize=15)
    #     axs[1,0].text(-0.3,0.4,'Length',rotation='vertical',fontsize=15)
    #     axs[2,0].text(-0.3,0.2,'Energy density',rotation='vertical',fontsize=15)
    #     axs[3,0].text(-0.3,0.2,'Absolute charge',rotation='vertical',fontsize=15)
    axs[0,ichi].set_title(f'chi={chi}',fontsize=15)

    ax=axs[0,ichi]
    pc=ax.pcolormesh(zpplot,zmplot,totale_all,vmin=totale_all.min()-1e-3,vmax=totale_all.max()+1e-3,cmap='viridis')
    ax.set_xlabel(r'$Z_c$')
    ax.set_ylabel(r'$-Z_a$')
    plt.colorbar(pc,ax=ax)

    ax=axs[1,ichi]
    pc=ax.pcolormesh(zpplot,zmplot,entropy_all,vmin=entropy_all.min(),vmax=entropy_all.max(),cmap='viridis')
    ax.set_xlabel(r'$Z_c$')
    ax.set_ylabel(r'$-Z_a$')
    plt.colorbar(pc,ax=ax)

    ax=axs[2,ichi]
    pc=ax.pcolormesh(zpplot,zmplot,chi_all,vmin=chi_all.min()-1e-3,vmax=chi_all.max()+1e-3,cmap='viridis')
    ax.set_xlabel(r'$Z_c$')
    ax.set_ylabel(r'$-Z_a$')
    plt.colorbar(pc,ax=ax)

    ax=axs[3,ichi]
    pc=ax.pcolormesh(zpplot,zmplot,psizsphi_all,vmin=psizsphi_all.min()-1e-5,vmax=psizsphi_all.max()+1e-5,cmap='viridis')
    ax.set_xlabel(r'$Z_c$')
    ax.set_ylabel(r'$-Z_a$')
    plt.colorbar(pc,ax=ax)

    ax=axs[4,ichi]
    pc=ax.pcolormesh(zpplot,zmplot,kappa_all,vmin=kappa_all.min(),vmax=kappa_all.max(),cmap='viridis')
    ax.set_xlabel(r'$Z_c$')
    ax.set_ylabel(r'$-Z_a$')
    plt.colorbar(pc,ax=ax)

    ax=axs[5,ichi]
    pc=ax.pcolormesh(zpplot,zmplot,nablapsi_all,vmin=nablapsi_all.min()-1e-5,vmax=nablapsi_all.max()+1e-5,cmap='viridis')
    ax.set_xlabel(r'$Z_c$')
    ax.set_ylabel(r'$-Z_a$')
    plt.colorbar(pc,ax=ax)

    ax=axs[6,ichi]
    pc=ax.pcolormesh(zpplot,zmplot,zero_all,vmin=zero_all.min(),vmax=zero_all.max(),cmap='viridis')
    ax.set_xlabel(r'$Z_c$')
    ax.set_ylabel(r'$-Z_a$')
    plt.colorbar(pc,ax=ax)

    ax=axs[7,ichi]
    pc=ax.pcolormesh(zpplot,zmplot,psizsphi_all+nablapsi_all,vmin=(psizsphi_all+nablapsi_all).min()-1e-5,vmax=(psizsphi_all+nablapsi_all).max()+1e-5,cmap='viridis')
    ax.set_xlabel(r'$Z_c$')
    ax.set_ylabel(r'$-Z_a$')
    plt.colorbar(pc,ax=ax)

    if(ichi==0):
        axs[0,0].text(-0.3,0.4,'total',rotation='vertical',fontsize=15)
        axs[1,0].text(-0.3,0.4,'entropy',rotation='vertical',fontsize=15)
        axs[2,0].text(-0.3,0.4,'chi',rotation='vertical',fontsize=15)
        axs[3,0].text(-0.3,0.2,'psi zs phi',rotation='vertical',fontsize=15)
        axs[4,0].text(-0.3,0.2,'kappa',rotation='vertical',fontsize=15)
        axs[5,0].text(-0.3,0.2,'nabla psi',rotation='vertical',fontsize=15)
        axs[6,0].text(-0.3,0.2,'phimean',rotation='vertical',fontsize=15)
        axs[7,0].text(-0.3,0.2,'total electric',rotation='vertical',fontsize=15)

fig.savefig('energy_contributions.pdf',dpi=300,bbox_inches='tight')













