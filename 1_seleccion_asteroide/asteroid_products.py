import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from prettytable import PrettyTable
import pandas as pd


def asteroid_out(H_min,H_max,date_app_min,date_app_max,asteroid,asteroid_removed):
    '''
    This function generates the following products
    1.- Generation of csv ouput file
    2.- Generation of decision table file (simplified)
    3.- Generation of plots
    
    return none
    '''
    asteroid.reset_index(drop=True, inplace=True)

    # 1.-
    archive=f'asteroid_output_{H_min}_{H_max}_{date_app_min}_{date_app_max}.csv'
    asteroid.to_csv(archive, index=False)
     
    # 2.- 
    summary_table = PrettyTable(['Asteroid', 'Accesibility [m/s]','Orbit Uncertainty','Synodic Period [y]','Additional Info'])
    for n_ast in range(len(asteroid)):
        
        additional_info=[]
        if asteroid.loc[n_ast,'SMASS taxonomy'] is not None: additional_info.append(f'SMASSII Taxonomy Known: {asteroid.loc[n_ast,"SMASS taxonomy"]}')
        if asteroid.loc[n_ast,'Spin period'] is not None: additional_info.append(f'Spin rate known {round(1/(float(asteroid.loc[n_ast,"Spin period"])*60),4)} rpm')
        if asteroid.loc[n_ast,'Satellites']==1: additional_info.append('Secondary body')
        if asteroid.loc[n_ast,'approaches']>=3: additional_info.append(f"{asteroid.loc[n_ast,'approaches']} close approaches from {date_app_min} to {date_app_max}")
        if asteroid.loc[n_ast,'is_NHATS']==True: additional_info.append('Included in NHATS database')
        if asteroid.loc[n_ast,'is_geometry']==True: additional_info.append('Geometry model available' )
            
        summary_table.add_row([asteroid.loc[n_ast,'ID'], round(asteroid.loc[n_ast,'delta_v_tot'],2), asteroid.loc[n_ast,'condition_code'], round(asteroid.loc[n_ast,'period_sin'],2),additional_info])
    print(summary_table)
    table = f'Decision_Table_{H_min}_{H_max}_{date_app_min}_{date_app_max}.txt'
    with open(table, 'w') as file:
        file.write(summary_table.get_string())
    
    # 3.-
    #Close approachers data
    a=np.zeros(len(asteroid))
    e=np.zeros(len(asteroid))
    i=np.zeros(len(asteroid))
    u_inf=np.zeros(len(asteroid))
    delta_v=np.zeros(len(asteroid))

    for n_ast in range(0,len(asteroid)):
        a[n_ast]=asteroid.loc[n_ast,'a'] #a
        e[n_ast]=asteroid.loc[n_ast,'e'] #e
        i[n_ast]=asteroid.loc[n_ast,'i'] #i
        u_inf[n_ast]=np.sqrt(abs(3-(1/a[n_ast]+2*np.sqrt(a[n_ast]*(1-e[n_ast]**2))*np.cos(i[n_ast]*np.pi/180)))) 
        delta_v[n_ast]=asteroid.loc[n_ast,'delta_v_tot'] #delta_v_tot

    #Removed asteroid data
    a_others=np.zeros(len(asteroid_removed))
    e_others=np.zeros(len(asteroid_removed))
    i_others=np.zeros(len(asteroid_removed))
    for n_ast in range(0,len(asteroid_removed)):
        a_others[n_ast]=asteroid_removed.loc[n_ast,'a'] #a
        e_others[n_ast]=asteroid_removed.loc[n_ast,'e'] #e
        i_others[n_ast]=asteroid_removed.loc[n_ast,'i'] #i

    #Plots
    fig=plt.figure(figsize=(15,15))

    ax1=plt.subplot(2,2,1) # plot a vs e
    ax1.plot(a_others,e_others,marker='+',markersize=4,color='b',linestyle='',label='NEAs analyzed')
    ax1.plot(a,e,marker='o',markeredgecolor='k',markersize=6,color='r',linestyle='',label='NEAs w. close approach')
    ax1.set_xlabel("a [au]")
    ax1.set_ylabel("e ")
    ax1.set_xlim(0,4)
    ax1.set_ylim(0,1)
    ax1.grid(True)
    ax1.legend()

    ax2=plt.subplot(2,2,2) # plot a vs i
    ax2.plot(a_others,i_others,marker='+',markersize=4,color='b',linestyle='',label='NEAs analized')
    ax2.plot(a,i,marker='o',markeredgecolor='k',markersize=6 ,color='r',linestyle='',label='NEAs w. close approach')
    ax2.set_xlabel("a [au]")
    ax2.set_ylabel("i [deg]")
    ax2.set_xlim(0,4)
    ax2.set_ylim(0,90)
    ax2.grid(True)
    ax2.legend()

    ax3=plt.subplot(2,2,3) # plot u_inf vs delta_v
    ax3.plot(u_inf,delta_v,marker='o',markeredgecolor='k',color='r',markersize=8 ,linestyle='')
    ax3.set_xlabel("U_inf")
    ax3.set_ylabel("Delta_v [m/s]")
    ax3.set_xlim(0,1.4)
    ax3.set_ylim(bottom=0)
    ax3.grid(True)

    # #d Accesibility
    ax4=fig.add_subplot(2,2,4,projection='3d') # plot u_inf vs delta_v
    norm = plt.Normalize(vmin=np.min(delta_v), vmax=np.max(delta_v))
    sc=ax4.scatter(a,e,i,c=delta_v,cmap='cool',alpha=1,norm=norm)
    #Projections
    ax4.scatter(a, e, np.zeros_like(i), color='gray',  marker='x')
    ax4.scatter(a, np.ones_like(e), i, color='gray',  marker='x')
    ax4.scatter(np.zeros_like(a), e, i, color='gray',  marker='x')
    ax4.set_xlabel('a [au]')
    ax4.set_ylabel('e')
    ax4.set_zlabel('i [deg]')
    ax4.set_xlim([0,2])
    ax4.set_ylim([0,1])
    ax4.set_zlim([0,30])
    cbar = plt.colorbar(sc, ax=ax4)
    cbar.set_label('delta_v [m/s]')

    plt.tight_layout()
    image_asteroid=f'Accesibility_{H_min}_{H_max}_{date_app_min}_{date_app_max}.png'
    plt.savefig(image_asteroid)
    
    return