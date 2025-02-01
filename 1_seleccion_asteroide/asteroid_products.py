import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from prettytable import PrettyTable


def asteroid_out(H_min,H_max,date_app_min,date_app_max,asteroid,asteroid_removed):
    '''
    This function generates the following products
    1.- Generation of csv ouput file
    2.- Generation of decision table file (simplified)
    3.- Generation of plots
    
    return none
    '''
    
    # 1.-
    archive=f'asteroid_output_{H_min}_{H_max}_{date_app_min}_{date_app_max}.csv'
    with open(archive, mode='w',newline='') as csvfile:
       writer=csv.writer(csvfile)
       writer.writerow(["0 ID, 1 H, 2 e, 3 a, 4 Perigee, 5 i, 6 RAAN, 7 arg_perig, 8 Apogee, 9 period [y], \
    10 MOID, 11 condition_code, 12 SMASS taxonomy, 13 Spin period,14 Satellites, 15 PHA,\
    16 del_v_tot, 17 period_sin_y, 18 g_2^2, 19 Familiar, 20 Approaches, 21 NHATs, 22 Geometry"])
       writer.writerows(asteroid)
     
    # 2.- 
    summary_table = PrettyTable(['Asteroid', 'Accesibility [m/s]','Orbit Uncertainty','Synodic Period [y]','Additional Info'])
    for n_ast in range(len(asteroid)):
        
        additional_info=[]
        if asteroid[n_ast][12] is not None: additional_info.append(f'SMASSII Taxonomy Known: {asteroid[n_ast][12]}')
        if asteroid[n_ast][13] is not None: additional_info.append(f'Spin rate known {round(1/(float(asteroid[n_ast][13])*60),4)} rpm')
        if asteroid[n_ast][14]==1: additional_info.append('Secondary body')
        if asteroid[n_ast][20]>=3: additional_info.append(f'{asteroid[n_ast][20]} close approaches from {date_app_min} to {date_app_max}')
        if asteroid[n_ast][21]==True: additional_info.append('Included in NHATS database')
        if asteroid[n_ast][22]==True: additional_info.append('Geometry model available' )
            
        summary_table.add_row([asteroid[n_ast][0], round(asteroid[n_ast][16],2), asteroid[n_ast][11], round(asteroid[n_ast][17],2),additional_info])
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
        a[n_ast]=asteroid[n_ast][3] #a
        e[n_ast]=asteroid[n_ast][2] #e
        i[n_ast]=asteroid[n_ast][5] #i
        u_inf[n_ast]=np.sqrt(abs(3-(1/a[n_ast]+2*np.sqrt(a[n_ast]*(1-e[n_ast]**2))*np.cos(i[n_ast]*np.pi/180)))) 
        delta_v[n_ast]=asteroid[n_ast][16] #delta_v_tot

    #Removed asteroid data
    a_others=np.zeros(len(asteroid_removed))
    e_others=np.zeros(len(asteroid_removed))
    i_others=np.zeros(len(asteroid_removed))
    for n_ast in range(0,len(asteroid_removed)):
        a_others[n_ast]=asteroid_removed[n_ast][3] #a
        e_others[n_ast]=asteroid_removed[n_ast][2] #e
        i_others[n_ast]=asteroid_removed[n_ast][5] #i

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