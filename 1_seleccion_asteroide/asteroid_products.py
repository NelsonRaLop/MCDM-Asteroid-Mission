import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from prettytable import PrettyTable

def asteroid_out(H_min,H_max,date_app_min,date_app_max,asteroid,asteroid_removed):
    '''
    This function generates the following products
    1.- Generation of csv ouput file
    2.- Generation of decision table file (simplified): Asteroid', 'Accesibility [m/s]',
        'Potential backup asteroids','Orbit Uncertainty','Synodic Period [y]','Spin Rate [rpm]','Additional Info' 
    3.- Generation of plots
    
    return none
    '''

    # 1.-
    archive=f'asteroid_output_{H_min}_{H_max}_{date_app_min}_{date_app_max}.csv'
    asteroid.to_csv(archive, index=False)
     
    # 2.- 
    summary_table = PrettyTable(['Asteroid', 'Accesibility [m/s]','Potential backup asteroids','Orbit Uncertainty','Synodic Period [y]','Spin Rate [rpm]','Additional Info'])
    for n_ast in range(len(asteroid)):
        
        # Spin period and backup asteroid required 'a priori' (not shown in decision table but computed)
        if asteroid.loc[n_ast,'Spin period'] is None or asteroid.loc[n_ast,'n_backup']==0: continue 
            
        else: 
            familiar_value = str(asteroid.loc[n_ast, 'familiar'])
            match = re.search(r'\((.*?)\)', familiar_value)
            familiar_value = match.group(1)

            asteroid.loc[n_ast,'Spin period']=round(1/(float(asteroid.loc[n_ast,"Spin period"])*60),4)
            additional_info=[]
            if asteroid.loc[n_ast,'SMASS taxonomy'] is not None: additional_info.append(f'SMASSII Taxonomy Known: {asteroid.loc[n_ast,"SMASS taxonomy"]}')
            if asteroid.loc[n_ast,'Satellites']==1: additional_info.append('Secondary body')
            if asteroid.loc[n_ast,'approaches']>=3: additional_info.append(f"{asteroid.loc[n_ast,'approaches']} close approaches from {date_app_min} to {date_app_max}")
            if asteroid.loc[n_ast,'is_NHATS']==True: additional_info.append('Included in NHATS database')
            if asteroid.loc[n_ast,'is_geometry']==True: additional_info.append('Geometry model available' )
            if asteroid.loc[n_ast,'PHA']=='Y': additional_info.append('PHA asteroid' )

            # Additional info (interest) is 'a priori' requirement (not shown in decision table but computed)
            if additional_info==[]: continue 
            asteroid.loc[n_ast, 'additional_info'] = ', '.join(additional_info)
            summary_table.add_row([asteroid.loc[n_ast,'ID'], round(asteroid.loc[n_ast,'delta_v_tot'],2), 
                                  f"{int(asteroid.loc[n_ast, 'n_backup'])} ({familiar_value}, delta_H={round(asteroid.loc[n_ast, 'delta_H'], 2)})",
                                   asteroid.loc[n_ast,'condition_code'], 
                                   round(asteroid.loc[n_ast,'period_sin'],2),asteroid.loc[n_ast,'Spin period'],
                                   asteroid.loc[n_ast,'additional_info']])
        
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
    a_candidate=np.array([])
    e_candidate=np.array([])
    i_candidate=np.array([])
    u_inf_candidate=np.array([])
    delta_v_candidate=np.array([])

    for n_ast in range(0,len(asteroid)):
        a[n_ast]=asteroid.loc[n_ast,'a'] #a
        e[n_ast]=asteroid.loc[n_ast,'e'] #e
        i[n_ast]=asteroid.loc[n_ast,'i'] #i
        u_inf[n_ast]=np.sqrt(abs(3-(1/a[n_ast]+2*np.sqrt(a[n_ast]*(1-e[n_ast]**2))*np.cos(i[n_ast]*np.pi/180)))) 
        delta_v[n_ast]=asteroid.loc[n_ast,'delta_v_tot'] #delta_v_tot
        if asteroid.loc[n_ast, 'n_backup'] >= 1 and asteroid.loc[n_ast, 'Spin period'] is not None and \
            pd.notna(asteroid.loc[n_ast, 'additional_info']):
            a_candidate = np.append(a_candidate,  a[n_ast])  # Añadir el valor
            e_candidate = np.append(e_candidate,  e[n_ast])
            i_candidate = np.append(i_candidate,  i[n_ast])
            u_inf_candidate = np.append(u_inf_candidate, u_inf[n_ast])
            delta_v_candidate = np.append(delta_v_candidate, delta_v[n_ast])
           
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
    ax1.plot(a,e,marker='o',markeredgecolor='k',markersize=6,color='c',linestyle='',label='NEAs w. close approach')
    ax1.plot(a_candidate,e_candidate,marker='o',markeredgecolor='k',markersize=6,color='r',linestyle='',label='NEAs candidates')
    ax1.set_xlabel("a [au]")
    ax1.set_ylabel("e ")
    ax1.set_xlim(0,4)
    ax1.set_ylim(0,1)
    ax1.grid(True)
    ax1.legend()

    ax2=plt.subplot(2,2,2) # plot a vs i
    ax2.plot(a_others,i_others,marker='+',markersize=4,color='b',linestyle='',label='NEAs analized')
    ax2.plot(a,i,marker='o',markeredgecolor='k',markersize=6 ,color='c',linestyle='',label='NEAs w. close approach')
    ax2.plot(a_candidate,i_candidate, marker='o',markeredgecolor='k',markersize=6 ,color='r',linestyle='',label='NEAs candidates')
    ax2.set_xlabel("a [au]")
    ax2.set_ylabel("i [deg]")
    ax2.set_xlim(0,4)
    ax2.set_ylim(0,90)
    ax2.grid(True)
    ax2.legend()

    ax3=plt.subplot(2,2,3) # plot u_inf vs delta_v
    ax3.plot(u_inf,delta_v,marker='o',markeredgecolor='k',color='c',markersize=8 ,linestyle='',label='NEAs w. close approach')
    ax3.plot(u_inf_candidate,delta_v_candidate,marker='o',markeredgecolor='k',color='r',markersize=8 ,linestyle='',label='NEAs candidates')
    ax3.set_xlabel("U_inf")
    ax3.set_ylabel("Delta_v [m/s]")
    ax3.set_xlim(0,1.4)
    ax3.set_ylim(bottom=0)
    ax3.grid(True)
    ax3.legend()

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