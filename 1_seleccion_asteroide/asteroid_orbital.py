import numpy as np

### Load constant values  ###
NU_SUN=1.32712440018e20 #Sun gravitational constant [SI]
NU_EARTH=3.986004418e14 #Earth gravitational constant [SI]
PERIOD_EARTH=365.25 #Earth period [Days]
AU=149597870700 # 1 Astronomic unit to meters

def asteroid_similarity(n_ast,asteroid,asteroid_all):
    '''
    This function checks the orbit similarity between the asteroid and other
    asteroids of the data base
    
    return g_2,familiar

    1.- Extract orbital data
    2.- For each asteroid in the DB compute similarity via g_2 parameter
    3.- Get the asteroid with an orbit most similar and record g_2
    
    '''

    # 1.- 
    ex=float(asteroid[n_ast][2])
    a=float(asteroid[n_ast][3]) #[ua]
    inc=float(asteroid[n_ast][5])*2*np.pi/360 #[rad]
    arg_per=float(asteroid[n_ast][7])*2*np.pi/360 #[rad]
    asc_node=float(asteroid[n_ast][6])*2*np.pi/360 #[rad]
    p=a*(1-ex**2) #[ua]

    # 2.-
    g_2=100 #Initialization
    familiar='none' #Initialization
    for n_ast_db in range(0,len(asteroid_all)):

        if asteroid[n_ast][0] == asteroid_all[n_ast_db][0]: continue #Skip itself

        # Orbit parameters of all asteroid in DB
        ex_2=float(asteroid_all[n_ast_db][2])
        a_2=float(asteroid_all[n_ast_db][3])
        inc_2=float(asteroid_all[n_ast_db][5])*2*np.pi/360
        arg_per_2=float(asteroid_all[n_ast_db][7])*2*np.pi/360
        asc_node_2=float(asteroid_all[n_ast_db][6])*2*np.pi/360
        p_2=a_2*(1-ex_2**2)

        # Compute afinity/similarity parameter
        cos_I=np.cos(inc)*np.cos(inc_2)+np.sin(inc)*np.sin(inc_2)*np.cos(asc_node-asc_node_2)
        cos_P=np.sin(inc)*np.sin(inc_2)*np.sin(arg_per)*np.sin(arg_per_2)+ \
                np.cos(arg_per)*np.cos(arg_per_2)*np.cos(asc_node-asc_node_2)+ \
                np.cos(inc)*np.cos(inc_2)*np.sin(arg_per)*np.sin(arg_per_2)*np.cos(asc_node-asc_node_2)+ \
                (np.cos(inc_2)*np.cos(arg_per)*np.sin(arg_per_2)-\
                np.cos(inc)*np.sin(arg_per)*np.cos(arg_per_2))*np.sin(asc_node-asc_node_2)  

    # 3.-
        g_2_ast=(1+ex**2)*p + (1+ex_2**2)*p_2 - 2*np.sqrt(p*p_2)*(cos_I+ex*ex_2*cos_P)
        if g_2_ast<g_2:
            g_2=g_2_ast #Update metric
            familiar=asteroid_all[n_ast_db][0] #Update the most similar asteroid
        

    return g_2,familiar


def asteroid_period(n_ast,asteroid):
    '''
    This function computes the synodic period of the asteroid [years]
    '''
    period=float(asteroid[n_ast][9]) #[y]
    period_sin=1/abs(1/period-1) #[y]

    return period_sin


def asteroid_accessibility(n_ast, asteroid):
    '''
    This function computes an approximation of the delta_v to capture the asteroid

    return delta_v_tot [m/s]

    1.- Extract orbital data
    2.- For crossing-Earth-orbit asteroids computes delta_v=min(change arg_perigee 
        or i to obtain an intersection 
    3.- For outer asteroids computes delta_v to increment a
    4.- For inner asteroids computes delta_v to reduce a
    5.- Estimates delta_v to achieve the parabolic capture limit
    6.- Estimates delta_v_total of al the maneuvre
    '''

    # 1.-
    ex=float(asteroid[n_ast][2])
    a=float(asteroid[n_ast][3]) #[au]
    a_m=a*AU #[m]
    inc=float(asteroid[n_ast][5])*2*np.pi/360 #[rad]
    arg_per=float(asteroid[n_ast][7])*2*np.pi/360 #[rad]
    p=a*(1-ex**2) #[au]
    p_m=a_m*(1-ex**2) #[m]
    q=float(asteroid[n_ast][4])#[au]
    Q=float(asteroid[n_ast][8])#[au]

    # 2.- 
    if q<=1 and Q>=1: #'crossing'
        # Change arg_perigee
        theta_enc=np.arccos((p-1)/ex) 
        w_enc=(np.pi-theta_enc,theta_enc+np.pi,theta_enc,2*np.pi-theta_enc)
        w_dif=(abs(arg_per-w_enc[0]),abs(arg_per-w_enc[1]),abs(arg_per-w_enc[2]),abs(arg_per-w_enc[3]))
        delta_dif=min(w_dif)
        delta_v_arg=2*np.sqrt(NU_SUN/p_m)*ex*np.sin(delta_dif/2)
        
        # Change inclination
        r_node_m=(p_m/(1+ex*np.cos(arg_per)),p_m/(1+ex*np.cos(arg_per-np.pi))) #Ascending and descending node
        r_node_max_m=max(r_node_m) #The node with the smallest speed
        v_nodes=np.sqrt(NU_SUN*(2/r_node_max_m-1/a_m))
        delta_v_inc=2*v_nodes*np.sin(inc/2)
        if delta_v_inc<delta_v_arg: inc=0 #Updates inclination if change i is cheaper

        # Crossing delta_v
        delta_v_cross=min(delta_v_inc,delta_v_arg)

    # 3.-
    elif q>=1: #'outer'
      r_apo=Q
      r_apo_m=r_apo*AU
      if np.pi/2 <= arg_per < 3*np.pi/2: #Closest node to the perigee
        theta_close = arg_per - np.pi
      else:
        theta_close = arg_per
      ex_tras=(Q-1)/(np.cos(theta_close)+Q)
      a_tras_m=a_m*(1+ex)/(1+ex_tras)
      delta_v_cross=np.sqrt(NU_SUN)*abs(np.sqrt(2/r_apo_m-1/a_tras_m)-np.sqrt(2/r_apo_m-1/a_m)) #Hohman
      ex=ex_tras #Update orbital parameters after the maneuvre
      a=a_tras_m/AU #Update orbital parameters after the maneuvre

    # 4.-
    elif Q<=1: #'inner'
      r_peri=q
      r_peri_m=r_peri*AU 
      if np.pi/2 <= arg_per < 3*np.pi/2: #Closest node to the apogee
        theta_close = arg_per
      else:
        theta_close = arg_per-np.pi
      ex_tras=(q-1)/(np.cos(theta_close)-q)
      a_tras_m=a_m*(1-ex)/(1-ex_tras)
      delta_v_cross=np.sqrt(NU_SUN)*abs(np.sqrt(2/r_peri_m-1/a_tras_m)-np.sqrt(2/r_peri_m-1/a_m)) #Hohman
      ex=ex_tras #Update orbital parameters after the maneuvre
      a=a_tras_m/AU #Update orbital parameters after the maneuvre

    # 5.- Parabolic capture with perigee of r_p [m]
    r_p=200000 #[m] based on LR
    tisse=1/a+2*np.sqrt(a*(1-ex**2))*np.cos(inc)
    U=np.sqrt(3-tisse) #Non-dymensional Hyperbolic speed 
    v_inf=U*np.sqrt(NU_SUN/AU)
    delta_v_cap=np.sqrt(2*NU_EARTH/r_p+v_inf**2)-np.sqrt(2*NU_EARTH/r_p)

    # 6.- Total budget = crossing + capture
    delta_v_tot=delta_v_cap+delta_v_cross 

    return delta_v_tot