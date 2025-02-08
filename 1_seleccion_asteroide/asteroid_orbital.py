import numpy as np
import astropy.units as u
from astropy.time import Time
from poliastro.bodies import Sun
from poliastro.twobody import Orbit

### Load constant values  ###
NU_SUN=1.32712440018e20 #Sun gravitational constant [SI]
NU_EARTH=3.986004418e14 #Earth gravitational constant [SI]
PERIOD_EARTH=365.25 #Earth period [Days]
AU=149597870700 # 1 Astronomic unit to meters

def asteroid_similarity(n_ast,asteroid,asteroid_all,g_2_threshold,rel_distance_threshold):
    '''
    This function calculates the number of possible back up asteroids once the mission is in flight 
    in case the spacecraft has to be redirected to another asteroid. This check will be performed exclusively 
    among all asteroids that have a natural approach to the Earth and that have not been eliminated by previous 
    filters. To do this:

    1.- Extract orbital data and calculate the astronomical position of the primary target at the reference time 
    (by default, the first approach to Earth).
    
    2.- Iterate over the remaining candidate asteroids by calculating the number of potential backup asteroids. For 
    an asteroid to be considered a backup asteroid, it must (1) have a sufficiently similar orbit - checked through 
    the limit of the orbital similarity metric (2) in the mean anomaly state corresponding to the reference time, 
    the geometric distance between primary asteroid and backup asteroid is small - checked through the limit of the distance.

    These requirements are a first approximation of how many asteroids the mission could deviate to with presumably small v deltas.

    return n_backup
    
    '''

    # 1.- 
    n_backup=0
    distance_threshold_m=rel_distance_threshold*2*np.pi*AU # To express threshold in meters
    approach_date=asteroid.loc[n_ast,'date_approach']
    target_date=Time(approach_date, scale="tdb")

    ex=float(asteroid.loc[n_ast,'e'])
    a=float(asteroid.loc[n_ast,'a']) #[ua]
    inc=float(asteroid.loc[n_ast,'i'])*2*np.pi/360 #[rad]
    arg_per=float(asteroid.loc[n_ast,'arg_perig'])*2*np.pi/360 #[rad]
    asc_node=float(asteroid.loc[n_ast,'RAAN'])*2*np.pi/360 #[rad]
    M_ast=float(asteroid.loc[n_ast,'M'])
    H_size=float(asteroid.loc[n_ast,'H'])
    p=a*(1-ex**2) #[ua]
    primary_orbit = Orbit.from_classical(Sun, a*u.AU, ex * u.dimensionless_unscaled, inc* u.deg,
             asc_node* u.deg, arg_per* u.deg, M_ast* u.deg, epoch=target_date)          
    pos1 = primary_orbit.propagate(target_date - primary_orbit.epoch).r
    
    # 2.-
    for n_ast_db in range(0,len(asteroid_all)):

        if asteroid.loc[n_ast,'ID'] == asteroid_all.loc[n_ast_db,'ID']: continue #Skip itself
        # Orbit parameters of all asteroid in DB
        ex_2=float(asteroid_all.loc[n_ast_db,'e'])
        a_2=float(asteroid_all.loc[n_ast_db,'a'])
        inc_2=float(asteroid_all.loc[n_ast_db,'i'])*2*np.pi/360
        arg_per_2=float(asteroid_all.loc[n_ast_db,'arg_perig'])*2*np.pi/360
        asc_node_2=float(asteroid_all.loc[n_ast_db,'RAAN'])*2*np.pi/360
        M_ast_2=float(asteroid_all.loc[n_ast_db,'M'])
        H_size_2=float(asteroid_all.loc[n_ast_db,'H'])
        p_2=a_2*(1-ex_2**2)

        # Compute afinity/similarity parameter
        cos_I=np.cos(inc)*np.cos(inc_2)+np.sin(inc)*np.sin(inc_2)*np.cos(asc_node-asc_node_2)
        cos_P=np.sin(inc)*np.sin(inc_2)*np.sin(arg_per)*np.sin(arg_per_2)+ \
                np.cos(arg_per)*np.cos(arg_per_2)*np.cos(asc_node-asc_node_2)+ \
                np.cos(inc)*np.cos(inc_2)*np.sin(arg_per)*np.sin(arg_per_2)*np.cos(asc_node-asc_node_2)+ \
                (np.cos(inc_2)*np.cos(arg_per)*np.sin(arg_per_2)-\
                np.cos(inc)*np.sin(arg_per)*np.cos(arg_per_2))*np.sin(asc_node-asc_node_2)  

        g_2_ast=np.sqrt((1+ex**2)*p + (1+ex_2**2)*p_2 - 2*np.sqrt(p*p_2)*(cos_I+ex*ex_2*cos_P))
        
        if g_2_ast <g_2_threshold:
            secondary_orbit = Orbit.from_classical(Sun, a_2*u.AU, ex_2 * u.dimensionless_unscaled, inc_2* u.deg,
             asc_node_2* u.deg, arg_per_2* u.deg, M_ast_2* u.deg, epoch=target_date) 
            pos2 = secondary_orbit.propagate(target_date - secondary_orbit.epoch).r
            distance = np.linalg.norm((pos1 - pos2).to(u.m).value)  # Distance [km]s
            if distance<distance_threshold_m: 
              n_backup=n_backup+1
            # # Save from the list of candidates the most similar one in size, 'familiar'
            #   delta_H_candidate=abs(H_size_2-H_size)
            #   if delta_H_final == []: # Initialization
            #      delta_H_final = delta_H_candidate
            #      familiar=asteroid_all.loc[n_ast_db,'ID']
            #   else:
            #      if delta_H_candidate<delta_H_final:
            #         delta_H_final=delta_H_candidate
            #         familiar=asteroid_all.loc[n_ast_db,'ID']
            #      else:continue
            
        else:continue
        
    return n_backup #, familiar, delta_H_final


def asteroid_period(n_ast,asteroid):
    '''
    This function computes the synodic period of the asteroid [years]
    '''
    period=float(asteroid.loc[n_ast,'period [y]']) #[y]
    period_sin=1/abs(1/period-1) #[y]

    return period_sin


def asteroid_accessibility(n_ast, asteroid):
    '''
    This function computes an approximation of the delta_v to capture the asteroid as described in 
    https://www.researchgate.net/publication236163825_Near-Earth_asteroid_resource_accessibility_and_future_capture_mission_opportunities

    
    return delta_v_tot [m/s]

    1.- Extract orbital data
    2.- For crossing-Earth-orbit asteroids computes delta_v=min(change arg_perigee 
        or i to obtain an intersection with Earth orbit
    3.- For outer asteroids computes delta_v to increment a to get an intersection with Earth orbit
    4.- For inner asteroids computes delta_v to reduce a to get an intersection with Earth orbit
    5.- Estimates delta_v to achieve the parabolic capture limit
    6.- Estimates delta_v_total of the complete maneuvre
    '''

    # 1.-
    ex=float(asteroid.loc[n_ast,'e'])
    a=float(asteroid.loc[n_ast,'a']) #[au]
    a_m=a*AU #[m]
    inc=float(asteroid.loc[n_ast,'i'])*2*np.pi/360 #[rad]
    arg_per=float(asteroid.loc[n_ast,'arg_perig'])*2*np.pi/360 #[rad]
    p=a*(1-ex**2) #[au]
    p_m=a_m*(1-ex**2) #[m]
    q=float(asteroid.loc[n_ast,'Perigee'])#[au]
    Q=float(asteroid.loc[n_ast,'Apogee'])#[au]

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