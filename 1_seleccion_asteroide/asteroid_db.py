### Import packages ###
import requests
import csv
import re

def asteroid_load(H_min,H_max,date_app_min,date_app_max,dist_app):
    """
    This function analyses diferent databases to obtain information of the asteroids analysed

    return (asteroid, asteroid_all,asteroid_app,asteroid_NHATS,asteroid_geometry,archive_in)

    1.- JPL - SmallBody Data Base Query: To obtain orbital and physical parameters of the asteroids 
    2.- JPL CNEOS - NEO Earth Close Approaches: To analyse earth close approches under the period of analysis
    3.- JPL CNEOS - Near-Earth Object Human Space Flight Accessible Targets Study (NHATS) 
    4.- DAMIT - Database of Asteroid Models from Inversion Techniques: To check is geometry model os available

    APIS will be applied when available

    """

    ### 1.-
    url=f'https://ssd-api.jpl.nasa.gov/sbdb_query.api?fields=full_name,H,e,a,q,i,om,w,ad,per_y,moid,condition_code,spec_B,rot_per,sats,pha&sb-kind=a&sb-group=neo&full-prec=true&sb-cdata=%7b%22AND%22%3a%5b%22H%7cRG%7c{H_min}%7c{H_max}%22%5d%7d'
    response=requests.get(url)
    response_json=response.json()
    asteroid = [
        [
            str(field).strip() if isinstance(field, str) else field
            for field in asteroid_data
        ]
        for asteroid_data in response_json.get("data", [])
    ]
    
    #Write a file with raw data
    archive_in=f'asteroid_input_{H_min}_{H_max}.csv'
    with open(archive_in, mode='w',newline='') as csvfile: 
        writer=csv.writer(csvfile)
        writer.writerow(["0 ID, 1 H, 2 e, 3 a, 4 Perigee, 5 i, 6 RAAN, 7 arg_perig, 8 Apogee, 9 period [y],10 MOID, 11 condition_code, 12 SMASS taconomy, 13 Spin period,14 Satellites, 15 PHA"])
        writer.writerows(asteroid)


    #Check all asteroid DB to compute Orbital Afinity and load data in asteroid_db
    url='https://ssd-api.jpl.nasa.gov/sbdb_query.api?fields=full_name,H,e,a,q,i,om,w,ad,per_y,moid,condition_code&sb-kind=a&sb-group=neo&full-prec=true'
    response=requests.get(url)
    response_json=response.json()
    asteroid_all = [
        [field.strip() if field is not None else '' for field in db_asteroid_data]
        for db_asteroid_data in response_json.get("data", [])
    ]


    ### 2.-
    url=f'https://ssd-api.jpl.nasa.gov/cad.api?dist-max={dist_app}&date-min={date_app_min}&date-max={date_app_max}&h-min={H_min}&h-max={H_max}&nea=true&sort=object'
    response=requests.get(url)
    response_json=response.json()
    asteroid_app = [
        [field.strip() if field is not None else '' for field in approaches_data]
        for approaches_data in response_json.get("data", [])
    ]


    ### 3.-
    url='https://ssd-api.jpl.nasa.gov/nhats.api?launch=2035-2040&h=30'
    response=requests.get(url)
    response_json=response.json()
    asteroid_NHATS = str(response.json())  # Dictionary to strin

    ### 4.-
    url='https://astro.troja.mff.cuni.cz/projects/damit/exports/table/asteroids'
    response=requests.get(url)
    asteroid_geometry = response.text

    return asteroid,asteroid_all,asteroid_app,asteroid_NHATS,asteroid_geometry,archive_in




def asteroid_filtering(filter,n_ast,asteroid,asteroid_removed,asteroid_app,asteroid_NHATS, asteroid_geometry):
    """
    This function checks for close approaches and other additional data

    If the asteroid wont have any close encounter in the period defined, simply remove it as candidate
    If filter=1 only PHA with spin OR taxonomy OR with secondary body are considered candidates

    return (cont,approaches,is_NHATS,is_geometry,asteroid_removed)

    1.- Check close approaches & check if geometry model is available 
    2.- Compute the number of close approaches in the period considered
    3.- Remove if NO approches
    4.- Apply additional filtering if filter==1 
    5.- Check if the asteroid is included in NHATS list

    """

 # 1.-
    match = re.match(r'^(\d+)', asteroid[n_ast][0].strip()) # Check if numbered
    if match:
        name = match.group(1)  # Extract number (ex, "756476")

        # Only few *NUMBERED* asteroid have geometry models available
        if name in asteroid_geometry: is_geometry=True
        else: is_geometry=False

    else:
        # Is not numbered, provisional id "(1993 BU3)")
        name = asteroid[n_ast][0].strip().strip("()")
        is_geometry=False

    # 2.- 
    approaches=0
    approaches = sum(1 for approach in asteroid_app if approach[0] == name)
    
    # 3.-
    if approaches == 0:
        asteroid_removed.append(asteroid[n_ast])
        del asteroid[n_ast]
        cont=0

    # 4.-  
    elif filter==1: 
      if asteroid[n_ast][12]==None and asteroid[n_ast][13]==None and \
      asteroid[n_ast][14]==0 or asteroid[n_ast][15]=='N':
        asteroid_removed.append(asteroid[n_ast])
        del asteroid[n_ast]
        cont=0     
    else: cont=1

    # 5.-
    if name in asteroid_NHATS: is_NHATS=True
    else: is_NHATS=False

    return cont,approaches,is_NHATS,is_geometry