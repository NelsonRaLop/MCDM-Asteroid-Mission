### Import packages ###
import requests
import re
import pandas as pd
import datetime

def asteroid_load(H_min,H_max,date_app_min,date_app_max,dist_app):
    """
    This function analyses diferent databases to obtain information of the asteroids analysed

    return df_asteroid,asteroid_app,asteroid_NHATS,asteroid_geometry,archive_in

    1.- JPL - SmallBody Data Base Query: To obtain orbital and physical parameters of the asteroids  
    2.- JPL CNEOS - NEO Earth Close Approaches: To analyse earth close approches under the period of analysis
    3.- JPL CNEOS - Near-Earth Object Human Space Flight Accessible Targets Study (NHATS) 
    4.- DAMIT - Database of Asteroid Models from Inversion Techniques: To check is geometry model os available

    APIS will be applied when available

    """

    ### 1.-
    url=f'https://ssd-api.jpl.nasa.gov/sbdb_query.api?fields=full_name,H,e,a,q,i,om,w,ad,ma,per_y,moid,condition_code,spec_B,rot_per,sats,pha&sb-kind=a&sb-group=neo&full-prec=true&sb-cdata=%7b%22AND%22%3a%5b%22H%7cRG%7c{H_min}%7c{H_max}%22%5d%7d'
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
    df_asteroid = pd.DataFrame(asteroid, columns=[
        "ID", "H", "e", "a", "Perigee", "i", "RAAN", "arg_perig", "Apogee", "M", "period [y]",
        "MOID", "condition_code", "SMASS taxonomy", "Spin period", "Satellites", "PHA"
    ])
    archive_in = f'asteroid_input_{H_min}_{H_max}.csv'
    df_asteroid.to_csv(archive_in, index=False)


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

    return df_asteroid,asteroid_app,asteroid_NHATS,asteroid_geometry,archive_in




def asteroid_filtering(n_ast,asteroid, asteroid_removed,asteroid_app,asteroid_NHATS, asteroid_geometry):
    """
    This function checks for close approaches and any other additional data

    If the asteroid wont have any close encounter in the period defined, simply remove it as candidate

    return cont,approaches,approach_date,is_NHATS,is_geometry,asteroid_removed

    1.- Check close approaches & check if geometry model is available and Check if the asteroid is included in NHATS list
    2.- Compute the number of close approaches in the period considered and get the date of first Earth approache as reference
    3.- Remove if NO approches

    """
 # 1.-
    match = re.match(r'^(\d+)', asteroid.loc[n_ast,'ID'].strip()) # Check if numbered
    if match:
        name = match.group(1)  # Extract number (ex, "756476")

        # Only few *NUMBERED* asteroid have geometry models available
        if name in asteroid_geometry: is_geometry=True
        else: is_geometry=False

    else:
        # Is not numbered, provisional id "(1993 BU3)")
        name = asteroid.loc[n_ast, 'ID'].strip().strip("()")
        is_geometry=False
    
    if name in asteroid_NHATS: is_NHATS=True
    else: is_NHATS=False
    
    # 2.- 
    approaches = 0 # Initialization
    approach_date_raw = approach_date = None  #  Initialization

    for approach in asteroid_app:
        if approach[0] == name:
            approaches += 1  # Increase the number of approaches in the period analysed
            if approach_date_raw is None:  # Save only the first close approach
                approach_date_raw = approach[3].split(" ")[0]  # Get the date of the first approach
                approach_date_time_format = datetime.datetime.strptime(approach_date_raw, "%Y-%b-%d") # Transform in time object
                approach_date = approach_date_time_format.strftime("%Y-%m-%d") # Correct format e.g. "2023-01-27" 
    
    # 3.-
    if approaches == 0 or asteroid.loc[n_ast,'condition_code']==9:
        asteroid_removed = pd.concat([asteroid_removed, asteroid.loc[[n_ast]]], ignore_index=True)
        asteroid.drop(n_ast, inplace=True)
        cont=0   
    else: cont=1

    return cont,approaches,approach_date,is_NHATS,is_geometry,asteroid_removed