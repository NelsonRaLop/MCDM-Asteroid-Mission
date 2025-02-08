def U_fuzzy(asteroid_candidates):
    '''
    This funcion computes the orbit uncertainty as a fuzzy number based on the
    ranges of the definition of U - Condition Code
    See http://www.minorplanetcenter.org/iau/info/UValue.html

    return none
    '''
    for n_ast in range(len(asteroid_candidates)):
        U_asteroid=int(asteroid_candidates.loc[n_ast,'condition_code'])
        # Note that asteroids with condition_code=9 have been removed
        if U_asteroid==0: 
            inf=0 
            sup=1 
        elif U_asteroid==1: 
            inf=1
            sup=4.4 
        elif U_asteroid==2: 
            inf=4.4 
            sup=19.6 
        elif U_asteroid==3: 
            inf=19.6 
            sup=86.5 
        elif U_asteroid==4: 
            inf=86.5 
            sup=382 
        elif U_asteroid==5: 
            inf=382 
            sup=1692 
        elif U_asteroid==6: 
            inf=1692 
            sup=7488 
        elif U_asteroid==7: 
            inf=7488 
            sup=33121 
        elif U_asteroid==8: 
            inf=33121 
            sup=146502   

        U_fuzzy=[inf,(sup-inf)/2,sup]
        asteroid_candidates.loc[n_ast,'U_fuzzy']=str(U_fuzzy)
    return

def weights_merec(decision_table,criterion_type):
    '''
    This funcion compute the weight using Merec Method 

    return merec weights
    '''
    from pyDecision.algorithm import merec_method
    
    weights=merec_method(decision_table,criterion_type)

    return weights