import pandas as pd
import numpy as np

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

        U_fuzzy=(inf,(sup-inf)/2,sup)
        asteroid_candidates.loc[n_ast,'U_fuzzy']=str(U_fuzzy)
    return

def decision_matrix(asteroid,criteria):
    '''
    This functions generates the decision matrix given the complete information
    of the asteroid analysed and the parameters selected as criteria

    return decision_matrix
    
    '''

    decision_matrix=pd.DataFrame(columns=criteria)
    decision_matrix=asteroid[criteria].copy()
    decision_matrix.loc[:, 'U_fuzzy'] = decision_matrix['U_fuzzy'].apply(lambda x: eval(x) if isinstance(x, str) else x)
    print (decision_matrix)
    return decision_matrix



def weights_fuzzy_merec(decision_matrix,criteria_type):
    '''
    This function computes the weights using Fuzzy Merec Method

    1.- Transform crisp values to triangular number
    2.- Implement Fuzzy Merec Weighting method as described in:

    https://www.journal-fea.com/article_206242_cecc86384c1b45f480d3af55d9507fcf.pdf

    To do so first transform crisp values to triangular number

    return fuzzy_merec_weights
    '''
    # 1.- Crisp values to triangular numbers a -> (a,a,a)
    decision_matrix_transformed = []
    
    for _, row in decision_matrix.iterrows():
        transformed_row = [
            (row['delta_v_tot'], row['delta_v_tot'], row['delta_v_tot']),  # Crisp to fuzzy
            (row['n_backup'], row['n_backup'], row['n_backup']),            # Crisp to fuzzy
            row['U_fuzzy'],                                                # Fuzzy
            (row['period_sin'], row['period_sin'], row['period_sin']),      # Crisp to fuzzy
            (row['Spin period'], row['Spin period'], row['Spin period'])    # Crisp to fuzzy
        ]
        # Create matrix with the correct format
        decision_matrix_transformed.append(transformed_row)
        decision_matrix_format=list(decision_matrix_transformed)

    m = len(decision_matrix_format[0]) # Number of criteria
    n = len(decision_matrix_format) # Number of alternatives
    
    # 2.- Fuzzy Merec method implementation
    fuzzy_merec_weights = np.zeros(m)
    X_a = np.zeros((n,m))
    X_a_norm = np.zeros((n,m))
    X_b = np.zeros((n,m))
    X_b_norm = np.zeros((n,m))
    X_c = np.zeros((n,m))
    X_c_norm = np.zeros((n,m))
    Q_dot_a  = np.zeros((n,m))
    Q_dot_b  = np.zeros((n,m))
    Q_dot_c = np.zeros((n,m))
    Q_a = np.zeros(n)
    Q_b = np.zeros(n)
    Q_c = np.zeros(n)
    V_dist = np.zeros((n,m))

    # Get triangular numbers
    for j in range(0, m):
        for i in range(0, n):
            a, b, c  = decision_matrix_format[i][j]
            X_a[i,j] = a
            X_b[i,j] = b
            X_c[i,j] = c

    # Normalization        
    for j in range(0, m):
        if (criteria_type[j] == 'max'):
            X_a_norm[:, j]  = np.min(X_a[:, j]) / X_c[:, j]
            X_b_norm[:, j]  = np.min(X_a[:, j]) / X_b[:, j]
            X_c_norm[:, j]  = np.min(X_a[:, j]) / X_a[:, j]
        else:
            X_a_norm[:, j]  = X_a[:, j]/np.max(X_c[:, j]) 
            X_b_norm[:, j]  = X_b[:, j]/np.max(X_c[:, j]) 
            X_c_norm[:, j]  = X_c[:, j]/np.max(X_c[:, j]) 


    # Overall performance Q and performance with criteria removed Q_dot
    for i in range(0, n):
        Q_a[i] = np.log(1+(1/m*np.sum(abs(np.log(X_c_norm[i,:])))))
        Q_b[i] = np.log(1+(1/m*np.sum(abs(np.log(X_b_norm[i,:])))))
        Q_c[i] = np.log(1+(1/m*np.sum(abs(np.log(X_a_norm[i,:])))))
        for j in range(0,m):
            #Remove index of the corresponding criteria    
            X_a_norm_rem = np.delete(X_a_norm[i, :], j)
            X_b_norm_rem = np.delete(X_b_norm[i, :], j)
            X_c_norm_rem = np.delete(X_c_norm[i, :], j)
            # Compute performance without the criteria
            Q_dot_a[i,j] = np.log(1+(1/m*np.sum(abs(np.log(X_c_norm_rem)))))
            Q_dot_b[i,j] = np.log(1+(1/m*np.sum(abs(np.log(X_b_norm_rem)))))
            Q_dot_c[i,j] = np.log(1+(1/m*np.sum(abs(np.log(X_a_norm_rem)))))

    # Euclidean deviation
    for i in range(0,n):
        for j in range(0,m):
            V_dist[i,j]=np.sqrt((Q_c[i]-Q_dot_c[i,j])**2+(Q_b[i]-Q_dot_b[i,j])**2+(Q_a[i]-Q_dot_a[i,j])**2)

    V_dist_total = np.sum(V_dist[:,:])

    # Final weights

    for j in range (0,m):
        fuzzy_merec_weights[j] = np.sum(V_dist[:,j])/V_dist_total

    return fuzzy_merec_weights