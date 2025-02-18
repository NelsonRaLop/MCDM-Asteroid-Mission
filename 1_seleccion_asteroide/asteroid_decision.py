import pandas as pd
import numpy as np
from scipy.optimize import minimize

#################################
### Functions for fuzzy MEREC ###
#################################


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

    return decision_matrix and decision_matrix_TFN (crisp numbers into TFN)
    
    '''

    decision_matrix=pd.DataFrame(columns=criteria)
    decision_matrix=asteroid[criteria].copy()
    decision_matrix.loc[:, 'U_fuzzy'] = decision_matrix['U_fuzzy'].apply(lambda x: eval(x) if isinstance(x, str) else x)
    
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
        decision_matrix_TFN=list(decision_matrix_transformed)
 
    return decision_matrix, decision_matrix_TFN


def weights_fuzzy_merec(decision_matrix_TFN,criteria_type):
    '''
    This function computes the weights using Fuzzy Merec Method

    Implement Fuzzy Merec Weighting method as described in:

    https://www.journal-fea.com/article_206242_cecc86384c1b45f480d3af55d9507fcf.pdf

    Note that to do so first transformation of crisp values into triangular number is required

    return fuzzy_merec_weights
    '''
    # 1.- Initialization

    m = len(decision_matrix_TFN[0]) # Number of criteria
    n = len(decision_matrix_TFN) # Number of alternatives
    
    # 2.- Fuzzy Merec method implementation
    weights = np.zeros(m)
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
            a, b, c  = decision_matrix_TFN[i][j]
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
        weights[j] = np.sum(V_dist[:,j])/V_dist_total

    return weights





#################################
### Functions for fuzzy FUCOM ###
#################################


# Objective function: to monimice chi
def objective(w, fi, phi_max):
    
    '''
    This function defines the objective functions to optimice (minimice chi) when aplying fuzzy fucom
    as depicted in eq 13 of:

    https://pdfs.semanticscholar.org/809b/d4cc55c66036264b463234d9e582471e9d31.pdf


    return none
    
    '''
        
    a = len(fi) // 3  # Number of criteria
    w_matrix = w.reshape((a, 3))
    chi = 0
    
    # Optimization of the value of the elements that are consecutives
    for j in range(a - 1):
        chi = max(chi, abs(w_matrix[j, 0] - w_matrix[j+1, 2] * fi[3*j + 3]))
        chi = max(chi, abs(w_matrix[j, 2] - w_matrix[j+1, 0] * fi[3*j + 5]))
        chi = max(chi, abs(w_matrix[j, 1] - w_matrix[j+1, 1] * fi[3*j + 4]))
    
    # Optimization of the value of the elements that are NO consecutives
    for j in range(a - 2):
        chi = max(chi, abs(w_matrix[j, 0] - w_matrix[j+2, 2] * fi[3*j + 3] * fi[3*j + 6]))
        chi = max(chi, abs(w_matrix[j, 2] - w_matrix[j+2, 0] * fi[3*j + 5] * fi[3*j + 8]))
        chi = max(chi, abs(w_matrix[j, 1] - w_matrix[j+2, 1] * fi[3*j + 4] * fi[3*j + 7]))
    
    return chi

# Constraints

    # Here the constraint of the optimization problem to be solved when aplying fuzzy fucom
    # are defined as depicted in eq 13 of:

    # https://pdfs.semanticscholar.org/809b/d4cc55c66036264b463234d9e582471e9d31.pdf


    # return none


def constraint_h(w, a):
    w_matrix = w.reshape((a, 3))
    h = np.sum((w_matrix[:, 0] + 4 * w_matrix[:, 1] + w_matrix[:, 2]) / 6)
    return h - 1  # sum of h[j] MUST be 1

def constraint_w0(w, a):
    w_matrix = w.reshape((a, 3))
    return w_matrix[:, 0]  # w[j,0] >= 0

def constraint_w1(w, a):
    w_matrix = w.reshape((a, 3))
    return w_matrix[:, 1] - w_matrix[:, 0]  # w[j,1] >= w[j,0]

def constraint_w2(w, a):
    w_matrix = w.reshape((a, 3))
    return w_matrix[:, 2] - w_matrix[:, 1]  # w[j,2] >= w[j,1]


# Funtion to compute fi vector from comparison vector
def compute_fi(comparison):
    '''
    This function computes the fi vector (comparation between weights) from the comparison
    vector obtained from linguistic comparison. See:

    https://pdfs.semanticscholar.org/809b/d4cc55c66036264b463234d9e582471e9d31.pdf


    return fi
    
    '''
    # Initialization
    fi = comparison[:6].copy()
    
    # Compute fi from comparison vector
    for i in range(6, len(comparison), 3):
        trio_actual = comparison[i:i+3]
        trio_anterior = comparison[i-3:i]
        
        fi.append(trio_actual[0] / trio_anterior[2])  # First value
        fi.append(trio_actual[1] / trio_anterior[1])  # Second value
        fi.append(trio_actual[2] / trio_anterior[0])  # Third value
    
    return fi


def initial_seeds(a):
    '''
    This function generates random initial conditions (seeds) to compute weight via optimization
    algorithm. The initial condition must be coherent with the afore mentioned constraints.

    return initial_guess
    
    '''
    initial_guess = np.zeros(3 * a)
    
    # Generates w[j, 1] (elemens 1, 4, 7, 10, ...)
    for j in range(a):
        initial_guess[3 * j + 1] = np.random.rand()  # From 0 to 1
    
    # Generates w[j, 0] (elements 0, 3, 6, 9, ...)
    for j in range(a):
        initial_guess[3 * j] = np.random.rand() * initial_guess[3 * j + 1]  # From 0 to w[j, 1]
    
    # Generates w[j, 2] (elements 2, 5, 8, 12, ...)
    for j in range(a):
        initial_guess[3 * j + 2] = initial_guess[3 * j + 1] + np.random.rand() * (1 - initial_guess[3 * j + 1])  # Entre w[j, 1] y 1
    
    # Normalize so sum h[j] = 1
    w_matrix = initial_guess.reshape((a, 3))
    h = np.sum((w_matrix[:, 0] + 4 * w_matrix[:, 1] + w_matrix[:, 2]) / 6)
    initial_guess /= h  
    
    return initial_guess

def transform_comparison(comparison):
    '''
    This function maps the linguistic label of comparison into triangular fuzzy numbers
    
    return comparison_num
    '''
    # Dictionary to map each type to its corresponding trio of values
    mapping = {
        'EI': (1, 1, 1),
        'WI': (2/3, 1, 3/2),
        'FI': (3/2, 2, 5/2),
        'VI': (5/2, 3, 7/2),
        'AI': (7/2, 4, 9/2)
    }
    
    # Initialize the comparison_num vector
    comparison_num = []
    
    # Convert each element in comparison to its corresponding trio
    for element in comparison:
        comparison_num.extend(mapping[element])  # Add the trio to comparison_num   
    return comparison_num


# Main function 
def weights_fuzzy_fucom(comparison_lin, num_seeds=100):

    '''
    Main function. It computes weights and fuzzy weights from fuzzy FUCOM algorithm. To do so,
    it applies SLSQP optimization method from different initial conditions. Chi is a metric 
    of the precision of the solution, the smaller, the better 

    return weigths, fuzzy_weights, chi
    
    '''
    comparison=transform_comparison(comparison_lin)
    fi = compute_fi(comparison)
    a = len(fi) // 3  # Number of criteria
    chi_max=0.1 #Precision threshold
    
    # Initialization to store the best solution
    best_w = None
    def_w = np.zeros(a)
    best_chi = np.inf
    
    # Applies different initial condition (seeds)
    for _ in range(num_seeds):
        
        initial_guess = initial_seeds(a)
        
        # Call the constraints
        constraints = [
            {'type': 'eq', 'fun': constraint_h, 'args': (a,)},
            {'type': 'ineq', 'fun': constraint_w0, 'args': (a,)},
            {'type': 'ineq', 'fun': constraint_w1, 'args': (a,)},
            {'type': 'ineq', 'fun': constraint_w2, 'args': (a,)}
        ]

        # Solve the optimization problem
        result = minimize(lambda w: objective(w, fi, chi_max), initial_guess, constraints=constraints, method='SLSQP')

        # Obtain fuzzy weights and precision metric
        w_matrix = result.x.reshape((a, 3))
        chi = objective(w_matrix.flatten(), fi, chi_max)

        # Updates the best solution
        if chi < best_chi:
            best_w = w_matrix
            best_chi = chi
    
    # Check the constraints
    for j in range(a - 1):
        diff1 = best_w[j, 0] - best_w[j+1, 2] * fi[3*j + 3]
        diff2 = best_w[j, 2] - best_w[j+1, 0] * fi[3*j + 5]
        diff3 = best_w[j, 1] - best_w[j+1, 1] * fi[3*j + 4]
        if abs(diff1) > chi_max:
            print(f"Constraint not satisfed: w[{j},0] - w[{j+1},2]*fi[{3*j + 3}] = {diff1} (must be ≤ {chi_max})")
        if abs(diff2) > chi_max:
            print(f"Constraint not satisfed: w[{j},2] - w[{j+1},0]*fi[{3*j + 5}] = {diff2} (must be ≤ {chi_max})")
        if abs(diff3) > chi_max:
            print(f"Constraint not satisfed: w[{j},1] - w[{j+1},1]*fi[{3*j + 4}] = {diff3} (must be ≤ {chi_max})")
        
    for j in range(a):    
        def_w[j]=(best_w[j, 0]+4*best_w[j, 1]+best_w[j, 2])/6

    # Print for testing
    print("TFN weight matrix:")
    print(best_w)
    print("Defuzzified weights:", def_w)
    print("Chi-metric:", best_chi)

    return def_w,best_w, best_chi   
