import numpy as np

def weights_merec(decision_table,criterion_type):
    '''
    This funcion compute the weight using Merec Method 

    return merec weights
    '''
    from pyDecision.algorithm import merec_method
    
    weights=merec_method(decision_table,criterion_type)

    return weights