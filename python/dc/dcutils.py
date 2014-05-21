#/usr/bin/python


def find_nearest(vector, value):
    '''
    Find value in \'vector\' nearest to \'value\'
    '''
    import numpy as np
    return np.abs(vector - value).argmin()
