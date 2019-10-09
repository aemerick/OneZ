#
# code a search function here too
#

import numpy as np

MAX_STARS = 100000

cdef int find_bin(double[:] array, double x, int l, int r, int total):

    cdef int w   = r - l
    cdef int mid = l + w/2

    if(x<array[0]):
        return 0

    if(w<=1):
        if(mid < total-1):
            if ((x>=array[mid]) and (x<array[mid+1])):
                return mid
        if(mid>0):
            if((x>=array[mid-1]) and (x<array[mid])):
                return mid-1    
    if (array[mid] > x):
        return find_bin(array,x,l,mid-1,total)
    elif (array[mid]<x):
        return find_bin(array,x,mid+1,r,total)
    else:
        return mid

    return -1

def sample_imf(double[:] table, double M, 
               int N,
               double IMF_dm,
               double IMF_start, double IMF_end,
               int table_points):

    cdef double[:] stars = np.zeros(int(M/IMF_start))
    cdef int i = -1
    cdef double total_mass = 0.0
    cdef double rnum = 0.0
    cdef int bin = 0
    cdef double m_o = np.log10(IMF_start)

    while total_mass < M:
        i = i + 1
        rnum = np.random.rand()
        bin  = find_bin(table, rnum, 0, table_points, table_points)
        stars[i]      = 10.0**(bin * IMF_dm + m_o)
        total_mass    += stars[i]

    return stars
