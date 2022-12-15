from mesh import Mesh
import numpy as np
import scipy.sparse
from scipy.sparse.linalg import lsmr
import importlib


model = "chevron"
m = Mesh(model+"/slice.obj")
attr = importlib.import_module(model + ".attributes")
#nboundary = sum(m.boundary)

dico_corner_horizon = {}
list_corner_horizon = []
for c in range(m.ncorners): # pour decouper par horizon
    id = attr.horizon_id[c]
    if id>=0:
        list_corner_horizon.append(c)
        if id not in dico_corner_horizon:
            dico_corner_horizon[id] = []
        dico_corner_horizon[id].append(c)

nhor = len(list_corner_horizon)
ids_horizon = dico_corner_horizon.keys()


################# Remplissage matrice#################################
A = scipy.sparse.lil_matrix((m.ncorners+nhor + m.nverts, m.nverts))
b = [0] * A.shape[0]

for row in range(m.ncorners):  ###### beau maillage
    i = m.org(row)
    j = m.dst(row)
    A[row, j] =  1
    A[row, i] = -1
    b[row] = m.V[j][1] - m.V[i][1]

cnt = 0
for id_hor in range(len(ids_horizon)):  #####horizons
    i_fix = m.org(dico_corner_horizon[id_hor][0])
#    y_hor = m.V[i_fix][1] 
#    print(id_hor, '    ', i_fix, '  ', y_hor)
    for corner in dico_corner_horizon[id_hor]:
        i = m.org(corner)
        j = m.dst(corner)
        if j==i_fix: continue
        A[row, j] =  -100
        A[row, i_fix] = 100
        b[row] = 0
        row += 1
        cnt += 1
"""       
for v in range(m.nverts): #pour bloquer les composantes connexes
    A[row, v] = 1*0.1
    b[row] = m.V[v][1]*0.1
    row += 1 """


""" 
PROF
for c in range(m.ncorners): # count the number of horizon corner
    id = attr.horizon_id[c]
    if id<0: continue
    i = m.org(c)
    j = m.dst(c)
    A[row, j] =  10
    A[row, i] = -10
    row += 1 
for v in range(m.nverts): #pour bloquer les composantes connexes
    A[row, v] = 1*0.1
    b[row] = m.V[v][1]*0.1
    row += 1

    
    """






'''
dico_corner_horizon = {}
list_corner_horizon = []
for c in range(m.ncorners): # count the number of horizon corner
    id = attr.horizon_id[c]
    if id>=0:
        list_corner_horizon.append(c)
        if id not in dico_corner_horizon:
            dico_corner_horizon[id] = []
        dico_corner_horizon[id].append(c)

list_corner_fault = []
for c in range(m.ncorners): # lower vertices in faults
    if attr.is_fault[c]:
        list_corner_fault.append(c)

ids_horizon = dico_corner_horizon.keys()

#traitement des y
A = scipy.sparse.lil_matrix(( len(list_corner_horizon) + len(ids_horizon) + m.ncorners, m.nverts))
b = [0] * A.shape[0]

row = 0

#laplace
for row in range(m.ncorners):
    i = m.org(row)
    j = m.dst(row)
    A[row, j] =  1
    A[row, i] = -1


# Horizon
list_i = []
for id_hor in ids_horizon:
    i =  m.dst(dico_corner_horizon[id_hor][-1])
    y_hor = m.V[i][1] 
    list_i.append(i)
    A[row,i] =  1 *10
    b[row] = y_hor *10

    row += 1
    for corner in dico_corner_horizon[id_hor]:
        i = m.org(corner)
        list_i.append(i)
        A[row, i] =  1 * 10
        b[row] = y_hor * 10
        row += 1
    
A = A.tocsr() # convert to compressed sparse row format for faster matrix-vector muliplications
x = lsmr(A, b)[0] # call the least squares solver
for i in range(m.nverts): # apply the computed flattening
        m.V[i][1] = x[i]
'''

'''
#traitemetn des x
A = scipy.sparse.lil_matrix(( len(list_corner_fault) + m.ncorners + nboundary, m.nverts))
b = [0] * A.shape[0]

row = 0

#laplace
for row in range(m.ncorners):
    i = m.org(row)
    j = m.dst(row)
    A[row, j] =  1
    A[row, i] = -1


# Fault
list_i = []
for corner in list_corner_fault:
    i = m.org(corner)
    list_i.append(i)
    A[row, i] =  1 * 100
    b[row] = 1 * 100
    row += 1
# Limit
for (i,v) in enumerate(m.V):
    if m.on_border(i):
        A[row, i] = 1 *10 # quadratic penalty to lock boundary vertices
        b[row] = v[0] *10
        row += 1
    
A = A.tocsr() # convert to compressed sparse row format for faster matrix-vector muliplications
x = lsmr(A, b)[0] # call the least squares solver
for i in range(m.nverts): # apply the computed flattening
        m.V[i][0] = x[i]


for c in range(m.ncorners): # lift all vertices of all horizons
    if attr.horizon_id[c]>=0:
        height = (1+attr.horizon_id[c]) # arbitrarily chosen coeff to get a visually nice result
        m.V[m.org(c)][2] = m.V[m.dst(c)][2] = height

""" for c in range(m.ncorners): # lower vertices in faults
    if attr.is_fault[c]:
        m.V[m.org(c)][2] -= 0.00431 # arbitrary scaling coefficent
        m.V[m.dst(c)][2] -= 0.00431 # to make the result look nice """
'''

A = A.tocsr() # convert to compressed sparse row format for faster matrix-vector muliplications
x = lsmr(A, b)[0] # call the least squares solver
for i in range(m.nverts): # apply the computed flattening
    m.V[i][1] = x[i]

for c in range(m.ncorners): # lift all vertices of all horizons
    if attr.horizon_id[c]>=0:
        height = (1+attr.horizon_id[c])/100. # arbitrarily chosen coeff to get a visually nice result
        m.V[m.org(c)][2] = m.V[m.dst(c)][2] = height

for id_hor in range(len(ids_horizon)):  #####horizons
    for i in range(10):
        i_fix = m.org(dico_corner_horizon[id_hor][i])
        y_hor = m.V[i_fix][1] 
        print("    ", id_hor, '    ', i_fix, '  ', y_hor)
        break

m.save("evolution/9_fault_horizon.obj")