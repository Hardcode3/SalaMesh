from mesh import Mesh
import numpy as np
import scipy.sparse
from scipy.sparse.linalg import lsmr
import importlib


model = "ifp1"
m = Mesh(model+"/slice.obj")
attr = importlib.import_module(model + ".attributes")
#nboundary = sum(m.boundary)

# horizons
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

# failles
list_corner_fault = []
nb_fault_opposite = 0
for c in range(m.ncorners):
    if attr.is_fault[c]:
        list_corner_fault.append(c)
    if attr.fault_opposite[c] != -1:
        nb_fault_opposite += 1


################# Remplissage matrice pour Y#################################
A = scipy.sparse.lil_matrix((m.ncorners+nhor, m.nverts))
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

A = A.tocsr() # convert to compressed sparse row format for faster matrix-vector muliplications
y = lsmr(A, b)[0] # call the least squares solver
for i in range(m.nverts): # apply the computed flattening
    m.V[i][1] = y[i]

####################### Remplissage matrice pour X ###################
A = scipy.sparse.lil_matrix((m.ncorners + len(list_corner_fault) + nb_fault_opposite + m.nverts, m.nverts))
b = [0] * A.shape[0]

# smooth mesh
for row in range(m.ncorners):  ###### beau maillage
    i = m.org(row)
    j = m.dst(row)
    A[row, j] =  1
    A[row, i] = -1
    b[row] = m.V[j][0] - m.V[i][0]

# verticalize faults


for c_fault in list_corner_fault:
    i = m.org(c_fault)
    j = m.dst(c_fault)
    
    A[row, j] =  1 * 50
    A[row, i] = -1 * 50
    b[row] = 0

    row += 1


# join connex components
for c_fault in list_corner_fault:
    i = m.dst(c_fault)
    opposite_edge = attr.fault_opposite[c_fault]
    j = m.org(opposite_edge)
    A[row, i] =  1 * 1
    A[row, j] = -1 * 1
    b[row] = 0
    row += 1

''' 
for v in range(m.nverts): #pour bloquer les composantes connexes
    A[row, v] = 1*0.05
    b[row] = m.V[v][0]*0.05
    row += 1
'''

A = A.tocsr() # convert to compressed sparse row format for faster matrix-vector muliplications
x = lsmr(A, b, atol=1e-20, btol=1e-20)[0] # call the least squares solver
for i in range(m.nverts): # apply the computed flattening
    m.V[i][0] = x[i]

for c in range(m.ncorners): # lower vertices in faults
    if attr.is_fault[c]:
        m.V[m.org(c)][2] -= 0.01 # arbitrary scaling coefficent
        m.V[m.dst(c)][2] -= 0.01 # to make the result look nice

for c in range(m.ncorners): # lift all vertices of all horizons
    if attr.horizon_id[c]>=0:
        height = (1+attr.horizon_id[c])/100. # arbitrarily chosen coeff to get a visually nice result
        m.V[m.org(c)][2] = m.V[m.dst(c)][2] = height

""" for id_hor in range(len(ids_horizon)):  #####horizons
    for i in range(10):
        i_fix = m.org(dico_corner_horizon[id_hor][i])
        y_hor = m.V[i_fix][1] 
        print("    ", id_hor, '    ', i_fix, '  ', y_hor)
        break """

m.save("evolution/11_fault_horizon.obj")