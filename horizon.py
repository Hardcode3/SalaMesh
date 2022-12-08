from mesh import Mesh
import numpy as np
import scipy.sparse
from scipy.sparse.linalg import lsmr
import importlib


model = "chevron"
m = Mesh(model+"/slice.obj")
attr = importlib.import_module(model + ".attributes")
nboundary = sum(m.boundary)

dico_corner_horizon = {}
list_corner_horizon = []
for c in range(m.ncorners): # count the number of horizon corner
    id = attr.horizon_id[c]
    if id>=0:
        list_corner_horizon.append(c)
        if id not in dico_corner_horizon:
            dico_corner_horizon[id] = []
        dico_corner_horizon[id].append(c)

ids_horizon = dico_corner_horizon.keys()
A = scipy.sparse.lil_matrix(( len(list_corner_horizon) + len(ids_horizon), m.nverts))
b = [0] * A.shape[0]

row = 0
list_i = []
for id_hor in ids_horizon:
    for corner in dico_corner_horizon[id_hor]:
        i = m.org(corner)
        list_i.append(i)
        A[row, i] =  1
        b[row] = id_hor
        row += 1
    i =  m.dst(dico_corner_horizon[id_hor][-1])
    list_i.append(i)
    A[row,i] =  1
    b[row] = id_hor
    row += 1
A = A.tocsr() # convert to compressed sparse row format for faster matrix-vector muliplications
x = lsmr(A, b)[0] # call the least squares solver
for i in list_i: # apply the computed flattening
        m.V[i][1] = x[i]

for c in range(m.ncorners): # lift all vertices of all horizons
    if attr.horizon_id[c]>=0:
        height = (1+attr.horizon_id[c]) # arbitrarily chosen coeff to get a visually nice result
        m.V[m.org(c)][2] = m.V[m.dst(c)][2] = height



print(A)
print(b)
print(x)
m.save("outout2.obj")