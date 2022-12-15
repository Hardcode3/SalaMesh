# SalaMesh
## Overview
A python geological mesh software to verticalize and horizontalize faults and horizons.
This work is based on the [least squares lecture of Dmitri Sokolov](https://github.com/ssloy/least-squares-course).
The work is done is described in the [following repository](https://github.com/ssloy/ENSG.git).

### Meshes
There are foud geological slices called [chevron](chevron), [ifp1](ifp1), [ifp2](ifp2) and [shell](shell).
These slices describe a 2D geological model having zero or more faults and some horizons.
Each mesh has different object files describing either the slice, the horizons or the set of faults.
The attributes.py file contains preprocessed data about horizon's and fault's edges.

### Files
