Currently the surface smoothing points (https://github.com/tkeskita/smoothMesh#boundary-point-smoothing-options) has been disabled:
the reason behind this is that currently, if one uses `surfaceFeatureExtract` on the .stl file, what one gets is the 'red' edges 
of the image on the github (ie., edges at final position/edges to project to), but the solver will read it as the 'blue' edges 
of the image on the github (ie., edges to move/edges to project). Since during the snap step, the solver (centroidalSmoothing) 
is only used when using buffer layer, the buffer layer currently does not support feature edge snappying. if instead the 
solver (centroidalSmoothing) is only used for the BL smoothing (which from tkeskita descrition is still early developement 
https://github.com/tkeskita/smoothMesh#boundary-layer-treatment-related-options), one could use classic snappying step (without buffer)
and use the solver only on the boundary layer step. this needs testing 