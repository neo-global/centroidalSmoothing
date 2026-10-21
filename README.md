This code consist of an rought quick & dirty implementation of excellent smoothing code from https://github.com/tkeskita/smoothMesh (all credits to original repo in regards of the algorithm itself).
This new implementation consist on enabling the smoothing with quality check inside snappyHexMesh for the existing MotionSolvers.
To use it, one needs to load the library in the controlDict file as:
```
libs (    "libnewSmoothSolver.so");
```
and in snappyHexMesh snap dictionary use:
```
solver newSmoothSolver;
newSmoothSolverCoeffs
    {
        centroidalIters 	1000;
	    // initEdgesFile		"constant/triSurface/initEdges.obj";
	    // targetEdgesFile		"constant/triSurface/initEdges.obj";
	    // targetSurfacesFile	"constant/triSurface/completeModel.stl";
    }
```
In the interest of maximum compatibility, and to utilize buffer layer snappying, we are currently developing and testing for ESI builds. Tested for v2512, 2606.
