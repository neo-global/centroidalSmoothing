#ifndef newSmoothSolverMotionSolver_H
#define newSmoothSolverMotionSolver_H

#include "displacementMotionSolver.H"
#include "volFields.H"
#include "fvMesh.H"
#include "fileName.H"
#include "word.H"
#include "string.H"
#include "labelList.H"

#include "argList.H"
#include "polyMesh.H"
#include "Time.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "ReadFields.H"
#include "syncTools.H"
#include "meshTools.H"
#include <bits/stdc++.h> // For std::stack
#include "vectorList.H"
#include "stringListOps.H"  // For stringListOps::findMatching();
#include "ListOps.H"
#include "wordReList.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "edgeMesh.H"
#include <fstream> // for fileExists

#include "smoothMeshCommon.H"
#include "boundaryPointSmoothing.C"  
#include "orthogonalBoundaryBlending.C"
namespace Foam
{

class newSmoothSolverMotionSolver
:
    public displacementMotionSolver
{
    // Private data

    // General/user-defined values
    int intValue_;
    fileName initEdgesFileName_;
    fileName targetEdgesFileName_;
    fileName targetSurfacesFileName_;
    scalar userScalar_;
    word userWord_;
    List<string> stringList_;

    // Smoothing parameters
    scalar minEdgeLength_;
    scalar maxStepLength_;
    scalar relStepFrac_;
    bool totalMinFreeze_;

    scalar minAngle_;
    scalar maxAngle_;

    bool edgeAngleConstraint_;
    bool faceAngleConstraint_;

    scalar layerMaxBlendingFraction_;
    scalar layerEdgeLength_;
    scalar layerExpansionRatio_;

    label minLayers_;
    label maxLayers_;

    scalar internalSmoothingBlendingFraction_;
    scalar relTol_;

    label centroidalIters_;
    label writeInterval_;

    scalar distanceTolerance_;

    // Private Member Functions

    static scalar meshMinEdgeLength(const polyMesh& mesh);

    void readCoeffs
    (
        const polyMesh& mesh,
        const dictionary& dict
    );


public:

    //- Runtime type information
    TypeName("newSmoothSolver");


    // Constructors

    //- Construct from mesh and dictionary
    newSmoothSolverMotionSolver
    (
        const polyMesh& mesh,
        const IOdictionary& dict
    );

    //- Construct from mesh, dictionary and displacement components
    newSmoothSolverMotionSolver
    (
        const polyMesh& mesh,
        const IOdictionary& dict,
        const pointVectorField& pointDisplacement,
        const pointIOField& points0
    );


    //- Destructor
    virtual ~newSmoothSolverMotionSolver() = default;


    // Member Functions

    //- Solve for motion
    virtual void solve();

    //- Return current points
    virtual tmp<pointField> curPoints() const;

    
/*Begin list of smoothMesh functions*/

// Help function to find internal points of the argument mesh.
// Updates argument boolList accordingly: true for internal points
// (including processor points); and false for boundary points.

int findInternalMeshPoints
(
    const fvMesh& mesh,
    boolList& isInternalPoint
);


// Function for centroidal smoothing of internal mesh points.

Foam::tmp<Foam::pointField> centroidalSmoothing
(
    const fvMesh& mesh,
    const boolList isInternalPoint,
    const bool doBoundarySmoothing
);


// Help function to return distance to current mesh point with index
// pointI from argument coords.

double getPointDistance
(
    const vector coords1,
    const vector coords2
);


/////////////////////////////////////////////
// Help functions for aspectRatioSmoothing //
/////////////////////////////////////////////

// Generate a lists of neighbour points from pointPoints where the
// point pair share a cell. Used in aspectRatioSmoothing.

int generatePointNeighPoints
(
    const fvMesh& mesh,
    labelListList& pointNeighPoints
);

// Help function which compares values of vectors element by element
// and returns true if first argument is smaller than second argument.

bool isSmallerByVectorElements
(
    const vector vec1,
    const vector vec2
);


// Help function to determine if point1 is closer to origin than
// point2. This function takes into account the special case where
// distance is equal but point coordinates are different.

bool isCloserPoint
(
    const vector point1,
    const vector point2
);

// Help function for aspectRatioSmoothing to find out Nth closest
// point index either among boundary points or all points

label findAppropriateClosestPointLabel
(
    const labelList& pointPoints,
    const labelList& sLabels,
    const label pointI,
    const boolList& isInternalPoint,
    const label stride
);

// Help function for aspectRatioSmoothing to find out closest
// three points (connected by edges); to each point

int findClosestPoints
(
    const fvMesh& mesh,
    const boolList& isInternalPoint,
    vectorList& closestPoints1,
    vectorList& closestPoints2,
    vectorList& closestPoints3,
    boolList& hasCommonCell,
    const labelListList& pointNeighPoints
);

// Help function for aspectRatioSmoothing to identify point locations
// of two short edge points which don't share cells.  Calculates also
// a blending fraction used for weighting the coordinates with
// centroidal smoothing coordinates. Blending fraction is set to zero
// if points share a cell, to prevent smoothing in that case.

double calcARSmoothingRatio
(
    const vector& closestPoint1,
    const vector& closestPoint2,
    const vector& closestPoint3,
    const bool hasCommonCell,
    const bool isInternalPoint
);

// Special smoothing algorithm for high aspect ratio cells. Blends
// resulting coordinates with centroidal point coordinates.

Foam::tmp<Foam::pointField> aspectRatioSmoothing
(
    const fvMesh& mesh,
    const boolList& isInternalPoint,
    const pointField& centroidalPoints,
    const labelListList& pointNeighPoints
);


// Quality control function which prohibits the decrease of edge
// length by freezing points to current locations, if edge length is
// below minEdgeLength and length is decreasing. This feature is used
// to limit the squishing or compression of cells near concave
// features.

int restrictEdgeShortening
(
    const fvMesh& mesh,
    pointField& origPoints,
    const double minEdgeLength,
    const bool totalMinFreeze,
    boolList& isFrozenPoint
);


// Calculate and return the maximum step length of the proposed new points

double getProposedMaxStepLength
(
    const fvMesh& mesh,
    const pointField& newPoints
);


// Quality control function to constrain the length of a step jump to
// new coordinates by an absolute length value. This increases the
// stability of the smoothing process, in case target coordinates are
// far off.

int constrainMaxStepLength
(
    const fvMesh& mesh,
    pointField& origPoints,
    const double maxStepLength,
    const double relStepFrac,
    const bool doGlobalScaling
);


/////////////////////////////////////////////////////
// Help functions for restrictMinEdgeAngleDecrease //
/////////////////////////////////////////////////////

// Calculate and return the edge-edge angle (in radians,
// 0 < angle < pi); of two edges which share a common point
// at coordinate cCoords. The two edge end point coordinates are
// p1Coords and p2Coords.

double edgeEdgeAngle
(
    const vector cCoords,
    const vector p1Coords,
    const vector p2Coords
);


// Find the neighbour mesh point indices neighPI1 and neighPI2 for the
// point with index pointI, where all points are part of face index
// faceI.

int getNeighbourPoints
(
    const fvMesh& mesh,
    const label centerPointI,
    const label faceI,
    label *neighPI1,
    label *neighPI2
);

// Calculate and store the minimum edge angles of point with index
// pointI for current mesh (minCAngleStorage); and a minimum for new
// mesh points (minNAngleStorage);.

int calc_min_edge_angles
(
    const fvMesh& mesh,
    const pointField& newPoints,
    const label pointI,
    double *minCAngleStorage,
    double *minNAngleStorage
);

// Restrict decrease of smallest edge-edge angles when angle is below
// minAngle (in degrees);. This is meant to avoid creation of
// self-intersections for concave features.

int restrictMinEdgeAngleDecrease
(
    const fvMesh& mesh,
    pointField& origPoints,
    const double minEdgeAngleInDegrees,
    boolList& isFrozenPoint
);

///////////////////////////////////////////////////////
// Help functions for restrictFaceAngleDeterioration //
///////////////////////////////////////////////////////

// Map face angle information from edges to points

int mapCurrentMinMaxFaceAnglesToPoints
(
    const fvMesh& mesh,
    const List<double>& minFaceAnglesForEdges,
    const List<double>& maxFaceAnglesForEdges,
    List<double>& minFaceAnglesForPoints,
    List<double>& maxFaceAnglesForPoints
);

// Calculate sum of angles between two edge unit vectors and a center
// unit vector. Assumes unit length for all vectors.

double calcEdgeCenterEdgeAngle
(
    const vector p0,
    const vector cC,
    const vector p1
);

// Calculate minimum and maximum face angles from the collection of
// precalculated point data projected on the edge normal plane

int calcMinMaxFinalProjectedAngle
(
    const label nCells,
    const vectorList& pVecs,
    const vectorList& cVecs,
    const labelList& faceIs,
    const labelList& f0Is,
    const labelList& f1Is,
    double& minFaceAngle,
    double& maxFaceAngle
);

// Find the pair of faces in face index list faceIs which are part of
// the cell cellI, and store the face label pair.

int findCellFacePair
(
    const fvMesh& mesh,
    const label cellI,
    const labelList& faceIs,
    const labelListList& cellFaces,
    label& f0I,
    label& f1I
);

// Calculate face center for face faceI as a weighted average of point coordinates.
// If point label pointI1 > -1, then that point is assumed to be at coords1.
// If point label pointI2 > -1, then that point is assumed to be at coords2.

vector calcFaceCenter
(
    const fvMesh& mesh,
    const label faceI,
    const label pointI1,
    const vector coords1,
    const label pointI2,
    const vector coords2
);

// Calculate minimum and maximum face angles for a single edge, with
// optional move of pointI1 to coords1 and pointI2 to coords2.

int calcMinMaxFaceAngleForEdge
(
    const fvMesh& mesh,
    const labelListList& cellFaces,
    const label edgeI,
    double& minFaceAngle,
    double& maxFaceAngle,
    const int pointI1,
    const vector coords1,
    const int pointI2,
    const vector coords2
);

// Wrapper for calcMinMaxEdgeAngle for current mesh points

int calcMinMaxFaceAngleForCurrentMeshEdge
(
    const fvMesh& mesh,
    const labelListList& cellFaces,
    const label edgeI,
    double& minFaceAngle,
    double& maxFaceAngle
);

// Calculate and store the minimum face angles for all mesh
// edges

int calcCurrentMinMaxFaceAnglesForEdges
(
    const fvMesh& mesh,
    const labelListList& cellFaces,
    List<double>& minFaceAnglesForEdges,
    List<double>& maxFaceAnglesForEdges
);

// Calculate and store minimum and maximum face angle for argument
// first point pointI1, with optional movement of pointI2 to given
// coordinates.

int calcMinMaxFaceAngleForPoint
(
    const fvMesh& mesh,
    const labelListList& cellFaces,
    const int pointI1,
    const vector coords1,
    const int pointI2,
    const vector coords2,
    double& minFaceAngle,
    double& maxFaceAngle
);


// Quality control function to restrict the decrease of smallest
// face-face angles when angle is below minAngle (in degrees);. This is
// meant to avoid creation of self-intersections. The face-face angle
// can be thought of a measure which quantifies squishing or folding
// of cell faces. In addition to prohibiting the deterioration of
// face-face angles of each point, it is also necessary to freeze
// movement of neighbouring points if their movement would deteriorate
// the angles at the current point.

int restrictFaceAngleDeterioration
(
    const fvMesh& mesh,
    const labelListList& cellFaces,
    pointField& origPoints,
    const double minFaceAngleInDegrees,
    const double maxFaceAngleInDegrees,
    boolList& isFrozenPoint
);


// Help function to get patch numbers from a given command line option

labelList getPatchIdsForOption
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& optionName,
    const bool defaultAll = false
);


// Calculate and return mesh statistics: the minimum and maximum edge
// length of the current mesh, and the total perimeter of the mesh,
// calculated as a sum of bounding box side lengths

int getMeshStats
(
    const fvMesh& mesh,
    double& minEdgeLength,
    double& maxEdgeLength,
    double& meshPerimeter
);

// Calculate the maximum relative change in the point positions, to be
// used as a measure of residual of smoothing

double calculateResidual
(
    const fvMesh& mesh,
    const pointField& newPoints,
    const boolList& isInternalPoint,
    const double maxStepLength
);


// Help function to generate lists of cell faces for all cells

int generateCellFaces
(
    const fvMesh& mesh,
    labelListList& cellFaces
);

// Help function to check does a file exists on the system

bool fileExists
(
    const std::string& fileName
);



// ************************************************************************* //


};

} // End namespace Foam

#endif

