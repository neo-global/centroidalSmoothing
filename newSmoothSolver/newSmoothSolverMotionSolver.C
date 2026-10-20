#include "newSmoothSolverMotionSolver.H"

#include "addToRunTimeSelectionTable.H"
#include "PstreamReduceOps.H"
#include "IOobject.H"
#include "dimensionedScalar.H"
#include "mathematicalConstants.H"

namespace Foam
{

defineTypeNameAndDebug(newSmoothSolverMotionSolver, 0);

addToRunTimeSelectionTable
(
    displacementMotionSolver,
    newSmoothSolverMotionSolver,
    dictionary 
);


addToRunTimeSelectionTable
(
    displacementMotionSolver,
    newSmoothSolverMotionSolver,
    displacement 
);

// * * * * * * * * Constructors * * * * * * * * //

newSmoothSolverMotionSolver::newSmoothSolverMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),

    intValue_(0),
    initEdgesFileName_("none"), 	//("constant/geometry/initEdges.obj"),
    targetEdgesFileName_("none"),	//("constant/geometry/targetEdges.obj"),
    targetSurfacesFileName_("none"),	//("constant/geometry/targetSurfaces.obj"),
    userScalar_(0),
    userWord_("none"),
    stringList_(),

    minEdgeLength_(0),
    maxStepLength_(0),
    relStepFrac_(0.5),
    totalMinFreeze_(false),

    minAngle_(35),
    maxAngle_(160),

    edgeAngleConstraint_(true),
    faceAngleConstraint_(true),

    layerMaxBlendingFraction_(0.3),
    layerEdgeLength_(0),
    layerExpansionRatio_(1.3),

    minLayers_(1),
    maxLayers_(4),

    internalSmoothingBlendingFraction_(0),
    relTol_(0.02),

    centroidalIters_(1000),
    writeInterval_(1000),

    distanceTolerance_(0)
{
    readCoeffs(mesh, coeffDict());
}


newSmoothSolverMotionSolver::newSmoothSolverMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointIOField& points0
)
:
    displacementMotionSolver
    (
        mesh,
        dict,
        pointDisplacement,
        points0,
        typeName
    ),

    intValue_(0),
    initEdgesFileName_("none"), 	//("constant/geometry/initEdges.obj"),
    targetEdgesFileName_("none"),	//("constant/geometry/targetEdges.obj"),
    targetSurfacesFileName_("none"),	//("constant/geometry/targetSurfaces.obj"),
    userScalar_(0),
    userWord_("none"),
    stringList_(),

    minEdgeLength_(0),
    maxStepLength_(0),
    relStepFrac_(0.5),
    totalMinFreeze_(false),

    minAngle_(35),
    maxAngle_(160),

    edgeAngleConstraint_(true),
    faceAngleConstraint_(true),

    layerMaxBlendingFraction_(0.3),
    layerEdgeLength_(0),
    layerExpansionRatio_(1.3),

    minLayers_(1),
    maxLayers_(4),

    internalSmoothingBlendingFraction_(0),
    relTol_(0.02),

    centroidalIters_(1000),
    writeInterval_(1000),

    distanceTolerance_(0)
{
    readCoeffs(mesh, coeffDict());
}


// * * * * * * * * Private functions * * * * * * * * //

scalar newSmoothSolverMotionSolver::meshMinEdgeLength
(
    const polyMesh& mesh
)
{
    const pointField& pts = mesh.points();
    const edgeList& edges = mesh.edges();

    scalar localMin = GREAT;

    forAll(edges, edgei)
    {
        const edge& e = edges[edgei];

        const scalar edgeLength =
            mag(pts[e[1]] - pts[e[0]]);

        localMin = min(localMin, edgeLength);
    }

    // Global minimum in parallel
    reduce(localMin, minOp<scalar>());

    if (localMin == GREAT)
    {
        FatalErrorInFunction
            << "Unable to calculate the minimum mesh edge length."
            << exit(FatalError);
    }

    return localMin;
}


void newSmoothSolverMotionSolver::readCoeffs
(
    const polyMesh& mesh,
    const dictionary& dict
)
{
    const scalar meshMin = meshMinEdgeLength(mesh);

    // User-defined data members
    intValue_ = dict.lookupOrDefault<label>("intValue", 0);

    initEdgesFileName_ =
        dict.lookupOrDefault<fileName>
        (
            "initEdgesFile",
            "none" //"constant/geometry/initEdges.obj"
        );

    targetEdgesFileName_ =
        dict.lookupOrDefault<fileName>
        (
            "targetEdgesFile",
            "none" //"constant/geometry/targetEdges.obj"
        );

    targetSurfacesFileName_ =
        dict.lookupOrDefault<fileName>
        (
            "targetSurfacesFile",
            "none" //"constant/geometry/targetSurfaces.obj"
        );

    userScalar_ =
        dict.lookupOrDefault<scalar>("userScalar", 0);

    userWord_ =
        dict.lookupOrDefault<word>("userWord", "none");

    stringList_ =
        dict.lookupOrDefault<List<string>>
        (
            "stringList",
            List<string>()
        );


    // Smoothing parameters
    minEdgeLength_ =
        dict.lookupOrDefault<scalar>
        (
            "minEdgeLength",
            0.5*meshMin
        );

    maxStepLength_ =
        dict.lookupOrDefault<scalar>
        (
            "maxStepLength",
            0.3*minEdgeLength_
        );

    if (maxStepLength_ > 0.5*minEdgeLength_)
    {
        WarningInFunction
            << "The maximum allowed step length is more than "
            << "half of the minimum edge length. This may cause "
            << "instability in smoothing."
            << nl << endl;
    }

    relStepFrac_ =
        dict.lookupOrDefault<scalar>("relStepFrac", 0.5);

    totalMinFreeze_ =
        dict.lookupOrDefault<bool>("totalMinFreeze", false);

    minAngle_ =
        dict.lookupOrDefault<scalar>("minAngle", 35.0);

    maxAngle_ =
        dict.lookupOrDefault<scalar>("maxAngle", 160.0);

    edgeAngleConstraint_ =
        dict.lookupOrDefault<bool>("edgeAngleConstraint", true);

    faceAngleConstraint_ =
        dict.lookupOrDefault<bool>("faceAngleConstraint", true);

    layerMaxBlendingFraction_ =
        dict.lookupOrDefault<scalar>
        (
            "layerMaxBlendingFraction",
            0.3
        );

    layerEdgeLength_ =
        dict.lookupOrDefault<scalar>
        (
            "layerEdgeLength",
            minEdgeLength_
        );

    layerExpansionRatio_ =
        dict.lookupOrDefault<scalar>
        (
            "layerExpansionRatio",
            1.3
        );

    minLayers_ =
        dict.lookupOrDefault<label>("minLayers", 1);

    maxLayers_ =
        dict.lookupOrDefault<label>("maxLayers", 4);

    internalSmoothingBlendingFraction_ =
        dict.lookupOrDefault<scalar>
        (
            "internalSmoothingBlendingFraction",
            0.0
        );

    relTol_ =
        dict.lookupOrDefault<scalar>("relTol", 0.02);

    centroidalIters_ =
        dict.lookupOrDefault<label>("centroidalIters", 1000);

    writeInterval_ =
        dict.lookupOrDefault<label>
        (
            "writeInterval",
            centroidalIters_
        );

    distanceTolerance_ =
        REL_TOL*min(meshMin, layerEdgeLength_);

    // Use the intended relative tolerance, if REL_TOL is defined
    // in your project, instead of the mathematical small constant.
    Info<< typeName << " coefficients:" << nl
        << "    meshMinEdgeLength = " << meshMin << nl
        << "    minEdgeLength = " << minEdgeLength_ << nl
        << "    maxStepLength = " << maxStepLength_ << nl
        << "    relStepFrac = " << relStepFrac_ << nl
        << "    totalMinFreeze = " << totalMinFreeze_ << nl
        << "    minAngle = " << minAngle_ << nl
        << "    maxAngle = " << maxAngle_ << nl
        << "    edgeAngleConstraint = " << edgeAngleConstraint_ << nl
        << "    faceAngleConstraint = " << faceAngleConstraint_ << nl
        << "    layerEdgeLength = " << layerEdgeLength_ << nl
        << "    minLayers = " << minLayers_ << nl
        << "    maxLayers = " << maxLayers_ << nl
        << "    centroidalIters = " << centroidalIters_ << nl
        << "    writeInterval = " << writeInterval_ << nl
        << "    distanceTolerance = " << distanceTolerance_
        << endl;
}


// * * * * * * * * Member functions * * * * * * * * //

void newSmoothSolverMotionSolver::solve()
{
    movePoints(curPoints());
    //-FO
    fvMesh& fvM = const_cast<fvMesh&>(refCast<const fvMesh>(this->mesh()));
    const fvMesh& mesh = fvM;
    //-FO

    double meshMinEdgeLength = 0, meshMaxEdgeLength = 0, meshPerimeter = 0;
    getMeshStats(mesh, meshMinEdgeLength, meshMaxEdgeLength, meshPerimeter);

    const labelList layerPatchIds = getPatchIdsForOption(mesh, coeffDict(), "layerPatches", false);
    const labelList smoothingPatchIds = getPatchIdsForOption(mesh, coeffDict(), "smoothingPatches", true);
    const fileName initEdgesFileName(initEdgesFileName_);
    const fileName targetEdgesFileName(targetEdgesFileName_);
    const fileName targetSurfacesFileName(targetSurfacesFileName_);
    const string initEdgesFileString(initEdgesFileName_.c_str());
    const string targetEdgesFileString(targetEdgesFileName_.c_str());
    const string targetSurfacesFileString(targetSurfacesFileName_.c_str());
    bool doLayerTreatment = false;
    bool doBoundarySmoothing = false;
    pointField finalPoints(mesh.points());

    // Storage for markers for internal points
    boolList isInternalPoint(mesh.nPoints(), false);
    const label nInternalPoints = findInternalMeshPoints(mesh, isInternalPoint);

    // Storage for number of edge hops to reach layer and free
    // boundaries for mesh points (for boundary layer treatment)
    labelList pointHopsToLayerBoundary(mesh.nPoints(), UNDEF_LABEL);
    labelList pointHopsToSmoothingBoundary(mesh.nPoints(), UNDEF_LABEL);

    // Storage for point normals (for boundary layer treatment)
    tmp<pointField> tPointNormals(new pointField(mesh.nPoints(), Zero));
    pointField& pointNormals = tPointNormals.ref();

    // Storage for neighbour point locations (for boundary layer treatment)
    tmp<pointField> tInnerNeighCoords(new pointField(mesh.nPoints(), UNDEF_VECTOR));
    pointField& innerNeighCoords = tInnerNeighCoords.ref();
    tmp<pointField> tOuterNeighCoords(new pointField(mesh.nPoints(), UNDEF_VECTOR));
    pointField& outerNeighCoords = tOuterNeighCoords.ref();

    // Storage for marking neighbour point being inside same processor
    // domain (for boundary layer treatment)
    boolList isInnerNeighInProc(mesh.nPoints(), false);
    boolList isOuterNeighInProc(mesh.nPoints(), false);

    // Storage for index map from point to neighbour point inside same
    // processor domain. One map points towards inner mesh, and
    // another map towards outer boundary. (for boundary layer treatment)
    labelList pointToInnerPointMap(mesh.nPoints(), UNDEF_LABEL);
    labelList pointToOuterPointMap(mesh.nPoints(), UNDEF_LABEL);

    // Boolean list for marking frozen points (points not allowed to
    // move during smoothing). This list is synced among processors.
    boolList isFrozenPoint(mesh.nPoints(), false);

    // A list of point label lists to indicate cell sharing. Used in
    // aspectRatioSmoothing.
    Info << "Starting to build pointNeighPoints (this may take some time)" << endl;
    labelListList pointNeighPoints(mesh.nPoints());
    generatePointNeighPoints(mesh, pointNeighPoints);
    Info << "Done building pointNeighPoints" << endl << endl;

    // Lists of all cell faces (includes both internal and boundary
    // faces) for all cells
    labelListList cellFaces(mesh.nCells());
    generateCellFaces(mesh, cellFaces);

    // Check prerequisites for carrying out boundary layer treatment
    doLayerTreatment = false;
    if ((layerPatchIds.size() > 0) and (layerMaxBlendingFraction_ > SMALL))
    {
        doLayerTreatment = true;
        Info << "Enabled boundary layer treatment" << endl << endl;
    }
    else
    {
        Info << "Boundary layer treatment is disabled. Either no layerPatches were specified or boundaryMaxBlendingFraction is zero" << endl << endl;
    }

    // IOLists for communication of isCornerPoint and
    // isFeatureEdgePoint between consequent runs.
    // Using labelIOList because there is no boolIOList.

    labelIOList isCornerPointIO
    (
        IOobject
        (
            "isCornerPoint",
            fvM.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            true
        ),
        labelList(mesh.nPoints(), 0)
    );

    labelIOList isFeatureEdgePointIO
    (
        IOobject
        (
            "isFeatureEdgePoint",
            fvM.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            true
        ),
        labelList(mesh.nPoints(), 0)
    );

    // Determine do the labelIOLists contain classification data
    bool labelIOListsHaveData = false;
    if ((findIndex(isCornerPointIO, 1) >= 0) or (findIndex(isFeatureEdgePointIO, 1) >= 0))
    {
        labelIOListsHaveData = true;
        Info << "Found corners and feature edges in isCornerPoint and isFeatureEdgePoint files" << endl << endl;
    }
    else
    {
        Info << "Did not find corners and feature edges in isCornerPoint and isFeatureEdgePoint files" << endl << endl;
    }

    // Check prerequisites for carrying out boundary point smoothing
    doBoundarySmoothing = false;
    bool noEdgeMesh	= false;
    //-RM: searches for both initEdges and targetSurfaces
    if ((fileExists(targetSurfacesFileString)) and
        ((fileExists(initEdgesFileString)) or (labelIOListsHaveData)) and
        (smoothingPatchIds.size() > 0))
    {
        doBoundarySmoothing = true;
        Info << "Enabled boundary point smoothing" << endl << endl;
    }
    else if ((fileExists(targetSurfacesFileString)) and
        !((fileExists(initEdgesFileString)) or (labelIOListsHaveData)) and
        (smoothingPatchIds.size() > 0))
    {
	
    	//-RM: placeholder of workaround to skip initEdges.
	//works if initEdges file is missing. 
        doBoundarySmoothing = true;
	noEdgeMesh	    = true;
        Info << "Enabled boundary point smoothing without edge mesh." << endl << endl;
    }
    else
    {
        Info << "Boundary point smoothing is disabled. Missing smoothingPatches, or one or both of files:" << endl
             << targetSurfacesFileString << endl
             << initEdgesFileString << endl << endl;
    }

    if ((doLayerTreatment) and (! doBoundarySmoothing))
    {
        Info << "WARNING: Boundary layer treatment will be done without boundary point smoothing. This can result in distorted boundary cells." << endl << endl;
    }

    // Objects for boundary point snapping to surfaces
    autoPtr<triSurface> surf(nullptr);
    autoPtr<triSurfaceSearch> searchSurfaces(nullptr);
    autoPtr<indexedOctree<treeDataTriSurface>> tree(nullptr);

    // Objects for boundary point snapping to feature edges
    autoPtr<edgeMesh> initEdges(nullptr);
    autoPtr<edgeMesh> targetEdges(nullptr);

    // Point classification lists
    boolList isConnectedToInternalPoint(mesh.nPoints(), false);
    boolList isFeatureEdgePoint(mesh.nPoints(), false);
    boolList isLayerSurfacePoint(mesh.nPoints(), false);
    boolList isSmoothingSurfacePoint(mesh.nPoints(), false);
    boolList isFrozenSurfacePoint(mesh.nPoints(), false);
    boolList isProcessorPoint(mesh.nPoints(), false);
    boolList isCornerPoint(mesh.nPoints(), false);
    boolList isSharpEdgePoint(mesh.nPoints(), false);

    // Storage for target corner point coordinates
    vectorList cornerPoints(mesh.nPoints(), UNDEF_VECTOR);

    // Storage for saving unique edge string indices of feature edges
    // for edges
    labelList targetEdgeStrings;

    // Storage for saving unique edge string indices of feature edges
    // for points
    labelList pointStrings(mesh.nPoints(), UNDEF_LABEL);

    // Preparations for boundary point smoothing
    if (doBoundarySmoothing)
    {
        // Target surface mesh, build search tree
        surf.reset(new triSurface(targetSurfacesFileName));
        Info << "Target surfaces file " << targetSurfacesFileName << " stats:" << endl;
        surf().writeStats(Info);
        Info << endl;

        searchSurfaces.reset(new triSurfaceSearch(surf));
        tree.reset(new indexedOctree<treeDataTriSurface>(searchSurfaces().tree()));

      if(noEdgeMesh)
      {
            // Initialize empty edge mesh, to avoid issues in
            // classifyBoundaryPoints()
            initEdges.reset(new edgeMesh());
            targetEdges.reset(new edgeMesh());
      }
      else
      {
        // Initial and target feature edge meshes
        initEdges.reset(new edgeMesh(initEdgesFileName));
        Info << "Initial feature edges file " << initEdgesFileName << " stats:" << endl;
        initEdges().writeStats(Info);
        Info << endl;

        if (fileExists(targetEdgesFileString))
        {
            targetEdges.reset(new edgeMesh(targetEdgesFileName));
            Info << "Target feature edges file " << targetEdgesFileName << " stats:" << endl;
            targetEdges().writeStats(Info);
        }
        else
        {
            targetEdges.reset(new edgeMesh(initEdgesFileName));
            Info << "WARNING: Initial feature edges will be used also as target edges, because" << endl

            << "did not find file " << targetEdgesFileString << "." << endl << endl;
        }
	// Check sanity of edge meshes
        Info << "Checking initial edge mesh sanity" << endl;
        checkEdgeMeshSanity(initEdges, meshMinEdgeLength, meshPerimeter);
        Info << "Checking target edge mesh sanity" << endl;
        checkEdgeMeshSanity(targetEdges, meshMinEdgeLength, meshPerimeter);

        // Generate indices for target edge mesh edge strings
        Info << endl << "Starting to build targetEdgeStrings" << endl;
        const label nStrings = findEdgeMeshStrings(targetEdgeStrings, targetEdges);
        Info << "Detected number of target edge mesh strings: " << nStrings + 1 << endl << endl;
      }
    }
    else
    {
        // Initialize empty edge mesh, to avoid issues in
        // classifyBoundaryPoints()
        initEdges.reset(new edgeMesh());
        targetEdges.reset(new edgeMesh());
    }

    const label nPoints = returnReduce(mesh.nPoints(), sumOp<label>());
    Info << "Mesh includes a total of " << nPoints << " points:" << endl
         << "  - " << nInternalPoints << " internal (non-boundary) points" << endl
         << "  - " << nPoints - nInternalPoints << " boundary points" << endl
         << "Mesh minimum edge length = " << meshMinEdgeLength << endl
         << "Mesh maximum edge length = " << meshMaxEdgeLength << endl
         << "Distance tolerance = " << distanceTolerance_ << endl << endl;

    // Classify boundary points and find target corner points for feature edge snapping
    classifyBoundaryPoints
    (
        mesh,
        initEdges,
        targetEdges,
        layerPatchIds,
        smoothingPatchIds,
        isInternalPoint,
        isProcessorPoint,
        isConnectedToInternalPoint,
        isFeatureEdgePoint,
        isFeatureEdgePointIO,
        isCornerPoint,
        isCornerPointIO,
        labelIOListsHaveData,
        cornerPoints,
        isLayerSurfacePoint,
        isSmoothingSurfacePoint,
        isFrozenSurfacePoint,
        targetEdgeStrings,
        doBoundarySmoothing,
        distanceTolerance_
    );

    // Preparations for optional smoothing and treatment
    if ((doBoundarySmoothing) or (doLayerTreatment))
    {
        calculatePointHopsToBoundary(mesh, layerPatchIds, isInternalPoint, isConnectedToInternalPoint, pointHopsToLayerBoundary, maxLayers_ + 1);
        calculatePointHopsToBoundary(mesh, smoothingPatchIds, isInternalPoint, isConnectedToInternalPoint, pointHopsToSmoothingBoundary, 2);
        calculateBoundaryPointNormals(mesh, pointNormals, isSharpEdgePoint);
        propagateOuterNeighInfo(mesh, isInternalPoint, isLayerSurfacePoint, isOuterNeighInProc, pointToOuterPointMap, pointNormals, pointHopsToLayerBoundary, maxLayers_ + 1);
        propagateInnerNeighInfo(mesh, isSmoothingSurfacePoint, isConnectedToInternalPoint, isInnerNeighInProc, pointToInnerPointMap, pointHopsToSmoothingBoundary);

        forAll (mesh.points(), pointI)
        {
            if (isSharpEdgePoint[pointI] == true)
            {
                Pout << "WARNING: Froze a sharp boundary edge point which is not on a feature edge: pointI " << pointI << " at " << mesh.points()[pointI] << endl;
            }
        }
    }

    // Find target edge string indices for points (used in feature
    // edge snapping)
    if (doBoundarySmoothing)
    {
        forAll(mesh.points(), pointI)
        {
            if (! isFeatureEdgePoint[pointI])
            {
                continue;
            }

            point dummyPoint;
            label dummy, dummy2;
            label pointStringI = UNDEF_LABEL;
            findClosestEdgeInfo(mesh.points()[pointI], targetEdges, -1, targetEdgeStrings, distanceTolerance_, dummyPoint, dummy, pointStringI, dummy2);
            pointStrings[pointI] = pointStringI;
        }
    }

    // Carry out smoothing iterations
    // ------------------------------

    bool stopIteration = false;

    for (label i = 0; i < centroidalIters_; ++i)
    {
        fvM. movePoints(curPoints());
        // Info << "Starting iteration " << i << endl;

        // Reset frozen points
        forAll(isFrozenPoint, pointI)
            isFrozenPoint[pointI] = false;

        // Recalculate point normals
        calculateBoundaryPointNormals(mesh, pointNormals, isSharpEdgePoint);

        // Calculate new point locations using centroidal smoothing
        tmp<pointField> tCentroidalPoints = centroidalSmoothing(mesh, isInternalPoint, doBoundarySmoothing);
        pointField& centroidalPoints = tCentroidalPoints.ref();

        // Constrain absolute length of jump to new coordinates, to stabilize smoothing
        // constrainMaxStepLength(mesh, centroidalPoints, maxStepLength_, relStepFrac_, false);

        // Blend centroidal points with points from aspect ratio smoothing
        tmp<pointField> tNewPoints = aspectRatioSmoothing(mesh, isInternalPoint, centroidalPoints, pointNeighPoints);
        pointField& newPoints = tNewPoints.ref();

        // Constrain absolute length of jump to new coordinates, to stabilize smoothing
        constrainMaxStepLength(mesh, newPoints, maxStepLength_, relStepFrac_, false);

        // Optional boundary layer treatment
        if (doLayerTreatment)
        {
            // Update neighbour coordinates and synchronize among processors
            updateNeighCoords(mesh, isOuterNeighInProc, pointToOuterPointMap, outerNeighCoords);
            // Blend orthogonal and centroidal coordinates to newPoints
            blendWithOrthogonalPoints
            (
                 mesh,
                 newPoints,
                 isInternalPoint,
                 pointHopsToLayerBoundary,
                 pointNormals,
                 outerNeighCoords,
                 layerMaxBlendingFraction_,
                 layerEdgeLength_,
                 layerExpansionRatio_,
                 minLayers_,
                 maxLayers_ + 1  // +1 for correct number of layers
            );

            // Constrain absolute length of jump to new coordinates, to stabilize smoothing
            constrainMaxStepLength(mesh, newPoints, maxStepLength_, relStepFrac_, false);
        }

        if (doBoundarySmoothing)
        {
            // Update neighbour coordinates and synchronize among processors
            updateNeighCoords(mesh, isInnerNeighInProc, pointToInnerPointMap, innerNeighCoords);
            // Project boundary points
            projectBoundaryPointsToEdgesAndSurfaces
            (
                mesh,
                newPoints,
                pointNormals,
                isInternalPoint,
                isSmoothingSurfacePoint,
                isFeatureEdgePoint,
                isCornerPoint,
                cornerPoints,
                targetEdges,
                surf,
                tree,
                meshMinEdgeLength,
                targetEdgeStrings,
                pointStrings,
                isSharpEdgePoint,
                distanceTolerance_,
                isFrozenPoint
            );

            // Constrain absolute length of jump to new coordinates, to stabilize smoothing
            // constrainMaxStepLength(mesh, newPoints, maxStepLength_, relStepFrac_, false);

            // Use the locations of first cell layer points for
            // projecting points to boundary surfaces
            projectPrismaticInternalPointsToSurfaces
            (
                mesh,
                newPoints,
                pointHopsToSmoothingBoundary,
                pointNormals,
                isSmoothingSurfacePoint,
                isConnectedToInternalPoint,
                isFeatureEdgePoint,
                isCornerPoint,
                pointToInnerPointMap,
                innerNeighCoords,
                internalSmoothingBlendingFraction_,
                isSharpEdgePoint
            );

            // Constrain absolute length of jump to new coordinates, to stabilize smoothing
            constrainMaxStepLength(mesh, newPoints, maxStepLength_, relStepFrac_, false);
        }

        // Avoid shortening of short edge length
        restrictEdgeShortening(mesh, newPoints, minEdgeLength_, totalMinFreeze_, isFrozenPoint);

        if (edgeAngleConstraint_)
        {
            // Restrict decrease of smallest edge-edge angle
            restrictMinEdgeAngleDecrease(mesh, newPoints, minAngle_, isFrozenPoint);
        }

        if (faceAngleConstraint_)
        {
            // Restrict deterioration of face-face angles
            restrictFaceAngleDeterioration(mesh, cellFaces, newPoints, minAngle_, maxAngle_, isFrozenPoint);
        }

        // Synchronize and combine the list of frozen points
        syncTools::syncPointList
        (
            mesh,
            isFrozenPoint,
            orEqOp<bool>(),
            false               // null value
        );

        // Restore original coordinates for frozen points and boundary
        // points which are not to be smoothed.
        label nFrozenPoints = 0;
        forAll(newPoints, pointI)
        {
            if
            (
                isFrozenPoint[pointI]
            || (!isInternalPoint[pointI] && !isSmoothingSurfacePoint[pointI])
            )
            {
                newPoints[pointI] = mesh.points()[pointI];
                ++nFrozenPoints;
            }
        }

        const scalar res = calculateResidual
        (
            mesh, newPoints, isInternalPoint, maxStepLength_
        );
        Info<< "Smoothing iteration=" << (i + 1)
            << " nFrozenPoints="
            << returnReduce(nFrozenPoints, sumOp<label>())
            << " residual=" << res << nl;

        // The motion solver must not move/write the mesh itself. Store
        // the candidate coordinates and let the motion framework apply
        // them through curPoints().
    
        //-FO
        finalPoints = newPoints;

	//-FO
        pointVectorField& disp = pointDisplacement();
        forAll(finalPoints, pointI)
        {
            disp.primitiveFieldRef()[pointI] = finalPoints[pointI] - points0()[pointI];
        }

        fvM.movePoints(curPoints());
        if (i == centroidalIters_ - 1)
        {
            Info<< "Maximum centroidalIters reached, stopping." << nl;
            break;
        }

    }
    
 	// Encode final coordinates as displacement from the reference mesh.

/*
	pointVectorField& disp = const_cast<pointVectorField&>(pointDisplacement_);
    	forAll(finalPoints, pointI)
    	{
            disp.primitiveFieldRef()[pointI] = finalPoints[pointI] - points0()[pointI];
    	}
*/
}



tmp<pointField>
newSmoothSolverMotionSolver::curPoints() const
{
    tmp<pointField> tPoints(new pointField(points0()));
    pointField& pts = tPoints.ref();
    const vectorField& disp = pointDisplacement_.primitiveField();
    forAll(pts, pointI)
    {
        pts[pointI] += disp[pointI];
    }
    return tPoints;
}

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

// Help function to find internal points of the argument mesh.
// Updates argument boolList accordingly: true for internal points
// (including processor points) and false for boundary points.

int newSmoothSolverMotionSolver::findInternalMeshPoints
(
    const fvMesh& mesh,
    boolList& isInternalPoint
)
{
    // Start from all points in
    isInternalPoint = true;

    // Remove points on boundary patches, except not
    // processor patches nor empty patches
    const faceList& faces = mesh.faces();
    forAll(mesh.boundary(), patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        // Skip processor and empty patches
        if (isA<processorPolyPatch>(pp))
            continue;

        // No smoothing for non-3D meshes
        if (isA<emptyPolyPatch>(pp))
        {
            FatalError
                << "Smoothing of non-3D meshes (meshes with type empty patches) is not supported"
                << endl << abort(FatalError);
        }

        const label startI = mesh.boundary()[patchI].start();
        const label endI = startI + mesh.boundary()[patchI].Cf().size();

        for (label faceI = startI; faceI < endI; faceI++)
        {
            const face& f = mesh.faces()[faceI];
            forAll (f, pointI)
            {
                const label i = faces[faceI][pointI];
                isInternalPoint[i] = false;
            }
        }
    }

    label nPoints = 0;
    forAll(isInternalPoint, pointI)
    {
        if (isInternalPoint[pointI])
            ++nPoints;
    }

    const label sumNPoints = returnReduce(nPoints, sumOp<label>());
    return sumNPoints;
}


// Function for centroidal smoothing of internal mesh points.

Foam::tmp<Foam::pointField> newSmoothSolverMotionSolver::centroidalSmoothing
(
    const fvMesh& mesh,
    const boolList isInternalPoint,
    const bool doBoundarySmoothing
)
{
    // Centroidal smoothing algorithm

    // Calculate number and sum of surrounding cell center
    // coordinates

    vectorField cellPoints(mesh.points().size(), Zero);

    labelList nPoints(mesh.points().size(), Zero);

    forAll(mesh.points(), pointI)
    {
        // Carry out centroidal smoothing only for internal points
        // if boundary point smoothing is disabled
        if ((! doBoundarySmoothing) and (! isInternalPoint[pointI]))
        {
            continue;
        }

        const labelList& pCells = mesh.pointCells(pointI);

        // Number of points
        nPoints[pointI] = pCells.size();

        // Sum of cell center coordinates
        for (const label celli : pCells)
        {
            cellPoints[pointI] += mesh.cellCentres()[celli];
        }
    }

    // Synchronize among processors
    syncTools::syncPointList
    (
        mesh,
        cellPoints,
        plusEqOp<vector>(),
        UNDEF_VECTOR               // null value
    );

    syncTools::syncPointList
    (
        mesh,
        nPoints,
        plusEqOp<label>(),
        label(0)                   // null value
    );

    // Calculate new point locations
    // - init with fallback value (the current mesh points)
    tmp<pointField> tnewPoints(new pointField(mesh.points()));
    auto& newPoints = tnewPoints.ref();

    forAll(newPoints, pointI)
    {
        // centroidal point
        if (nPoints[pointI])
        {
            newPoints[pointI] =
                cellPoints[pointI] / nPoints[pointI];
        }
    }

    return tnewPoints;
}


// Help function to return distance to current mesh point with index
// pointI from argument coords.

double newSmoothSolverMotionSolver::getPointDistance
(
    const vector coords1,
    const vector coords2
)
{
    const vector v = coords2 - coords1;
    return mag(v);
}


/////////////////////////////////////////////
// Help functions for aspectRatioSmoothing //
/////////////////////////////////////////////

// Generate a lists of neighbour points from pointPoints where the
// point pair share a cell. Used in aspectRatioSmoothing.

int newSmoothSolverMotionSolver::generatePointNeighPoints
(
    const fvMesh& mesh,
    labelListList& pointNeighPoints
)
{
    forAll (mesh.points(), pointI)
    {
        forAll (mesh.pointCells(pointI), pointCellI)
        {
            const label cellI = mesh.pointCells(pointI)[pointCellI];
            forAll (mesh.cellPoints(cellI), pointCPI)
            {
                const label pointPointI = mesh.cellPoints(cellI)[pointCPI];
                if (pointI == pointPointI)
                    continue;

                // Add if not already added
                if (findIndex(pointNeighPoints[pointI], pointPointI) == -1)
                {
                    pointNeighPoints[pointI].append(pointPointI);
                }
            }
        }
    }

    return 0;
}

// Help function which compares values of vectors element by element
// and returns true if first argument is smaller than second argument.

bool newSmoothSolverMotionSolver::isSmallerByVectorElements
(
    const vector vec1,
    const vector vec2
)
{
    if (vec1.size() != vec2.size())
        FatalError << vec1 << " size does not equal size of " << vec2 << endl << abort(FatalError);

    forAll(vec1, i)
    {
        if (vec1[i] < vec2[i])
            return true;
        else if (vec1[i] > vec2[i])
            return false;
    }
    return false;
}


// Help function to determine if point1 is closer to origin than
// point2. This function takes into account the special case where
// distance is equal but point coordinates are different.

bool newSmoothSolverMotionSolver::isCloserPoint
(
    const vector point1,
    const vector point2
)
{
    if (point1 == point2)
    {
        return false;
    }

    // point1 is clearly closer than point2
    const double deltaDistance = mag(point1) - mag(point2);
    if (deltaDistance < VSMALL)
    {
        return true;
    }

    // point1 and point2 are at same distance, prefer the point with
    // smaller x, y, or z coordinate
    else if ((abs(deltaDistance) < VSMALL) and (isSmallerByVectorElements(point1, point2)))
    {
        return true;
    }

    return false;
}

// Help function for aspectRatioSmoothing to find out Nth closest
// point index either among boundary points or all points

label newSmoothSolverMotionSolver::findAppropriateClosestPointLabel
(
    const labelList& pointPoints,
    const labelList& sLabels,
    const label pointI,
    const boolList& isInternalPoint,
    const label stride
)
{
    const bool isThisInternalPoint = isInternalPoint[pointI];
    label counter = 0;

    forAll(sLabels, i)
    {
        const label labelI = sLabels[i];

        // For boundary points, consider only other boundary points
        if ((! isThisInternalPoint) and (isInternalPoint[pointPoints[labelI]]))
        {
            continue;
        }

        if (counter == stride)
        {
            return labelI;
        }

        ++counter;
    }

    return UNDEF_LABEL;
}

// Help function for aspectRatioSmoothing to find out closest
// three points (connected by edges) to each point

int newSmoothSolverMotionSolver::findClosestPoints
(
    const fvMesh& mesh,
    const boolList& isInternalPoint,
    vectorList& closestPoints1,
    vectorList& closestPoints2,
    vectorList& closestPoints3,
    boolList& hasCommonCell,
    const labelListList& pointNeighPoints
)
{
    // Initialize with local information
    forAll (mesh.points(), pointI)
    {
        const vector cCoords = mesh.points()[pointI];
        const labelList& pointPoints = mesh.pointPoints()[pointI];
        const label nPoints = pointPoints.size();
        List<scalar> edgeLengths(nPoints, 0.0);

        // Calculate edge lengths and put into list
        forAll(pointPoints, neighPointI)
        {
            const label pointI = pointPoints[neighPointI];
            const double edgeLength = getPointDistance(mesh.points()[pointI], cCoords);
            edgeLengths[neighPointI] = edgeLength;
        }

        // Get order of sorted edgeLength
        // openfoam.com :
        // const labelList sLabels(Foam::sortedOrder(edgeLengths));
        //
        // openfoam.com and openfoam.org :
        labelList sLabels;
        Foam::sortedOrder(edgeLengths, sLabels);

        // Find correct closest point indices depending on whether
        // this is a boundary point or an internal point
        const label cLabel1 = findAppropriateClosestPointLabel(pointPoints, sLabels, pointI, isInternalPoint, 0);
        const label cLabel2 = findAppropriateClosestPointLabel(pointPoints, sLabels, pointI, isInternalPoint, 1);
        const label cLabel3 = findAppropriateClosestPointLabel(pointPoints, sLabels, pointI, isInternalPoint, 2);

        if (cLabel1 == UNDEF_LABEL)
        {
            FatalError << "Failed to find cLabel1 for pointI " << pointI << " at " << cCoords << endl << abort(FatalError);
        }

        if (cLabel2 == UNDEF_LABEL)
        {
            FatalError << "Failed to find cLabel2 for pointI " << pointI << " at " << cCoords << endl << abort(FatalError);
        }

        // cLabel3 may be undefined for processor interface cells at boundary
        //if (cLabel3 == UNDEF_LABEL)
        //{
        //    FatalError << "Failed to find cLabel3 for pointI " << pointI << " at " << cCoords << endl << abort(FatalError);
        //}

        // Save local values
        closestPoints1[pointI] = mesh.points()[pointPoints[cLabel1]] - cCoords;
        closestPoints2[pointI] = mesh.points()[pointPoints[cLabel2]] - cCoords;
        if (cLabel3 == UNDEF_LABEL)
        {
            closestPoints3[pointI] = UNDEF_VECTOR;
        }
        else
        {
            closestPoints3[pointI] = mesh.points()[pointPoints[cLabel3]] - cCoords;
        }

        // Check if closest two points share a cell
        if (findIndex(pointNeighPoints[pointPoints[cLabel1]], pointPoints[cLabel2]) >= 0)
            hasCommonCell[pointI] = true;
        else
            hasCommonCell[pointI] = false;
    }

    // For each position, deduce the globally closest points among processors

    // Position 1
    // ----------

    // Storage for sharing relative point locations among processors
    vectorField syncPoints(mesh.points().size(), Zero);

    forAll (mesh.points(), pointI)
    {
        syncPoints[pointI] = closestPoints1[pointI];
    }

    syncTools::syncPointList
    (
        mesh,
        syncPoints,
        minMagSqrEqOp<vector>(),
        UNDEF_VECTOR               // null value
    );

    forAll (mesh.points(), pointI)
    {
        if (isCloserPoint(syncPoints[pointI], closestPoints1[pointI]))
        {
            closestPoints3[pointI] = closestPoints2[pointI];
            closestPoints2[pointI] = closestPoints1[pointI];
            closestPoints1[pointI] = syncPoints[pointI];
            hasCommonCell[pointI] = false;
        }
    }

    // Position 2
    // ----------

    forAll (mesh.points(), pointI)
    {
        syncPoints[pointI] = closestPoints2[pointI];
    }

    syncTools::syncPointList
    (
        mesh,
        syncPoints,
        minMagSqrEqOp<vector>(),
        UNDEF_VECTOR               // null value
    );

    forAll (mesh.points(), pointI)
    {
        if (isCloserPoint(syncPoints[pointI], closestPoints2[pointI]))
        {
            closestPoints3[pointI] = closestPoints2[pointI];
            closestPoints2[pointI] = syncPoints[pointI];
            hasCommonCell[pointI] = false;
        }
    }

    // Position 3
    // ----------

    forAll (mesh.points(), pointI)
    {
        syncPoints[pointI] = closestPoints3[pointI];
    }

    syncTools::syncPointList
    (
        mesh,
        syncPoints,
        minMagSqrEqOp<vector>(),
        UNDEF_VECTOR               // null value
    );

    forAll (mesh.points(), pointI)
    {
        if (isCloserPoint(syncPoints[pointI], closestPoints3[pointI]))
        {
            closestPoints3[pointI] = syncPoints[pointI];
        }
    }

    // Synchronize cell sharing boolean
    syncTools::syncPointList
    (
         mesh,
         hasCommonCell,
         orEqOp<bool>(),
         false                      // null value
    );

    return 0;
}

// Help function for aspectRatioSmoothing to identify point locations
// of two short edge points which don't share cells.  Calculates also
// a blending fraction used for weighting the coordinates with
// centroidal smoothing coordinates. Blending fraction is set to zero
// if points share a cell, to prevent smoothing in that case.

double newSmoothSolverMotionSolver::calcARSmoothingRatio
(
    const vector& closestPoint1,
    const vector& closestPoint2,
    const vector& closestPoint3,
    const bool hasCommonCell,
    const bool isInternalPoint
)
{
    // Do nothing if the points are part of a common cell
    if (hasCommonCell)
    {
        return 0.0;
    }

    if ((closestPoint1 == ZERO_VECTOR) or (closestPoint2 == ZERO_VECTOR))
    {
        return 0.0;
    }

    const double lengthRatio1 = mag(closestPoint2) / mag(closestPoint1);
    const double lengthRatio2 = mag(closestPoint3) / mag(closestPoint2);

    // Aspect ratio treatment for internal points
    if (isInternalPoint)
    {
        // Minimum and maximum edge length ratios for detecting and
        // blending high aspect ratio.
        const double minRatio = 1.5;
        const double maxRatio = 3.0;

        // Blending fraction is applied if two closest points have similar
        // length, and if the third closest point is clearly farther away.
        if ((lengthRatio1 < minRatio) and (lengthRatio2 > minRatio))
        {
            const double frac = (lengthRatio2 - minRatio) / (maxRatio - minRatio);
            const double blendFrac = min(1.0, max(0.0, frac));
            return blendFrac;
        }
    }

    // Aspect ratio treatment for boundary points, to help
    // with edge folding when mesh is contracting against feature
    // edges
    else
    {
        const double minRatio = 1.0;
        const double maxRatio = 2.0;
        const double frac = (lengthRatio1 - minRatio) / (maxRatio - minRatio);
        const double blendFrac = min(1.0, max(0.0, frac));
        return blendFrac;
    }

    return 0.0;
}

// Special smoothing algorithm for high aspect ratio cells. Blends
// resulting coordinates with centroidal point coordinates.

Foam::tmp<Foam::pointField> newSmoothSolverMotionSolver::aspectRatioSmoothing
(
    const fvMesh& mesh,
    const boolList& isInternalPoint,
    const pointField& centroidalPoints,
    const labelListList& pointNeighPoints
)
{
    // Storage for surrounding point locations (relative to current
    // point) of closest edge points
    const label nPoints = mesh.nPoints();
    vectorList closestPoints1(nPoints, Zero);
    vectorList closestPoints2(nPoints, Zero);
    vectorList closestPoints3(nPoints, Zero);

    // Boolean to mark that the two shortest local edge points share a cell
    List<bool> hasCommonCell(nPoints, false);

    // New point locations storage
    tmp<pointField> tNewPoints(new pointField(nPoints, Zero));
    auto& newPoints = tNewPoints.ref();

    // Initialize with centroidal coordinates
    forAll(centroidalPoints, pointI)
    {
        newPoints[pointI] = centroidalPoints[pointI];
    }

    // Find closest points among processors
    findClosestPoints(mesh, isInternalPoint, closestPoints1, closestPoints2, closestPoints3, hasCommonCell, pointNeighPoints);

    // Calculate closest middle point coordinates and blend with centroidal coordinates
    forAll(mesh.points(), pointI)
    {
        const double blendFrac = calcARSmoothingRatio(closestPoints1[pointI], closestPoints2[pointI], closestPoints3[pointI], hasCommonCell[pointI], isInternalPoint[pointI]);

        if (blendFrac > 0.0)
        {
            const vector aCoords = mesh.points()[pointI] + (closestPoints1[pointI] + closestPoints2[pointI]) / 2.0;
            const vector newCoords = (1.0 - blendFrac) * centroidalPoints[pointI] + blendFrac * aCoords;
            newPoints[pointI] = newCoords;
        }
    }

    return tNewPoints;
}


// Quality control function which prohibits the decrease of edge
// length by freezing points to current locations, if edge length is
// below minEdgeLength and length is decreasing. This feature is used
// to limit the squishing or compression of cells near concave
// features.

int newSmoothSolverMotionSolver::restrictEdgeShortening
(
    const fvMesh& mesh,
    pointField& origPoints,
    const double minEdgeLength,
    const bool totalMinFreeze,
    boolList& isFrozenPoint
)
{
    forAll(origPoints, pointI)
    {
        if (isFrozenPoint[pointI])
            continue;

        const vector cCoords = mesh.points()[pointI];
        const vector nCoords = origPoints[pointI];

        // Calculate shortest edge length from current mesh point and
        // new mesh point
        double shortestCurrentEdgeLength = GREAT;
        double shortestNewEdgeLength = GREAT;
        forAll(mesh.pointPoints(pointI), pointPpI)
        {
            const label neighI = mesh.pointPoints(pointI)[pointPpI];
            const double testCurrentLength = getPointDistance(mesh.points()[neighI], cCoords);
            if (testCurrentLength < shortestCurrentEdgeLength)
                shortestCurrentEdgeLength = testCurrentLength;
            const double testNewLength = getPointDistance(mesh.points()[neighI], nCoords);
            if (testNewLength < shortestNewEdgeLength)
                shortestNewEdgeLength = testNewLength;
        }

        // Unconditional freeze below minEdgeLength (totalMinFreeze option)
        const double shortestLength = min(shortestNewEdgeLength, shortestCurrentEdgeLength);
        if (totalMinFreeze and (shortestLength < minEdgeLength))
        {
            isFrozenPoint[pointI] = true;
        }

        // Otherwise freeze only if edge length decreases and length
        // is below threshold value
        else if ((shortestNewEdgeLength < minEdgeLength) and
                 (shortestNewEdgeLength < shortestCurrentEdgeLength))
        {
            isFrozenPoint[pointI] = true;
        }

    }

    return 0;
}


// Calculate and return the maximum step length of the proposed new points

double newSmoothSolverMotionSolver::getProposedMaxStepLength
(
    const fvMesh& mesh,
    const pointField& newPoints
)
{
    double maxLength = 0.0;

    forAll(mesh.points(), pointI)
    {
        const double length = mag(vector(newPoints[pointI] - mesh.points()[pointI]));

        if (length > maxLength)
            maxLength = length;
    }

    const double proposedMaxLength = returnReduce(maxLength, maxOp<double>());

    return proposedMaxLength;
}


// Quality control function to constrain the length of a step jump to
// new coordinates by an absolute length value. This increases the
// stability of the smoothing process, in case target coordinates are
// far off.

int newSmoothSolverMotionSolver::constrainMaxStepLength
(
    const fvMesh& mesh,
    pointField& origPoints,
    const double maxStepLength,
    const double relStepFrac,
    const bool doGlobalScaling
)
{
    // Copy original points for temporary working point field
    tmp<pointField> tNewPoints(new pointField(mesh.nPoints(), Zero));
    auto& newPoints = tNewPoints.ref();

    forAll(newPoints, pointI)
        newPoints[pointI] = origPoints[pointI];

    // Global scaling factor is a means to scale down all point
    // movements with a common factor. It is meant to allow centroidal
    // smoothing to accelerate movement of points which need larger
    // step sizes than surrounding points for stable smoothing. If
    // global scaling factor is 1.0, then movement is scaled down by
    // the relative step fraction only.

    // Calculate global scaling factor from maximum proposed step size
    // and maximum allowed step size.
    const double proposedMaxLength = getProposedMaxStepLength(mesh, origPoints);

    double globalScale = 1.0;
    if (doGlobalScaling)
    {
        globalScale = min(1.0, maxStepLength / (proposedMaxLength * relStepFrac));
    }

    // Info << "Proposed maximum step length = " << proposedMaxLength << endl
    //      << "Maximum allowed step length = " << maxStepLength << endl
    //      << "Global scaling factor = " << globalScale << endl
    //      << endl << endl;

    forAll(origPoints, pointI)
    {
        // Scale down the length of the jump from current coordinates
        // towards new coordinates if jump would be too long
        const vector cCoords = mesh.points()[pointI];
        const vector stepDir = origPoints[pointI] - cCoords;

        // Max step length constraining if global scaling is disabled
        if (! doGlobalScaling)
        {
            if (mag(stepDir) > maxStepLength)
            {
                globalScale = maxStepLength / (mag(stepDir) * relStepFrac);
            }
            else
            {
                globalScale = 1.0;
            }
        }

        // Save the constrained point
        const vector nCoords = cCoords + relStepFrac * globalScale * stepDir;
        newPoints[pointI] = nCoords;
    }

    // Save new point coordinates to original point field
    forAll(origPoints, pointI)
    {
        origPoints[pointI] = newPoints[pointI];
    }

    return 0;
}


/////////////////////////////////////////////////////
// Help functions for restrictMinEdgeAngleDecrease //
/////////////////////////////////////////////////////

// Calculate and return the edge-edge angle (in radians,
// 0 < angle < pi) of two edges which share a common point
// at coordinate cCoords. The two edge end point coordinates are
// p1Coords and p2Coords.

double newSmoothSolverMotionSolver::edgeEdgeAngle
(
    const vector cCoords,
    const vector p1Coords,
    const vector p2Coords
)
{
    vector vec1 = (p1Coords - cCoords);
    vector vec2 = (p2Coords - cCoords);
    vec1 /= mag(vec1);
    vec2 /= mag(vec2);

    const double cosA = vec1 & vec2;

    // Ensure cos angle is in sane range before calling arc cos
    const double MAX = 0.99999;
    const double cosAlpha = std::max(-MAX, std::min(MAX, cosA));
    const double angle = std::acos(cosAlpha);

    return angle;
}


// Find the neighbour mesh point indices neighPI1 and neighPI2 for the
// point with index pointI, where all points are part of face index
// faceI.

int newSmoothSolverMotionSolver::getNeighbourPoints
(
    const fvMesh& mesh,
    const label centerPointI,
    const label faceI,
    label *neighPI1,
    label *neighPI2
)
{
    int i = 0;
    const face facePoints = mesh.faces()[faceI];
    const label nPoints = facePoints.size();

    forAll(facePoints, pointI)
    {
        if (facePoints[pointI] == centerPointI)
        {
            // Handle index wrapping at beginning and at end
            label prevI = pointI - 1;
            if (pointI == 0)
                prevI = nPoints - 1;

            label nextI = pointI + 1;
            if (pointI == nPoints - 1)
                nextI = 0;

            *neighPI1 = facePoints[prevI];
            *neighPI2 = facePoints[nextI];
            return 0;
        }
        i++;
    }

    FatalError << "Sanity broken, didn't find neighbour points for center point "
               << centerPointI << " and face " << faceI << endl
               << abort(FatalError);

    return 0;
}

// Calculate and store the minimum edge angles of point with index
// pointI for current mesh (minCAngleStorage) and a minimum for new
// mesh points (minNAngleStorage).

int newSmoothSolverMotionSolver::calc_min_edge_angles
(
    const fvMesh& mesh,
    const pointField& newPoints,
    const label pointI,
    double *minCAngleStorage,
    double *minNAngleStorage
)
{
    // Current and new minimum angles
    double minCAngle = DBL_MAX;
    double minNAngle = DBL_MAX;

    forAll(mesh.pointFaces()[pointI], pointFI)
    {
        // For both edges connected at this point, find the other point index
        label neighPI1 = 0;
        label neighPI2 = 0;
        const label faceI = mesh.pointFaces()[pointI][pointFI];
        getNeighbourPoints(mesh, pointI, faceI, &neighPI1, &neighPI2);

        // Current angle from current mesh
        const vector cp0 = mesh.points()[pointI];
        const vector cp1 = mesh.points()[neighPI1];
        const vector cp2 = mesh.points()[neighPI2];
        const double cAngle = edgeEdgeAngle(cp0, cp1, cp2);

        // New angle from current mesh
        const vector np0 = newPoints[pointI];
        const double nAngle0 = edgeEdgeAngle(np0, cp1, cp2);

        // New angle from new mesh points
        const vector np1 = newPoints[neighPI1];
        const vector np2 = newPoints[neighPI2];
        const double nAngle1 = edgeEdgeAngle(np0, np1, np2);

        // New angle from combination of current and new mesh points
        // Note: These might not be really necessary, but do it for now
        // for completeness sake.
        const double nAngle2 = edgeEdgeAngle(np0, cp1, np2);
        const double nAngle3 = edgeEdgeAngle(np0, np1, cp2);

        // Use smallest of new angles.
        const double nAngle = min(min(min(nAngle0, nAngle1), nAngle2), nAngle3);

        // Update minimum values
        if (cAngle < minCAngle)
            minCAngle = cAngle;
        if (nAngle < minNAngle)
            minNAngle = nAngle;
    }

    // Save minimum angles to storage
    *minCAngleStorage = minCAngle;
    *minNAngleStorage = minNAngle;

    return 0;
}

// Restrict decrease of smallest edge-edge angles when angle is below
// minAngle (in degrees). This is meant to avoid creation of
// self-intersections for concave features.

int newSmoothSolverMotionSolver::restrictMinEdgeAngleDecrease
(
    const fvMesh& mesh,
    pointField& origPoints,
    const double minEdgeAngleInDegrees,
    boolList& isFrozenPoint
)
{
    forAll(origPoints, pointI)
    {
        if (isFrozenPoint[pointI])
            continue;

        // Calculate current minimum angle and new minimum angle
        double minCAngle;
        double minNAngle;
        calc_min_edge_angles(mesh, origPoints, pointI, &minCAngle, &minNAngle);

        // If minimum angle is below threshold and would decrease,
        // then freeze the point. Note: This allows angle to
        // increase, so points are not frozen permanently.
        const double smallAngle = M_PI * minEdgeAngleInDegrees / 180.0;

        if ((minNAngle < smallAngle) and (minNAngle < minCAngle))
        {
            isFrozenPoint[pointI] = true;
        }
    }

    return 0;
}

///////////////////////////////////////////////////////
// Help functions for restrictFaceAngleDeterioration //
///////////////////////////////////////////////////////

// Map face angle information from edges to points

int newSmoothSolverMotionSolver::mapCurrentMinMaxFaceAnglesToPoints
(
    const fvMesh& mesh,
    const List<double>& minFaceAnglesForEdges,
    const List<double>& maxFaceAnglesForEdges,
    List<double>& minFaceAnglesForPoints,
    List<double>& maxFaceAnglesForPoints
)
{
    forAll(mesh.points(), pointI)
    {
        minFaceAnglesForPoints[pointI] = 2.0 * M_PI;
        maxFaceAnglesForPoints[pointI] = 0.0;
    }

    forAll(mesh.edges(), edgeI)
    {
        const edge e = mesh.edges()[edgeI];

        {
            const label pointI = e.start();
            if (minFaceAnglesForPoints[pointI] > minFaceAnglesForEdges[edgeI])
                minFaceAnglesForPoints[pointI] = minFaceAnglesForEdges[edgeI];
            if (maxFaceAnglesForPoints[pointI] < maxFaceAnglesForEdges[edgeI])
                maxFaceAnglesForPoints[pointI] = maxFaceAnglesForEdges[edgeI];
        }

        {
            const label pointI = e.end();
            if (minFaceAnglesForPoints[pointI] > minFaceAnglesForEdges[edgeI])
                minFaceAnglesForPoints[pointI] = minFaceAnglesForEdges[edgeI];
            if (maxFaceAnglesForPoints[pointI] < maxFaceAnglesForEdges[edgeI])
                maxFaceAnglesForPoints[pointI] = maxFaceAnglesForEdges[edgeI];
        }
    }

    return 0;
}

// Calculate sum of angles between two edge unit vectors and a center
// unit vector. Assumes unit length for all vectors.

double newSmoothSolverMotionSolver::calcEdgeCenterEdgeAngle
(
    const vector p0,
    const vector cC,
    const vector p1
)
{
    const double cosA0 = p0 & cC;
    const double cosA1 = cC & p1;

    // Ensure cos angle is in sane range before calling arc cos
    static const double MAX = 0.99999;
    const double cosAlpha0 = std::max(-MAX, std::min(MAX, cosA0));
    const double angle0 = std::acos(cosAlpha0);
    const double cosAlpha1 = std::max(-MAX, std::min(MAX, cosA1));
    const double angle1 = std::acos(cosAlpha1);

    return angle0 + angle1;
}

// Calculate minimum and maximum face angles from the collection of
// precalculated point data projected on the edge normal plane

int newSmoothSolverMotionSolver::calcMinMaxFinalProjectedAngle
(
    const label nCells,
    const vectorList& pVecs,
    const vectorList& cVecs,
    const labelList& faceIs,
    const labelList& f0Is,
    const labelList& f1Is,
    double& minFaceAngle,
    double& maxFaceAngle
)
{
    double minAngle = 2.0 * M_PI;
    double maxAngle = 0.0;

    for (label i = 0; i < nCells; ++i)
    {
        const label f0I = f0Is[i];
        const label f1I = f1Is[i];
        const vector p0 = pVecs[f0I];
        const vector p1 = pVecs[f1I];
        const vector cC = cVecs[i];
        const double angle = calcEdgeCenterEdgeAngle(p0, cC, p1);

        if (angle < minAngle)
            minAngle = angle;
        if (angle > maxAngle)
            maxAngle = angle;
    }

    minFaceAngle = minAngle;
    maxFaceAngle = maxAngle;

    return 0;
}

// Find the pair of faces in face index list faceIs which are part of
// the cell cellI, and store the face label pair.

int newSmoothSolverMotionSolver::findCellFacePair
(
    const fvMesh& mesh,
    const label cellI,
    const labelList& faceIs,
    const labelListList& cellFaces,
    label& f0I,
    label& f1I
)
{
    label face0I = UNDEF_LABEL;
    label face1I = UNDEF_LABEL;
    static const int nInternalFaces = mesh.owner().size();

    forAll (cellFaces[cellI], cellFaceI)
    {
        const label faceI = cellFaces[cellI][cellFaceI];
        const label faceIsI = findIndex(faceIs, faceI);

        if (faceIsI >= 0)
        {
            if (face0I == UNDEF_LABEL)
            {
                face0I = faceIsI;
            }
            else if (face1I == UNDEF_LABEL)
            {
                face1I = faceIsI;
            }
            else
            {
                FatalError << "Sanity broken, more than two edge faces belong to same cell" << endl << abort(FatalError);
            }
        }
    }

    // Sanity check
    if ((face0I == UNDEF_LABEL) or (face1I == UNDEF_LABEL) or (face0I == face1I) or (cellI == UNDEF_LABEL))
    {
        Info << "faceIs";
        for (int i = 0; i < faceIs.size(); ++i)
            Info << " " << faceIs[i];
        Info << " nInternalFaces " << nInternalFaces;
        Info << " cellI " << cellI << endl;
        Info << " cellFaces[cellI] " << cellFaces[cellI] << endl;
        FatalError << "Sanity broken, didn't find face pairs for cell "
                   << cellI << ". Face indices: " << face0I << " " << face1I << endl
                   << abort(FatalError);
    }

    // Save the face pair
    f0I = face0I;
    f1I = face1I;

    return 0;
}

// Calculate face center for face faceI as a weighted average of point coordinates.
// If point label pointI1 > -1, then that point is assumed to be at coords1.
// If point label pointI2 > -1, then that point is assumed to be at coords2.

vector newSmoothSolverMotionSolver::calcFaceCenter
(
    const fvMesh& mesh,
    const label faceI,
    const label pointI1,
    const vector coords1,
    const label pointI2,
    const vector coords2
)
{
    vector center = vector(0, 0, 0);

    forAll(mesh.faces()[faceI], facePointI)
    {
        const label pointI = mesh.faces()[faceI][facePointI];

        if ((pointI1 >= 0) and (pointI == pointI1))
            center += coords1;
        else if ((pointI2 >= 0) and (pointI == pointI2))
            center += coords2;
        else
            center += mesh.points()[pointI];
    }

    center /= double(mesh.faces()[faceI].size());

    return center;
}

// Calculate minimum and maximum face angles for a single edge, with
// optional move of pointI1 to coords1 and pointI2 to coords2.

int newSmoothSolverMotionSolver::calcMinMaxFaceAngleForEdge
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
)
{
    // Find all faces of this edge
    const labelList edgeFaces = mesh.edgeFaces(edgeI);
    const label nFaces = edgeFaces.size();

    const edge e = mesh.edges()[edgeI];

    // Edge start point coordinates
    const label e0I = e.start();
    vector e0 = mesh.points()[e0I];
    if ((pointI1 >= 0) and (e0I == pointI1))
        e0 = coords1;
    else if ((pointI2 >= 0) and (e0I == pointI2))
        e0 = coords2;

    // Edge end point coordinates
    const label e1I = e.end();
    vector e1 = mesh.points()[e1I];
    if ((pointI1 >= 0) and (e1I == pointI1))
        e1 = coords1;
    else if ((pointI2 >= 0) and (e1I == pointI2))
        e1 = coords2;

    // Edge center coordinates
    const vector cCoords = 0.5 * (e0 + e1);

    // Edge normal vector
    const vector eVec = (e1 - e0) / mag(e1 - e0);

    // Note: Edge center and edge normal vector defines the plane
    // where points are projected to, prior to angle calculation.

    // Calculate the projected face center vectors
    vectorList pVecs(nFaces, UNDEF_VECTOR);  // projected coordinates
    labelList faceIs(nFaces, UNDEF_LABEL);  // face indices of edge faces

    forAll(edgeFaces, edgeFaceI)
    {
        // Get face center coordinates
        const label faceI = edgeFaces[edgeFaceI];
        const vector fCoords = calcFaceCenter(mesh, faceI, pointI1, coords1, pointI2, coords2);

        // Project face center to the plane
        const vector cf = cCoords - fCoords;
        const double dotProd = cf & eVec;
        const vector pCoords = fCoords + dotProd * eVec;

        // Save projected face center coordinate vector
        const vector cp = (pCoords - cCoords) / mag(pCoords - cCoords);
        pVecs[edgeFaceI] = cp;

        // Save face index for pairing below
        faceIs[edgeFaceI] = faceI;
    }

    // Face-face angles are to be calculated for face pairs
    // belonging to all cells of the edge. There are always
    // exactly two faces for each cell, connected at the
    // edge. Search face indices for the face pairs.

    const labelList edgeCells = mesh.edgeCells(edgeI);
    const int nCells = edgeCells.size();
    labelList f0Is(nCells, UNDEF_LABEL);  // index of first face of the cell
    labelList f1Is(nCells, UNDEF_LABEL);  // index of second face of the cell
    vectorList cVecs(nCells, UNDEF_VECTOR);  // projected cell center coordinates

    // Calculate, project and save projected cell center coordinates
    forAll(edgeCells, i)
    {
        findCellFacePair(mesh, edgeCells[i], faceIs, cellFaces, f0Is[i], f1Is[i]);
        const label cellI = edgeCells[i];
        const vector cellCenter = mesh.C()[cellI];
        const vector cf = cCoords - cellCenter;
        const double dotProd = cf & eVec;
        const vector pCoords = cellCenter + dotProd * eVec;
        const vector cp = (pCoords - cCoords) / mag(pCoords - cCoords);
        cVecs[i] = cp;
    }

    // Finally calculate the minimum and maximum angles using the
    // projected point data and save to storage
    calcMinMaxFinalProjectedAngle(nCells, pVecs, cVecs, faceIs, f0Is, f1Is, minFaceAngle, maxFaceAngle);

    return 0;
}

// Wrapper for calcMinMaxEdgeAngle for current mesh points

int newSmoothSolverMotionSolver::calcMinMaxFaceAngleForCurrentMeshEdge
(
    const fvMesh& mesh,
    const labelListList& cellFaces,
    const label edgeI,
    double& minFaceAngle,
    double& maxFaceAngle
)
{
    static const vector dummy = vector(0, 0, 0);
    calcMinMaxFaceAngleForEdge(mesh, cellFaces, edgeI, minFaceAngle, maxFaceAngle, -1, dummy, -1, dummy);
    return 0;
}

// Calculate and store the minimum face angles for all mesh
// edges

int newSmoothSolverMotionSolver::calcCurrentMinMaxFaceAnglesForEdges
(
    const fvMesh& mesh,
    const labelListList& cellFaces,
    List<double>& minFaceAnglesForEdges,
    List<double>& maxFaceAnglesForEdges
)
{
    forAll(mesh.edges(), edgeI)
    {
        double minAngle;
        double maxAngle;
        calcMinMaxFaceAngleForCurrentMeshEdge(mesh, cellFaces, edgeI, minAngle, maxAngle);
        minFaceAnglesForEdges[edgeI] = minAngle;
        maxFaceAnglesForEdges[edgeI] = maxAngle;
    }

    return 0;
}

// Calculate and store minimum and maximum face angle for argument
// first point pointI1, with optional movement of pointI2 to given
// coordinates.

int newSmoothSolverMotionSolver::calcMinMaxFaceAngleForPoint
(
    const fvMesh& mesh,
    const labelListList& cellFaces,
    const int pointI1,
    const vector coords1,
    const int pointI2,
    const vector coords2,
    double& minFaceAngle,
    double& maxFaceAngle
)
{
    // Initialize
    minFaceAngle = 2.0 * M_PI;
    maxFaceAngle = 0.0;

    // Loop over all point edges to find min and max angle

    forAll(mesh.pointEdges()[pointI1], pointEdgeI)
    {
        double minAngle;
        double maxAngle;
        const label edgeI = mesh.pointEdges()[pointI1][pointEdgeI];
        calcMinMaxFaceAngleForEdge(mesh, cellFaces, edgeI, minAngle, maxAngle, pointI1, coords1, pointI2, coords2);

        if (minFaceAngle > minAngle)
            minFaceAngle = minAngle;
        if (maxFaceAngle < maxAngle)
            maxFaceAngle = maxAngle;
    }

    return 0;
}


// Quality control function to restrict the decrease of smallest
// face-face angles when angle is below minAngle (in degrees). This is
// meant to avoid creation of self-intersections. The face-face angle
// can be thought of a measure which quantifies squishing or folding
// of cell faces. In addition to prohibiting the deterioration of
// face-face angles of each point, it is also necessary to freeze
// movement of neighbouring points if their movement would deteriorate
// the angles at the current point.

int newSmoothSolverMotionSolver::restrictFaceAngleDeterioration
(
    const fvMesh& mesh,
    const labelListList& cellFaces,
    pointField& origPoints,
    const double minFaceAngleInDegrees,
    const double maxFaceAngleInDegrees,
    boolList& isFrozenPoint
)
{
    // Calculate minimum and maximum face angles for all edges from
    // current mesh
    static const size_t nEdges = size_t(mesh.nEdges());
    List<double> currentMinAnglesForEdges(nEdges, GREAT);
    List<double> currentMaxAnglesForEdges(nEdges, GREAT);
    calcCurrentMinMaxFaceAnglesForEdges(mesh, cellFaces, currentMinAnglesForEdges, currentMaxAnglesForEdges);

    // Map face angle information from edges to points
    static const size_t nPoints = size_t(mesh.nPoints());
    List<double> currentMinAnglesForPoints(nPoints, GREAT);
    List<double> currentMaxAnglesForPoints(nPoints, GREAT);
    mapCurrentMinMaxFaceAnglesToPoints(mesh, currentMinAnglesForEdges, currentMaxAnglesForEdges, currentMinAnglesForPoints, currentMaxAnglesForPoints);

    // Use a stack to walk through all points. If a point is frozen by
    // it's neighbour then the point must be processed again to allow
    // recursive neighbour freezing. Use of stack is needed, since
    // list of items to be processed is modified on the run.
    std::stack<label> pointStack;

    // Initialize stack with all points. Angle calculation must be
    // made for all points, also boundary points, in order to be able
    // to stop angle deterioration at boundary by neighbour point
    // movement
    forAll(origPoints, pointI)
        pointStack.push(pointI);

    while (! pointStack.empty())
    {
        // Get and remove a point label from stack.
        const label pointI = pointStack.top();
        pointStack.pop();

        // 1. Check nothing for points whose face angles are in good range

        const double smallAngle = M_PI * minFaceAngleInDegrees / 180.0;
        const double largeAngle = M_PI * maxFaceAngleInDegrees / 180.0;

        if ((currentMinAnglesForPoints[pointI] > smallAngle) and
            (currentMaxAnglesForPoints[pointI] < largeAngle))
            continue;

        const vector cCoords = mesh.points()[pointI];
        vector nCoords = origPoints[pointI];

        // If this point is already frozen, set coordinates to current
        // mesh coordinates
        if (isFrozenPoint[pointI])
            nCoords = cCoords;

        // 2. Calculate new face angles for this point, if it's not
        // frozen. This point is hypothetically moved to new
        // coordinates, while surrounding points are kept at current
        // locations. Freeze this point if angle change with the move
        // of this point deteriorates the angles.

        if (nCoords != cCoords)
        {
            double newMinFaceAngle;
            double newMaxFaceAngle;
            calcMinMaxFaceAngleForPoint(mesh, cellFaces, pointI, nCoords, -1, nCoords, newMinFaceAngle, newMaxFaceAngle);

            if (((newMinFaceAngle < smallAngle) and
                (newMinFaceAngle < currentMinAnglesForPoints[pointI])) or
                ((newMaxFaceAngle > largeAngle) and
                (newMaxFaceAngle > currentMaxAnglesForPoints[pointI])))
            {
                // Freeze this point (self freeze)
                nCoords = cCoords;
                isFrozenPoint[pointI] = true;
            }
        }

        // 3. Calculate the effect from all neighbouring point movements
        // to face angles at this point. Freeze the neighbour point if
        // angle becomes worse than it currently is.

        forAll(mesh.pointPoints()[pointI], pointNI)
        {
            const label neighPointI = mesh.pointPoints()[pointI][pointNI];
            const vector neighCoords = origPoints[neighPointI];

            // Skip this neighbour point if it's not moving
            if (isFrozenPoint[neighPointI])
                continue;
            if (neighCoords == mesh.points()[neighPointI])
                continue;

            double newMinFaceAngle;
            double newMaxFaceAngle;
            calcMinMaxFaceAngleForPoint(mesh, cellFaces, pointI, nCoords, neighPointI, neighCoords, newMinFaceAngle, newMaxFaceAngle);

            if (((newMinFaceAngle < smallAngle) and
                (newMinFaceAngle < currentMinAnglesForPoints[pointI])) or
                ((newMaxFaceAngle > largeAngle) and
                (newMaxFaceAngle > currentMaxAnglesForPoints[pointI])))
            {
                // Freeze the neighbour point (neighbour freeze)
                isFrozenPoint[neighPointI] = true;

                // Add neighbour point index to stack list, as it
                // needs to be (re)checked after neighbour freezing
                pointStack.push(neighPointI);
            }
        }
    }

    return 0;
}


// Help function to get patch numbers from a given command line option

labelList newSmoothSolverMotionSolver::getPatchIdsForOption
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& optionName,
    const bool defaultAll
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    labelHashSet patchSet(patches.size());

    if (dict.found(optionName))
    {
        const wordReList patterns(dict.lookup(optionName));
        patchSet = patches.patchSet(patterns);
    }
    else if (defaultAll)
    {
        forAll(patches, patchi)
        {
            patchSet.insert(patchi);
        }
    }

    labelList patchIds(patchSet.size());
    label i = 0;
    forAllConstIter(labelHashSet, patchSet, iter)
    {
        patchIds[i++] = iter.key();
    }
    return patchIds;
}

// Calculate and return mesh statistics: the minimum and maximum edge
// length of the current mesh, and the total perimeter of the mesh,
// calculated as a sum of bounding box side lengths

int newSmoothSolverMotionSolver::getMeshStats
(
    const fvMesh& mesh,
    double& minEdgeLength,
    double& maxEdgeLength,
    double& meshPerimeter
)
{
    double minLength = VGREAT;
    double maxLength = 0.0;
    double bbMinX = VGREAT;
    double bbMaxX = -VGREAT;
    double bbMinY = VGREAT;
    double bbMaxY = -VGREAT;
    double bbMinZ = VGREAT;
    double bbMaxZ = -VGREAT;

    forAll(mesh.edges(), edgeI)
    {
        const edge e = mesh.edges()[edgeI];

        const label startI = e.start();
        const vector startCoords = mesh.points()[startI];
        const label endI = e.end();
        const vector endCoords = mesh.points()[endI];

        const double length = mag(vector(endCoords - startCoords));

        if (length < minLength)
            minLength = length;
        if (length > maxLength)
            maxLength = length;

        // Bounding box update
        if (startCoords[0] < bbMinX) { bbMinX = startCoords[0]; }
        if (startCoords[1] < bbMinY) { bbMinY = startCoords[1]; }
        if (startCoords[2] < bbMinZ) { bbMinZ = startCoords[2]; }
        if (startCoords[0] > bbMaxX) { bbMaxX = startCoords[0]; }
        if (startCoords[1] > bbMaxY) { bbMaxY = startCoords[1]; }
        if (startCoords[2] > bbMaxZ) { bbMaxZ = startCoords[2]; }

        if (endCoords[0] < bbMinX) { bbMinX = endCoords[0]; }
        if (endCoords[1] < bbMinY) { bbMinY = endCoords[1]; }
        if (endCoords[2] < bbMinZ) { bbMinZ = endCoords[2]; }
        if (endCoords[0] > bbMaxX) { bbMaxX = endCoords[0]; }
        if (endCoords[1] > bbMaxY) { bbMaxY = endCoords[1]; }
        if (endCoords[2] > bbMaxZ) { bbMaxZ = endCoords[2]; }
    }

    const double meshMinLength = returnReduce(minLength, minOp<double>());
    const double meshMaxLength = returnReduce(maxLength, maxOp<double>());
    const double meshBbMinX = returnReduce(bbMinX, minOp<double>());
    const double meshBbMaxX = returnReduce(bbMaxX, maxOp<double>());
    const double meshBbMinY = returnReduce(bbMinY, minOp<double>());
    const double meshBbMaxY = returnReduce(bbMaxY, maxOp<double>());
    const double meshBbMinZ = returnReduce(bbMinZ, minOp<double>());
    const double meshBbMaxZ = returnReduce(bbMaxZ, maxOp<double>());

    minEdgeLength = meshMinLength;
    maxEdgeLength = meshMaxLength;
    meshPerimeter = meshBbMaxX - meshBbMinX + meshBbMaxY - meshBbMinY + meshBbMaxZ + meshBbMinZ;

    return 0;
}

// Calculate the maximum relative change in the point positions, to be
// used as a measure of residual of smoothing

double newSmoothSolverMotionSolver::calculateResidual
(
    const fvMesh& mesh,
    const pointField& newPoints,
    const boolList& isInternalPoint,
    const double maxStepLength
)
{
    scalar maxStep = ZERO;

    forAll(isInternalPoint, pointI)
    {
        const scalar distance =
            mag(newPoints[pointI] - mesh.points()[pointI]) / maxStepLength;

        if (distance > maxStep)
        {
            maxStep = distance;
        }
    }

    const scalar maxMaxStep = returnReduce(maxStep, maxOp<scalar>());

    return maxMaxStep;
}


// Help function to generate lists of cell faces for all cells

int newSmoothSolverMotionSolver::generateCellFaces
(
    const fvMesh& mesh,
    labelListList& cellFaces
)
{
    // Traverse internal faces
    forAll (mesh.owner(), faceI)
    {
        const label cellI = mesh.owner()[faceI];
        cellFaces[cellI].append(faceI);
    }

    forAll (mesh.neighbour(), faceI)
    {
        const label cellI = mesh.neighbour()[faceI];
        cellFaces[cellI].append(faceI);

        // Sanity check, can be removed after confirming this never triggers
        if (cellI < 0)
        {
            FatalError << "Found undefined neighbour cellI for faceI " << faceI << endl << abort(FatalError);
        }
    }

    // Traverse boundary faces
    forAll (mesh.boundary(), patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        const label startI = mesh.boundary()[patchI].start();
        // const label endI = startI + mesh.boundary()[patchI].Cf().size();

        forAll (mesh.boundary()[patchI].Cf(), patchFaceI)
        {
            const label faceI = startI + patchFaceI;
            const label& cellI = pp.faceCells()[patchFaceI];
            if (cellI >= 0)
            {
                cellFaces[cellI].append(faceI);
            }
        }
    }

    return 0;
}

// Help function to check does a file exists on the system

bool newSmoothSolverMotionSolver::fileExists
(
    const std::string& fileName
)
{
    std::ifstream file(fileName);
    return file.good();
}



/*End list of sample function definitions*/

} // End namespace Foam

