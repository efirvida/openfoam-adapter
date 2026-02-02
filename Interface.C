#include "Interface.H"
#include "Utilities.H"
#include "faceTriangulation.H"
#include "cellSet.H"


using namespace Foam;

preciceAdapter::Interface::Interface(
    precice::Participant& precice,
    const fvMesh& mesh,
    std::string meshName,
    std::string locationsType,
    std::vector<std::string> patchNames,
    std::vector<std::string> cellSetNames,
    bool meshConnectivity,
    bool restartFromDeformed,
    const std::string& namePointDisplacement,
    const std::string& nameCellDisplacement)
: precice_(precice),
  meshName_(meshName),
  patchNames_(patchNames),
  cellSetNames_(cellSetNames),
  meshConnectivity_(meshConnectivity),
  restartFromDeformed_(restartFromDeformed)
{
    dim_ = precice_.getMeshDimensions(meshName);

    if (dim_ == 2 && meshConnectivity_ == true)
    {
        DEBUG(adapterInfo("meshConnectivity is currently only supported for 3D cases. \n"
                          "You might set up a 3D case and restrict the 3rd dimension by z-dead = true. \n"
                          "Have a look in the adapter documentation for detailed information.",
                          "warning"));
    }

    if (locationsType == "faceCenters" || locationsType == "faceCentres")
    {
        locationType_ = LocationType::faceCenters;
    }
    else if (locationsType == "faceNodes")
    {
        locationType_ = LocationType::faceNodes;
    }
    else if (locationsType == "volumeCenters" || locationsType == "volumeCentres")
    {
        locationType_ = LocationType::volumeCenters;
    }
    else
    {
        adapterInfo("Interface points location type \""
                    "locations = "
                        + locationsType + "\" is invalid.",
                    "error-deferred");
    }


    // For every patch that participates in the coupling
    for (uint j = 0; j < patchNames.size(); j++)
    {
        // Get the patchID
        int patchID = mesh.boundaryMesh().findPatchID(patchNames.at(j));

        // Throw an error if the patch was not found
        if (patchID == -1)
        {
            adapterInfo("Patch \""
                            + patchNames.at(j) + "\" does not exist and therefore cannot be used as a coupling interface for mesh \""
                            + meshName + "\". Check the system/preciceDict.",
                        "error");
        }

        // Add the patch in the list
        patchIDs_.push_back(patchID);
    }

    // Check if this is a global data interface (no patches)
    isGlobalDataInterface_ = patchNames.empty();

    // Configure the mesh (set the data locations)
    configureMesh(mesh, namePointDisplacement, nameCellDisplacement);
}

void preciceAdapter::Interface::configureMesh(const fvMesh& mesh, const std::string& namePointDisplacement, const std::string& nameCellDisplacement)
{
    // Handle global data interface (single vertex at origin)
    if (isGlobalDataInterface_)
    {
        DEBUG(adapterInfo("Configuring global data interface with single vertex at origin"));

        numDataLocations_ = 1;
        vertexIDs_.resize(1);

        // Single vertex at origin (0, 0, 0)
        std::vector<double> vertices(dim_, 0.0);

        // Pass the single vertex to preCICE
        precice_.setMeshVertices(meshName_, vertices, vertexIDs_);

        DEBUG(adapterInfo("Global data interface configured: 1 vertex at (0,0,0)"));
        return;
    }

    // The way we configure the mesh differs between meshes based on face centers
    // and meshes based on face nodes.
    // TODO: Reduce code duplication. In the meantime, take care to update
    // all the branches.

    if (locationType_ == LocationType::faceCenters)
    {
        // Count the data locations for all the patches
        for (uint j = 0; j < patchIDs_.size(); j++)
        {
            numDataLocations_ +=
                mesh.boundaryMesh()[patchIDs_.at(j)].faceCentres().size();
        }
        DEBUG(adapterInfo("Number of face centres: " + std::to_string(numDataLocations_)));

        // In case we want to perform the reset later on, look-up the corresponding data field name
        Foam::volVectorField const* cellDisplacement = nullptr;
        if (mesh.foundObject<volVectorField>(nameCellDisplacement))
            cellDisplacement =
                &mesh.lookupObject<volVectorField>(nameCellDisplacement);

        // Array of the mesh vertices.
        // One mesh is used for all the patches and each vertex has 3D coordinates.
        std::vector<double> vertices(dim_ * numDataLocations_);

        // Array of the indices of the mesh vertices.
        // Each vertex has one index, but three coordinates.
        vertexIDs_.resize(numDataLocations_);

        // Initialize the index of the vertices array
        int verticesIndex = 0;

        // Get the locations of the mesh vertices (here: face centers)
        // for all the patches
        for (uint j = 0; j < patchIDs_.size(); j++)
        {
            // Get the face centers of the current patch
            vectorField faceCenters =
                mesh.boundaryMesh()[patchIDs_.at(j)].faceCentres();

            // Move the interface according to the current values of the cellDisplacement field,
            // to account for any displacements accumulated before restarting the simulation.
            // This is information that OpenFOAM reads from its result/restart files.
            // If the simulation is not restarted, the displacement should be zero and this line should have no effect.
            if (cellDisplacement != nullptr && !restartFromDeformed_)
                faceCenters -= cellDisplacement->boundaryField()[patchIDs_.at(j)];

            // Assign the (x,y,z) locations to the vertices
            for (int i = 0; i < faceCenters.size(); i++)
                for (unsigned int d = 0; d < dim_; ++d)
                    vertices[verticesIndex++] = faceCenters[i][d];

            // Check if we are in the right layer in case of preCICE dimension 2
            // If there is at least one node with a different z-coordinate, then the (2D) geometry is not on the xy-plane, as required.
            if (dim_ == 2)
            {
                const pointField faceNodes =
                    mesh.boundaryMesh()[patchIDs_.at(j)].localPoints();
                const auto faceNodesSize = faceNodes.size();
                // Allocate memory for z-coordinates
                std::array<double, 2> z_location({0, 0});
                constexpr unsigned int z_axis = 2;

                // Find out about the existing planes
                // Store z-coordinate of the first layer
                if (faceNodesSize > 0)
                {
                    z_location[0] = faceNodes[0][z_axis];
                }
                // Go through the remaining points until we find the second z-coordinate
                // and store it (there are only two allowed in case we are in the xy-layer)
                for (int i = 0; i < faceNodesSize; i++)
                {
                    if (z_location[0] == faceNodes[i][z_axis])
                    {
                        continue;
                    }
                    else
                    {
                        z_location[1] = faceNodes[i][z_axis];
                        break;
                    }
                }

                // Check if the z-coordinates of all nodes match the z-coordinates we have collected above
                for (int i = 0; i < faceNodesSize; i++)
                {
                    if (z_location[0] == faceNodes[i][z_axis] || z_location[1] == faceNodes[i][z_axis])
                    {
                        continue;
                    }
                    else
                    {
                        adapterInfo("It seems like you are using preCICE in 2D and your geometry is not located int the xy-plane. "
                                    "The OpenFOAM adapter implementation supports preCICE 2D cases only with the z-axis as out-of-plane direction."
                                    "Please rotate your geometry so that the geometry is located in the xy-plane."
                                    "If you are running a 2D axisymmetric case just ignore this.",
                                    "warning");
                    }
                }
            }
        }

        // === GATHER-TO-MASTER: Collect all vertices at rank 0 ===

        // Step 1: Gather counts from all ranks
        int localVertexCount = numDataLocations_;
        int localScalarCount = numDataLocations_ * dim_;

        gatherCounts_.setSize(Pstream::nProcs());
        gatherCounts_ = 0; // Initialize all to 0
        gatherCounts_[Pstream::myProcNo()] = localScalarCount;
        Pstream::gatherList(gatherCounts_);
        Pstream::broadcast(gatherCounts_);

        // Step 2: Calculate global total and displacements
        globalNumDataLocations_ = 0;
        gatherDisplacements_.setSize(Pstream::nProcs());
        gatherDisplacements_ = 0;
        for (int i = 0; i < Pstream::nProcs(); i++)
        {
            gatherDisplacements_[i] = globalNumDataLocations_ * dim_;
            globalNumDataLocations_ += gatherCounts_[i] / dim_;
        }

        DEBUG(Pout << "Adapter [Procid " << Pstream::myProcNo() << "]: Local vertices (faceCenters): "
                   << localVertexCount << ", Global: " << globalNumDataLocations_ << endl);

        // Step 3: Gather all vertices to rank 0
        if (Pstream::master())
        {
            globalDataBuffer_.resize(globalNumDataLocations_ * dim_);
            globalVertexIDs_.resize(globalNumDataLocations_);
        }

        List<double> localVerticesList(vertices.size());
        forAll(localVerticesList, i)
        {
            localVerticesList[i] = vertices[i];
        }

        List<List<double>> allVerticesLists(Pstream::nProcs());
        allVerticesLists[Pstream::myProcNo()] = localVerticesList;
        Pstream::gatherList(allVerticesLists);

        // Step 4: All ranks call setMeshVertices (master with data, others with empty)
        if (Pstream::master())
        {
            int globalIndex = 0;
            for (int rank = 0; rank < Pstream::nProcs(); rank++)
            {
                const List<double>& rankVertices = allVerticesLists[rank];
                forAll(rankVertices, i)
                {
                    globalDataBuffer_[globalIndex++] = rankVertices[i];
                }
            }

            DEBUG(Pout << "Adapter [Master]: Registering " << globalNumDataLocations_
                       << " global vertices (faceCenters) with preCICE" << endl);

            precice_.setMeshVertices(meshName_, globalDataBuffer_, globalVertexIDs_);
        }
        else
        {
            std::vector<double> emptyVertices;
            std::vector<int> dummyIDs;
            precice_.setMeshVertices(meshName_, emptyVertices, dummyIDs);
        }
    }
    else if (locationType_ == LocationType::faceNodes)
    {
        // Count the data locations for all the patches
        for (uint j = 0; j < patchIDs_.size(); j++)
        {
            numDataLocations_ +=
                mesh.boundaryMesh()[patchIDs_.at(j)].localPoints().size();
        }
        DEBUG(adapterInfo("Number of face nodes: " + std::to_string(numDataLocations_)));

        // In case we want to perform the reset later on, look-up the corresponding data field name
        Foam::pointVectorField const* pointDisplacement = nullptr;
        if (mesh.foundObject<pointVectorField>(namePointDisplacement))
            pointDisplacement =
                &mesh.lookupObject<pointVectorField>(namePointDisplacement);

        // Array of the mesh vertices.
        // One mesh is used for all the patches and each vertex has 3D coordinates.
        std::vector<double> vertices(dim_ * numDataLocations_);

        // Array of the indices of the mesh vertices.
        // Each vertex has one index, but three coordinates.
        vertexIDs_.resize(numDataLocations_);

        // Initialize the index of the vertices array
        int verticesIndex = 0;

        // Get the locations of the mesh vertices (here: face nodes)
        // for all the patches
        for (uint j = 0; j < patchIDs_.size(); j++)
        {
            // Get the face nodes of the current patch
            // TODO: Check if this is correct.
            // TODO: Check if this behaves correctly in parallel.
            // TODO: Check if this behaves correctly with multiple, connected patches.
            // TODO: Maybe this should be a pointVectorField?
            pointField faceNodes =
                mesh.boundaryMesh()[patchIDs_.at(j)].localPoints();

            // Similar to the cell displacement above:
            // Move the interface according to the current values of the cellDisplacement field,
            // to account for any displacements accumulated before restarting the simulation.
            // This is information that OpenFOAM reads from its result/restart files.
            // If the simulation is not restarted, the displacement should be zero and this line should have no effect.
            if (pointDisplacement != nullptr && !restartFromDeformed_)
            {
                const vectorField& resetField = refCast<const vectorField>(
                    pointDisplacement->boundaryField()[patchIDs_.at(j)]);
                faceNodes -= resetField;
            }

            // Assign the (x,y,z) locations to the vertices
            // TODO: Ensure consistent order when writing/reading
            for (int i = 0; i < faceNodes.size(); i++)
            {
                for (unsigned int d = 0; d < dim_; ++d)
                {
                    vertices[verticesIndex++] = faceNodes[i][d];
                }
            }
        }

        // === GATHER-TO-MASTER: Collect all vertices at rank 0 ===

        // Step 1: Gather counts from all ranks
        int localVertexCount = numDataLocations_;        // Number of vertices (not scalar values)
        int localScalarCount = numDataLocations_ * dim_; // Number of scalar values

        // Use Pstream to gather counts
        gatherCounts_.setSize(Pstream::nProcs());
        gatherCounts_ = 0;
        gatherCounts_[Pstream::myProcNo()] = localScalarCount;
        Pstream::gatherList(gatherCounts_);
        Pstream::broadcast(gatherCounts_); // All ranks need this info for later scatter

        // Step 2: Calculate global total and displacements
        globalNumDataLocations_ = 0;
        gatherDisplacements_.setSize(Pstream::nProcs());
        gatherDisplacements_ = 0;
        for (int i = 0; i < Pstream::nProcs(); i++)
        {
            gatherDisplacements_[i] = globalNumDataLocations_ * dim_;
            globalNumDataLocations_ += gatherCounts_[i] / dim_; // Convert back to vertex count
        }

        DEBUG(Pout << "Adapter [Procid " << Pstream::myProcNo() << "]: Local vertices: " << localVertexCount
                   << ", Global vertices: " << globalNumDataLocations_ << endl);

        // Step 3: Gather all vertices to rank 0
        if (Pstream::master())
        {
            // Master allocates global buffers
            globalDataBuffer_.resize(globalNumDataLocations_ * dim_);
            globalVertexIDs_.resize(globalNumDataLocations_);
        }

        // Gather vertices using OpenFOAM's Pstream (wraps MPI)
        // Convert to List for Pstream compatibility
        List<double> localVerticesList(vertices.size());
        forAll(localVerticesList, i)
        {
            localVerticesList[i] = vertices[i];
        }

        List<List<double>> allVerticesLists(Pstream::nProcs());
        allVerticesLists[Pstream::myProcNo()] = localVerticesList;
        Pstream::gatherList(allVerticesLists);

        // Step 4: All ranks call setMeshVertices (master with data, others with empty)
        if (Pstream::master())
        {
            int globalIndex = 0;
            for (int rank = 0; rank < Pstream::nProcs(); rank++)
            {
                const List<double>& rankVertices = allVerticesLists[rank];
                forAll(rankVertices, i)
                {
                    globalDataBuffer_[globalIndex++] = rankVertices[i];
                }
            }

            DEBUG(Pout << "Adapter [Master]: Registering " << globalNumDataLocations_
                       << " global vertices with preCICE" << endl);

            // Only master registers mesh with preCICE
            precice_.setMeshVertices(meshName_, globalDataBuffer_, globalVertexIDs_);
        }
        else
        {
            std::vector<double> emptyVertices;
            std::vector<int> dummyIDs;
            precice_.setMeshVertices(meshName_, emptyVertices, dummyIDs);
        }

        // Broadcast globalVertexIDs_ is not needed - only master uses them for preCICE calls
        // But we keep local vertexIDs_ for reference (may be useful for debugging)

        if (meshConnectivity_)
        {
            // === GATHER-TO-MASTER: Collect all triangles at rank 0 ===

            // Each rank calculates its triangles as indices into its LOCAL vertices.
            // We then offset these indices to point to the correct position in the GLOBAL master buffer.
            int localVertexOffset = 0;
            int rankVertexOffset = gatherDisplacements_[Pstream::myProcNo()] / dim_;

            std::vector<int> localTriVertIDs;
            const int triaPerQuad = 2;
            const int nodesPerTria = 3;

            for (uint j = 0; j < patchIDs_.size(); j++)
            {
                const List<face> faceField = mesh.boundaryMesh()[patchIDs_.at(j)].localFaces();
                Field<point> pointCoords = mesh.boundaryMesh()[patchIDs_.at(j)].localPoints();
                const int numPatchPoints = pointCoords.size();

                // Triangulate faces
                forAll(faceField, facei)
                {
                    const face& faceQuad = faceField[facei];
                    faceTriangulation faceTri(pointCoords, faceQuad, false);

                    for (uint triIndex = 0; triIndex < triaPerQuad; triIndex++)
                    {
                        for (uint nodeIndex = 0; nodeIndex < nodesPerTria; nodeIndex++)
                        {
                            int patchPointIndex = faceTri[triIndex][nodeIndex];
                            // The global index in the master's vertex buffer is:
                            // Rank's Starting Index + Patch's Offset in Rank + Point's Index in Patch
                            localTriVertIDs.push_back(rankVertexOffset + localVertexOffset + patchPointIndex);
                        }
                    }
                }
                localVertexOffset += numPatchPoints;
            }

            // Gather all triangle lists to master
            List<int> localTriList(localTriVertIDs.size());
            forAll(localTriList, i)
            {
                localTriList[i] = localTriVertIDs[i];
            }

            List<List<int>> allTriLists(Pstream::nProcs());
            allTriLists[Pstream::myProcNo()] = localTriList;
            Pstream::gatherList(allTriLists);

            if (Pstream::master())
            {
                std::vector<int> globalTriVertIDs;
                for (int rank = 0; rank < Pstream::nProcs(); rank++)
                {
                    const List<int>& rankTris = allTriLists[rank];
                    forAll(rankTris, i)
                    {
                        // Indices already offset, so we only need the vertex ID from the master's registration
                        globalTriVertIDs.push_back(globalVertexIDs_[rankTris[i]]);
                    }
                }

                DEBUG(Pout << "Adapter [Master]: Registering " << globalTriVertIDs.size() / 3
                           << " global triangles with preCICE" << endl);
                precice_.setMeshTriangles(meshName_, globalTriVertIDs);
            }
        }
    }
    else if (locationType_ == LocationType::volumeCenters)
    {
        // The volume coupling implementation considers the mesh points in the volume and
        // on the boundary patches in order to take the boundary conditions into account

        // Get the cell labels of the overlapping region
        std::vector<labelList> overlapCells;

        if (!cellSetNames_.empty())
        {
            // For every cellSet that participates in the coupling
            for (uint j = 0; j < cellSetNames_.size(); j++)
            {
                // Create a cell set
                cellSet overlapRegion(mesh, cellSetNames_[j]);

                // Add the cells IDs to the vector and count how many overlap cells the interface has
                overlapCells.push_back(overlapRegion.toc());
                numDataLocations_ += overlapCells[j].size();
            }
        }
        else
        {
            numDataLocations_ = mesh.C().size();
        }

        // Count the data locations for all the patches
        // and add those to the previously determined number of mesh points in the volume
        for (uint j = 0; j < patchIDs_.size(); j++)
        {
            numDataLocations_ +=
                mesh.boundaryMesh()[patchIDs_.at(j)].faceCentres().size();
        }
        DEBUG(adapterInfo("Number of coupling volumes: " + std::to_string(numDataLocations_)));

        // Array of the mesh vertices.
        // One mesh is used for all the patches and each vertex has 3D coordinates.
        std::vector<double> vertices(dim_ * numDataLocations_);

        // Array of the indices of the mesh vertices.
        // Each vertex has one index, but three coordinates.
        vertexIDs_.resize(numDataLocations_);

        // Initialize the index of the vertices array
        int verticesIndex = 0;

        if (!cellSetNames_.empty())
        {
            // for all the overlapping cells (cellSets)
            for (uint j = 0; j < cellSetNames_.size(); j++)
            {
                // Get the cell centres of the current cellSet.
                const labelList& cells = overlapCells.at(j);

                // Get the coordinates of the cells of the current cellSet.
                for (int i = 0; i < cells.size(); i++)
                {
                    vertices[verticesIndex++] = mesh.C().internalField()[cells[i]].x();
                    vertices[verticesIndex++] = mesh.C().internalField()[cells[i]].y();
                    if (dim_ == 3)
                    {
                        vertices[verticesIndex++] = mesh.C().internalField()[cells[i]].z();
                    }
                }
            }
        }
        else
        {
            const vectorField& CellCenters = mesh.C();

            for (int i = 0; i < CellCenters.size(); i++)
            {
                vertices[verticesIndex++] = CellCenters[i].x();
                vertices[verticesIndex++] = CellCenters[i].y();
                if (dim_ == 3)
                {
                    vertices[verticesIndex++] = CellCenters[i].z();
                }
            }
        }

        // Get the locations of the mesh vertices (here: face centers)
        // for all the patches
        for (uint j = 0; j < patchIDs_.size(); j++)
        {
            // Get the face centers of the current patch
            const vectorField faceCenters =
                mesh.boundaryMesh()[patchIDs_.at(j)].faceCentres();

            // Assign the (x,y,z) locations to the vertices
            for (int i = 0; i < faceCenters.size(); i++)
            {
                vertices[verticesIndex++] = faceCenters[i].x();
                vertices[verticesIndex++] = faceCenters[i].y();
                if (dim_ == 3)
                {
                    vertices[verticesIndex++] = faceCenters[i].z();
                }
            }
        }

        // Pass the mesh vertices information to preCICE
        precice_.setMeshVertices(meshName_, vertices, vertexIDs_);
    }
}


void preciceAdapter::Interface::addCouplingDataWriter(
    std::string dataName,
    CouplingDataUser* couplingDataWriter)
{
    // Set the data name (from preCICE)
    couplingDataWriter->setDataName(dataName);

    // Set the patchIDs of the patches that form the interface
    couplingDataWriter->setPatchIDs(patchIDs_);

    // Set the names of the cell sets to be coupled (for volume coupling)
    couplingDataWriter->setCellSetNames(cellSetNames_);

    // Set the location type in the CouplingDataUser class
    couplingDataWriter->setLocationsType(locationType_);

    // Set the location type in the CouplingDataUser class
    couplingDataWriter->checkDataLocation(meshConnectivity_);

    // Initilaize class specific data
    couplingDataWriter->initialize();

    // Add the CouplingDataUser to the list of writers
    couplingDataWriters_.push_back(couplingDataWriter);
}


void preciceAdapter::Interface::addCouplingDataReader(
    std::string dataName,
    preciceAdapter::CouplingDataUser* couplingDataReader)
{
    // Set the patchIDs of the patches that form the interface
    couplingDataReader->setDataName(dataName);

    // Add the CouplingDataUser to the list of readers
    couplingDataReader->setPatchIDs(patchIDs_);

    // Set the location type in the CouplingDataUser class
    couplingDataReader->setLocationsType(locationType_);

    // Set the names of the cell sets to be coupled (for volume coupling)
    couplingDataReader->setCellSetNames(cellSetNames_);

    // Check, if the current location type is supported by the data type
    couplingDataReader->checkDataLocation(meshConnectivity_);

    // Initilaize class specific data
    couplingDataReader->initialize();

    // Add the CouplingDataUser to the list of readers
    couplingDataReaders_.push_back(couplingDataReader);
}

void preciceAdapter::Interface::createBuffer()
{
    // Will the interface buffer need to store 3D vector data?
    bool needsVectorData = false;
    int dataBufferSize = 0;

    // Check all the coupling data readers
    for (uint i = 0; i < couplingDataReaders_.size(); i++)
    {
        if (couplingDataReaders_.at(i)->hasVectorData())
        {
            needsVectorData = true;
        }
    }

    // Check all the coupling data writers
    for (uint i = 0; i < couplingDataWriters_.size(); i++)
    {
        if (couplingDataWriters_.at(i)->hasVectorData())
        {
            needsVectorData = true;
        }
    }

    // Set the appropriate buffer size
    if (needsVectorData)
    {
        dataBufferSize = dim_ * numDataLocations_;
    }
    else
    {
        dataBufferSize = numDataLocations_;
    }

    // Create the data buffer
    // An interface has only one data buffer, which is shared between several
    // CouplingDataUsers.
    // TODO: Check (write tests) if this works properly when we have multiple
    // scalar and vector coupling data users in an interface. With the current
    // preCICE implementation, it should work as, when writing scalars,
    // it should  only use the first 1/3 elements of the buffer.
    dataBuffer_.resize(dataBufferSize);
}

void preciceAdapter::Interface::readCouplingData(double relativeReadTime)
{
    // Make every coupling data reader read
    for (uint i = 0; i < couplingDataReaders_.size(); i++)
    {
        preciceAdapter::CouplingDataUser* couplingDataReader = couplingDataReaders_.at(i);

        int dataDim = precice_.getDataDimensions(meshName_, couplingDataReader->dataName());

        // === SCATTER-FROM-MASTER: Only master reads from preCICE, then scatter ===

        // Step 1: Master reads global data from preCICE
        if (Pstream::master())
        {
            std::size_t globalReadSize = globalNumDataLocations_ * dataDim;
            globalDataBuffer_.resize(globalReadSize);

            precice_.readData(
                meshName_,
                couplingDataReader->dataName(),
                globalVertexIDs_,
                relativeReadTime,
                {globalDataBuffer_.data(), globalReadSize});

            DEBUG(Pout << "Adapter [Master]: Read " << globalReadSize << " values from preCICE" << endl);
        }

        // Step 2: Scatter data from master to all ranks
        // Prepare receive buffer on each rank
        int localSize = numDataLocations_ * dataDim;
        dataBuffer_.resize(localSize);

        // Use Pstream gatherList in reverse (build lists on master, scatter)
        List<List<double>> allDataLists(Pstream::nProcs());

        if (Pstream::master())
        {
            // Master splits global buffer into per-rank lists
            int globalIdx = 0;
            for (int rank = 0; rank < Pstream::nProcs(); rank++)
            {
                int rankCount = gatherCounts_[rank] / dim_ * dataDim; // Adjust for dataDim
                allDataLists[rank].resize(rankCount);
                for (int j = 0; j < rankCount; j++)
                {
                    allDataLists[rank][j] = globalDataBuffer_[globalIdx++];
                }
            }
        }

        // Scatter the lists (broadcast from master)
        Pstream::broadcast(allDataLists);

        // Each rank copies its data from the list to local buffer
        const List<double>& myData = allDataLists[Pstream::myProcNo()];
        for (int j = 0; j < myData.size() && j < localSize; j++)
        {
            dataBuffer_[j] = myData[j];
        }

        DEBUG(Pout << "Adapter [Procid " << Pstream::myProcNo() << "]: Received "
                   << myData.size() << " values via scatter" << endl);

        // Step 3: Apply data to OpenFOAM fields
        couplingDataReader->read(dataBuffer_.data(), dim_);
    }
}

void preciceAdapter::Interface::writeCouplingData()
{
    // TODO: wrap around isWriteDataRequired
    // Does the participant need to write data or is it subcycling?
    // if (precice_.isWriteDataRequired(computedTimestepLength))
    // {

    // Make every coupling data writer write
    for (uint i = 0; i < couplingDataWriters_.size(); i++)
    {
        preciceAdapter::CouplingDataUser* couplingDataWriter = couplingDataWriters_.at(i);

        int dataDim = precice_.getDataDimensions(meshName_, couplingDataWriter->dataName());

        // === GATHER-TO-MASTER: Each rank fills local buffer, gather, master writes ===

        // Step 1: Each rank fills its local buffer
        auto nWrittenData = couplingDataWriter->write(dataBuffer_.data(), meshConnectivity_, dim_);

        DEBUG(Pout << "Adapter [Procid " << Pstream::myProcNo() << "]: Local buffer has "
                   << nWrittenData << " values to send" << endl);

        // Step 2: Gather all local buffers to master
        // Convert to List for Pstream compatibility
        List<double> localDataList(nWrittenData);
        for (std::size_t j = 0; j < nWrittenData; j++)
        {
            localDataList[j] = dataBuffer_[j];
        }

        List<List<double>> allDataLists(Pstream::nProcs());
        allDataLists[Pstream::myProcNo()] = localDataList;
        Pstream::gatherList(allDataLists);

        // Step 3: Master flattens gathered data and writes to preCICE
        if (Pstream::master())
        {
            // Calculate total size and flatten
            std::size_t globalSize = 0;
            for (int rank = 0; rank < Pstream::nProcs(); rank++)
            {
                globalSize += allDataLists[rank].size();
            }

            globalDataBuffer_.resize(globalSize);
            int globalIdx = 0;
            for (int rank = 0; rank < Pstream::nProcs(); rank++)
            {
                const List<double>& rankData = allDataLists[rank];
                forAll(rankData, j)
                {
                    globalDataBuffer_[globalIdx++] = rankData[j];
                }
            }

            DEBUG(Pout << "Adapter [Master]: Writing " << globalSize << " values to preCICE" << endl);

            // Master writes to preCICE
            precice_.writeData(
                meshName_,
                couplingDataWriter->dataName(),
                globalVertexIDs_,
                {globalDataBuffer_.data(), globalSize});
        }
    }
    // }
}

preciceAdapter::Interface::~Interface()
{
    // Delete all the coupling data readers
    for (uint i = 0; i < couplingDataReaders_.size(); i++)
    {
        delete couplingDataReaders_.at(i);
    }
    couplingDataReaders_.clear();

    // Delete all the coupling data writers
    for (uint i = 0; i < couplingDataWriters_.size(); i++)
    {
        delete couplingDataWriters_.at(i);
    }
    couplingDataWriters_.clear();
}
