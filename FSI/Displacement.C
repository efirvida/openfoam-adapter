#include "Displacement.H"

using namespace Foam;


preciceAdapter::FSI::Displacement::Displacement
(
    const Foam::fvMesh& mesh,
    const std::string namePointDisplacement
    /* TODO: We should add any required field names here.
    /  They would need to be vector fields.
    /  See CHT/Temperature.C for details.
    */
)
:
pointDisplacement_(
    const_cast<pointVectorField*>
    (
        &mesh.lookupObject<pointVectorField>(namePointDisplacement)
    )
)
{
    dataType_ = vector;
}



void preciceAdapter::FSI::Displacement::write(double * buffer)
{
    /* TODO: Implement
    * We need two nested for-loops for each patch,
    * the outer for the locations and the inner for the dimensions.
    * See the preCICE writeBlockVectorData() implementation.
    */
    FatalErrorInFunction
        << "Writing displacements is not supported."
        << exit(FatalError);
}

// return the displacement to use later in the velocity?
void preciceAdapter::FSI::Displacement::read(double * buffer)
{
    /* TODO: Implement
    * We need two nested for-loops for each patch,
    * the outer for the locations and the inner for the dimensions.
    * See the preCICE readBlockVectorData() implementation.
    */

    // For every element in the buffer
    int bufferIndex = 0;


    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);


        // get the displacement on the patch
        fixedValuePointPatchVectorField& pointDisplacementFluidPatch =
            refCast<fixedValuePointPatchVectorField>
            (
                pointDisplacement_->boundaryFieldRef()[patchID]
            );

        // Info << nl << "Displacement of  old point: " << pointDisplacementFluidPatch[160] << endl; //33 for flap_perp., 153
        // For every cell of the patch
        forAll(pointDisplacement_->boundaryFieldRef()[patchID], i)
        {
            // Set the buffer as the pointdisplacment value
            pointDisplacementFluidPatch[i][0] = buffer[bufferIndex++];
            pointDisplacementFluidPatch[i][1] = buffer[bufferIndex++];
            pointDisplacementFluidPatch[i][2] = buffer[bufferIndex++];
        }

        // Info << nl << "Displacement of watchpoint: " << pointDisplacementFluidPatch[160] << endl; //33 for flap_perp., 153
        // For every cell of the patch
        // Info << nl << "old value pointDisplacement: " << max(pointDisplacement_->oldTime().primitiveField()) << endl;
    }
}
