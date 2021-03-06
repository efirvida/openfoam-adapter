#ifndef FSI_DISPLACEMENTDELTA_H
#define FSI_DISPLACEMENTDELTA_H

#include "CouplingDataUser.H"

#include "fvCFD.H"
#include "fixedValuePointPatchFields.H"

#include "pointFields.H"
#include "vectorField.H"

namespace preciceAdapter
{
namespace FSI
{

//- Class that writes and reads displacement
class DisplacementDelta : public CouplingDataUser
{

private:

    /* TODO: We probably need to declare some fields here,
    /  see CHT/Temperature.H.
    */

    // Displacement pointVectorField
    Foam::pointVectorField * pointDisplacement_;

public:

    //- Constructor
    DisplacementDelta
    (
        const Foam::fvMesh& mesh,
        const std::string namePointDisplacement
        /* TODO: We should add any required field names here.
        /  They would need to be vector fields.
        /  See CHT/Temperature.H for details.
        /  Apply any changes also to DisplacementDelta.C.
        */
    );

    //- Write the displacementDelta values into the buffer
    void write(double * buffer);

    //- Read the displacementDelta values from the buffer
    void read(double * buffer);

};

}
}

#endif
