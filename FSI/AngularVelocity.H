#ifndef FSI_ANGULARVELOCITY_H
#define FSI_ANGULARVELOCITY_H

#include "CouplingDataUser.H"
#include "fvCFD.H"

namespace preciceAdapter
{
namespace FSI
{

//- Class that reads angular velocity (Omega) from preCICE and stores it in a scalar field
class AngularVelocity : public CouplingDataUser
{
private:
    //- Reference to mesh
    const Foam::fvMesh& mesh_;

    //- Angular velocity field (UniformDimensionedScalarField)
    //  Stored in the object registry for use by other models
    Foam::uniformDimensionedScalarField* omegaField_;

    //- Owning pointer if we create the field
    std::unique_ptr<Foam::uniformDimensionedScalarField> omegaFieldOwning_;

    //- Name of the angular velocity field
    std::string omegaFieldName_;

    //- Current angular velocity value (rad/s)
    Foam::scalar currentOmega_;

    //- Axis of rotation (for vector omega calculation if needed)
    Foam::vector rotationAxis_;

public:
    //- Constructor
    AngularVelocity(
        const Foam::fvMesh& mesh,
        const std::string omegaFieldName,
        const Foam::vector& rotationAxis);

    //- Destructor
    ~AngularVelocity();

    //- Write angular velocity (not supported)
    std::size_t write(double* buffer, bool meshConnectivity, const unsigned int dim) final;

    //- Read angular velocity from preCICE and update the field
    void read(double* buffer, const unsigned int dim) final;

    bool isLocationTypeSupported(const bool meshConnectivity) const final;

    //- Get the name of the current data field
    std::string getDataName() const final;

    //- Initialize
    void initialize() final;

    //- Get current omega value
    Foam::scalar getCurrentOmega() const { return currentOmega_; }

    //- Get omega field name
    const std::string& getOmegaFieldName() const { return omegaFieldName_; }
};

}
}

#endif
