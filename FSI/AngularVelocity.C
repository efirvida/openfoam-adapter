#include "AngularVelocity.H"
#include "mathematicalConstants.H"

using namespace Foam;
using namespace Foam::constant::mathematical;

preciceAdapter::FSI::AngularVelocity::AngularVelocity(
    const Foam::fvMesh& mesh,
    const std::string omegaFieldName,
    const Foam::vector& rotationAxis)
: mesh_(mesh),
  omegaFieldName_(omegaFieldName),
  rotationAxis_(rotationAxis / mag(rotationAxis))
{
    dataType_ = scalar;
    
    const word fieldName(omegaFieldName_);
    
    // Check if an omega field with the requested name exists.
    // If yes, bind omegaField_ to that field.
    // If not, create it.
    if (mesh_.foundObject<uniformDimensionedScalarField>(fieldName))
    {
        omegaField_ = 
            &const_cast<uniformDimensionedScalarField&>(
                mesh_.lookupObject<uniformDimensionedScalarField>(fieldName));
    }
    else
    {
        omegaFieldOwning_.reset(new uniformDimensionedScalarField(
            IOobject(
                fieldName,
                mesh_.time().constant(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            dimensionedScalar(
                fieldName,
                dimensionSet(0, 0, -1, 0, 0, 0, 0),
                0.0
            )
        ));

        omegaField_ = omegaFieldOwning_.get();
    }
}

std::size_t preciceAdapter::FSI::AngularVelocity::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    adapterInfo("Writing angular velocity is not supported.", "error");
    return 0;
}

void preciceAdapter::FSI::AngularVelocity::read(double* buffer, const unsigned int dim)
{
    if (dim != 1)
    {
        adapterInfo("AngularVelocity expects a scalar value (omega in rad/s), but dim=" + std::to_string(dim), "error");
        return;
    }
    
    // Read the first value from the buffer
    omegaField_->value() = buffer[0];
    
    DEBUG(adapterInfo("Received angular velocity: " + std::to_string(buffer[0]) + " rad/s ("
                      + std::to_string(buffer[0] * 60.0 / (2.0 * pi)) + " RPM)"));
}

bool preciceAdapter::FSI::AngularVelocity::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::faceCenters || 
            this->locationType_ == LocationType::faceNodes);
}

std::string preciceAdapter::FSI::AngularVelocity::getDataName() const
{
    return "AngularVelocity";
}
