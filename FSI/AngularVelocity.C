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
    // For global data (single vertex), buffer contains exactly one scalar value
    // The dim parameter is the spatial dimension (2D/3D), not relevant for scalars
    // For scalar data, preCICE provides 1 value per vertex
    
    // Read the first (and only) value from the buffer
    omegaField_->value() = buffer[0];
    
    DEBUG(adapterInfo("Received angular velocity: " + std::to_string(buffer[0]) + " rad/s ("
                      + std::to_string(buffer[0] * 60.0 / (2.0 * pi)) + " RPM)"));
}

bool preciceAdapter::FSI::AngularVelocity::isLocationTypeSupported(const bool meshConnectivity) const
{
    // AngularVelocity is global data - it's a single scalar value not associated
    // with specific mesh locations. Support all location types since we only
    // use a single vertex at origin for global data exchange.
    return true;
}

std::string preciceAdapter::FSI::AngularVelocity::getDataName() const
{
    return "AngularVelocity";
}
