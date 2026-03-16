#include "AngularVelocity.H"
#include "Pstream.H"

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
    
    // Register the omega field in Time, not in mesh.
    // This ensures that preciceOmega (which reads from Time) can find it.
    const Time& runTime = mesh_.time();
    
    // Check if an omega field with the requested name exists in Time registry.
    // If yes, bind omegaField_ to that field.
    // If not, create it.
    if (runTime.foundObject<uniformDimensionedScalarField>(fieldName))
    {
        omegaField_ = 
            &const_cast<uniformDimensionedScalarField&>(
                runTime.lookupObject<uniformDimensionedScalarField>(fieldName));
        
        DEBUG(adapterInfo("Found existing omega field '" + omegaFieldName_ + "' in Time registry"));
    }
    else
    {
        omegaFieldOwning_.reset(new uniformDimensionedScalarField(
            IOobject(
                fieldName,
                runTime.constant(),
                runTime,  // Register in Time, not mesh
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
        
        DEBUG(adapterInfo("Created omega field '" + omegaFieldName_ + "' in Time registry"));
    }
}

std::size_t preciceAdapter::FSI::AngularVelocity::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    adapterInfo("Writing angular velocity is not supported.", "error");
    return 0;
}

void preciceAdapter::FSI::AngularVelocity::read(double* buffer, const unsigned int dim)
{
    // Only the master rank owns the single global vertex and has a valid buffer.
    if (!Pstream::parRun() || Pstream::master())
    {
        omegaField_->value() = buffer[0];
    }
    // Broadcast the value to all secondary ranks.
    if (Pstream::parRun())
    {
        Pstream::broadcast(omegaField_->value());
    }
    
    DEBUG(adapterInfo("Received angular velocity: " + std::to_string(omegaField_->value()) + " rad/s ("
                      + std::to_string(omegaField_->value() * 60.0 / (2.0 * pi)) + " RPM)"));
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
