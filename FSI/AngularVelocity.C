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
  currentOmega_(0.0),
  rotationAxis_(rotationAxis / mag(rotationAxis)) // Normalizar
{
    dataType_ = scalar;
    
    // Buscar o crear el campo uniformDimensionedScalarField
    const word fieldName(omegaFieldName_);
    
    if (mesh_.foundObject<uniformDimensionedScalarField>(fieldName))
    {
        // El campo ya existe, usarlo
        omegaField_ = 
            &const_cast<uniformDimensionedScalarField&>(
                mesh_.lookupObject<uniformDimensionedScalarField>(fieldName));
        
        DEBUG(adapterInfo("Found existing angular velocity field: " + omegaFieldName_));
    }
    else
    {
        // Crear nuevo campo uniformDimensionedScalarField
        omegaFieldOwning_.reset(new uniformDimensionedScalarField(
            IOobject(
                fieldName,
                mesh_.time().constant(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            dimensionedScalar(
                fieldName,
                dimensionSet(0, 0, -1, 0, 0, 0, 0), // rad/s [1/T]
                0.0
            )
        ));
        
        omegaField_ = omegaFieldOwning_.get();
        
        // Registrar en el object registry para que otros modelos puedan acceder
        omegaField_->store();
        
        DEBUG(adapterInfo("Created angular velocity field: " + omegaFieldName_));
    }
    
    DEBUG(adapterInfo("AngularVelocity coupling initialized:"));
    DEBUG(adapterInfo("  Field name: " + omegaFieldName_));
    DEBUG(adapterInfo("  Rotation axis: (" + std::to_string(rotationAxis_.x()) + ", "
                      + std::to_string(rotationAxis_.y()) + ", "
                      + std::to_string(rotationAxis_.z()) + ")"));
}

preciceAdapter::FSI::AngularVelocity::~AngularVelocity()
{
    // Si creamos el campo, se limpiará automáticamente por unique_ptr
    // pero si está registrado, debemos tener cuidado
}

void preciceAdapter::FSI::AngularVelocity::initialize()
{
    // No se necesita inicialización especial
}

std::size_t preciceAdapter::FSI::AngularVelocity::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    adapterInfo("Writing angular velocity is not supported.", "error");
    return 0;
}

void preciceAdapter::FSI::AngularVelocity::read(double* buffer, const unsigned int dim)
{
    // Leer la velocidad angular desde el buffer
    // preCICE envía omega en rad/s como un valor escalar
    
    if (dim != 1)
    {
        adapterInfo("AngularVelocity expects a scalar value (omega in rad/s), but dim=" + std::to_string(dim), "error");
        return;
    }
    
    // Leer el primer valor del buffer (asumimos que es el mismo para toda la interfaz)
    currentOmega_ = buffer[0];
    
    DEBUG(adapterInfo("Received angular velocity: " + std::to_string(currentOmega_) + " rad/s ("
                      + std::to_string(currentOmega_ * 60.0 / (2.0 * pi)) + " RPM)"));
    
    // Actualizar el campo uniformDimensionedScalarField
    omegaField_->value() = currentOmega_;
    
    DEBUG(adapterInfo("Updated field '" + omegaFieldName_ + "' = " + std::to_string(currentOmega_) + " rad/s"));
}

bool preciceAdapter::FSI::AngularVelocity::isLocationTypeSupported(const bool meshConnectivity) const
{
    // Para velocidad angular, soportamos cualquier tipo de localización
    // ya que es un valor escalar uniforme
    return (this->locationType_ == LocationType::faceCenters || 
            this->locationType_ == LocationType::faceNodes);
}

std::string preciceAdapter::FSI::AngularVelocity::getDataName() const
{
    return "AngularVelocity";
}
