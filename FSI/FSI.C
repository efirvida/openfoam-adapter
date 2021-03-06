#include "FSI.H"

#include "Utilities.H"

using namespace Foam;

preciceAdapter::FSI::FluidStructureInteraction::FluidStructureInteraction
(
    const Foam::fvMesh& mesh,
    const Foam::Time& runTime
)
:
mesh_(mesh),
runTime_(runTime)
{}

bool preciceAdapter::FSI::FluidStructureInteraction::configure(const YAML::Node adapterConfig)
{
    DEBUG(adapterInfo("Configuring the FSI module..."));

    // Read the FSI-specific options from the adapter's configuration file
    if (!readConfig(adapterConfig)) return false;

    /* TODO: If we need different solver types,
    /  here is the place to determine it.
    */

    return true;
}

bool preciceAdapter::FSI::FluidStructureInteraction::readConfig(const YAML::Node adapterConfig)
{
    /* TODO: Read the solver type, if needed.
    /  If you want to determine it automatically, implement a method
    /  as in CHT/CHT.C
    */

    /* TODO: Read the names of any needed fields and parameters.
    * Include the force here?
    */

    // Read the name of the pointDisplacement field (if different)
    if (adapterConfig["namePointDisplacement"])
    {
        namePointDisplacement_ = adapterConfig["namePointDisplacement"].as<std::string>();
    }
    DEBUG(adapterInfo("    pointDisplacement field name : " + namePointDisplacement_));

    if (adapterConfig["p"])
    {
        nameP_ = adapterConfig["p"].as<std::string>();
    }
    DEBUG(adapterInfo("    Pressure field name : " + nameP_));

    if (adapterConfig["U"])
    {
        nameU_ = adapterConfig["U"].as<std::string>();
    }
    DEBUG(adapterInfo("    pointDisplacement field name : " + nameU_));

    if (adapterConfig["porosity"])
    {
        porosity_ = adapterConfig["porosity"].as<bool>();
    }
    DEBUG(adapterInfo("    porosity enabled : " + porosity_));

    if (adapterConfig["directForceDensity"] && adapterConfig["fD"])
    {
        directForceDensity_ = adapterConfig["directForceDensity"].as<bool>();
        namefD_ = adapterConfig["fD"].as<std::string>();
        DEBUG(adapterInfo("    force density field being supplied: " + directForceDensity_));
        DEBUG(adapterInfo("    force density (fD) field name : " + namefD_));
    }
    else if (adapterConfig["directForceDensity"] && !adapterConfig["fD"])
    {
        FatalErrorInFunction << "    force density field name has to be provided" << exit(FatalError);
    } 

    return true;
}

void preciceAdapter::FSI::FluidStructureInteraction::addWriters(std::string dataName, Interface * interface)
{
    /* TODO: Add writers. See CHT/CHT.C for reference.
    /  We probably need to do this for displacements and forces.
    /  If different coupling data users per solver type are defined,
    /  we need to check for that here.
    */
    if (dataName.find("Force") == 0)
    {
        interface->addCouplingDataWriter
        (
            dataName, new Force(
                mesh_, 
                runTime_.timeName(),
                nameP_,
                nameU_,
                namefD_,
                directForceDensity_,
                porosity_
            ) 
            /* TODO: Add any other arguments here */
        );
        DEBUG(adapterInfo("Added writer: Force."));
    }
    // TODO MOVE THIS UP
    else if (dataName.find("DisplacementDelta") == 0)
    {
        interface->addCouplingDataWriter
        (
            dataName,
            new DisplacementDelta(mesh_, namePointDisplacement_) /* TODO: Add any other arguments here */
        );
        DEBUG(adapterInfo("Added writer: DisplacementDelta."));
    }
    // TODO Do we need to include the displacement and velocity? They will never be written...
    else if (dataName.find("Displacement") == 0)
    {
        interface->addCouplingDataWriter
        (
            dataName,
            new Displacement(mesh_, namePointDisplacement_) /* TODO: Add any other arguments here */
        );
        DEBUG(adapterInfo("Added writer: Displacement."));
    }

    // NOTE: If you want to couple another variable, you need
    // to add your new coupling data user as a coupling data
    // writer here (and as a reader below).
    // The argument of the dataName.compare() needs to match
    // the one provided in the adapter's configuration file.
}

void preciceAdapter::FSI::FluidStructureInteraction::addReaders(std::string dataName, Interface * interface)
{
    /* TODO: Add readers. See CHT/CHT.C for reference.
    /  We probably need to do this for displacements and forces.
    /  If different coupling data users per solver type are defined,
    /  we need to check for that here.
    */

    // TODO do we need to include the force here, since it will not be read by openFOAM?
    if (dataName.find("Force") == 0)
    {
        interface->addCouplingDataReader
        (
            dataName, new Force(
                mesh_, 
                runTime_.timeName(),
                nameP_,
                nameU_,
                namefD_,
                directForceDensity_,
                porosity_
            ) /* TODO: Add any other arguments here */
        );
        DEBUG(adapterInfo("Added reader: Force."));
    }
    else if (dataName.find("DisplacementDelta") == 0)
    {
        interface->addCouplingDataReader
        (
            dataName,
            new DisplacementDelta(mesh_, namePointDisplacement_) /* TODO: Add any other arguments here */
        );
        DEBUG(adapterInfo("Added reader: DisplacementDelta."));
    }
    else if (dataName.find("Displacement") == 0)
    {
        interface->addCouplingDataReader
        (
            dataName,
            new Displacement(mesh_, namePointDisplacement_) /* TODO: Add any other arguments here */
        );
        DEBUG(adapterInfo("Added reader: Displacement."));
    }
    // NOTE: If you want to couple another variable, you need
    // to add your new coupling data user as a coupling data
    // writer here (and as a writer above).
    // The argument of the dataName.compare() needs to match
    // the one provided in the adapter's configuration file.
}
