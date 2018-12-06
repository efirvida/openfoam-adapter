#include "Force.H"
#include "fvcGrad.H"
#include "porosityModel.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam;

preciceAdapter::FSI::Force::Force
(
    const Foam::fvMesh& mesh,
    const fileName& timeName,
    const std::string nameP,
    const std::string nameU,
    const std::string namefD,
    const bool directForceDensity,
    const bool porosity 
    /* TODO: We should add any required field names here.
    /  They would need to be vector fields.
    /  See CHT/Temperature.C for details.
    */
)
:
mesh_(mesh),
nameP_(nameP),
nameU_(nameU),
namefD_(namefD),
directForceDensity_(directForceDensity),
porosity_(porosity)
{
    dataType_ = vector;

    // TODO: Is this ok?
    Force_ = new volVectorField(
        IOobject(
            "Force",
            timeName,
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        mesh,
        dimensionedVector("fdim", dimForce, Zero));

}

Foam::tmp<Foam::volSymmTensorField> preciceAdapter::FSI::Force::devRhoReff() const
{
        typedef compressible::turbulenceModel cmpTurbModel;
        typedef incompressible::turbulenceModel icoTurbModel;

        if (mesh_.foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
        {
                const cmpTurbModel &turb =
                    mesh_.lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

                return turb.devRhoReff();
        }
        else if (mesh_.foundObject<icoTurbModel>(icoTurbModel::propertiesName))
        {
                const incompressible::turbulenceModel &turb =
                    mesh_.lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

                return rho() * turb.devReff();
        }
        else if (mesh_.foundObject<fluidThermo>(fluidThermo::dictName))
        {
                const fluidThermo &thermo =
                    mesh_.lookupObject<fluidThermo>(fluidThermo::dictName);

                const volVectorField &U = mesh_.lookupObject<volVectorField>(nameU_);

                return -thermo.mu() * dev(twoSymm(fvc::grad(U)));
        }
        else if (
            mesh_.foundObject<transportModel>(nameTransportProperties_))
        {
                const transportModel &laminarT =
                    mesh_.lookupObject<transportModel>(nameTransportProperties_);

                const volVectorField &U = mesh_.lookupObject<volVectorField>(nameU_);

                return -rho() * laminarT.nu() * dev(twoSymm(fvc::grad(U)));
        }
        else if (mesh_.foundObject<dictionary>(nameTransportProperties_))
        {
                const dictionary &transportProperties =
                    mesh_.lookupObject<dictionary>(nameTransportProperties_);

                dimensionedScalar nu(
                    nameNu_,
                    dimViscosity,
                    transportProperties.lookup(nameNu_));

                const volVectorField &U = mesh_.lookupObject<volVectorField>(nameU_);

                return -rho() * nu * dev(twoSymm(fvc::grad(U)));
        }
        else
        {
                FatalErrorInFunction
                    << "No valid model for viscous stress calculation"
                    << exit(FatalError);

                return volSymmTensorField::null();
        }
}

Foam::tmp<Foam::volScalarField> preciceAdapter::FSI::Force::mu() const
{
    if (mesh_.foundObject<fluidThermo>(basicThermo::dictName))
    {
        const fluidThermo &thermo = mesh_.lookupObject<fluidThermo>(basicThermo::dictName);

        return thermo.mu();
    }
    else if (mesh_.foundObject<transportModel>(nameTransportProperties_))
    {
        const transportModel &laminarT = mesh_.lookupObject<transportModel>(nameTransportProperties_);

        return rho() * laminarT.nu();
    }
    else if (mesh_.foundObject<dictionary>(nameTransportProperties_))
    {
        const dictionary &transportProperties = mesh_.lookupObject<dictionary>(nameTransportProperties_);

        dimensionedScalar nu(nameNu_, dimViscosity, transportProperties.lookup(nameNu_));

        return rho() * nu;
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for dynamic viscosity calculation"
            << exit(FatalError);

        return volScalarField::null();
    }
}

Foam::tmp<Foam::volScalarField> preciceAdapter::FSI::Force::rho() const
{
    if (nameRho_ == nameRhoInf_)
    {
        return tmp<volScalarField>(
            new volScalarField(
                IOobject(
                    nameRho_,
                    mesh_.time().timeName(),
                    mesh_),
                mesh_,
                dimensionedScalar(nameRho_, dimDensity, rhoRef_)));
    }
    else
    {
        return (mesh_.lookupObject<volScalarField>(nameRho_));
    }
}

Foam::scalar preciceAdapter::FSI::Force::rho(const volScalarField &p) const
{
        if (p.dimensions() == dimPressure)
        {
            return 1.0;
        }
        else
        {
            if (nameRho_ != nameRhoInf_)
            {
                FatalErrorInFunction
                    << "Dynamic pressure is expected but kinematic is provided."
                    << exit(FatalError);
            }

            return rhoRef_;
        }
}

void preciceAdapter::FSI::Force::addToFields(
    const label patchID,
    const vectorField &fN,
    const vectorField &fT,
    const vectorField &fP)
{
        volVectorField &force = const_cast<volVectorField &>(mesh_.lookupObject<volVectorField>("Force"));

        vectorField &pf = force.boundaryFieldRef()[patchID];
        pf += fN + fT + fP;

}

void preciceAdapter::FSI::Force::addToFields(
    const labelList &cellIDs,
    const vectorField &fN,
    const vectorField &fT,
    const vectorField &fP)
{
    volVectorField &force =
        const_cast<volVectorField &>(
            mesh_.lookupObject<volVectorField>("Force"));


    forAll(cellIDs, i)
    {
        label celli = cellIDs[i];
        force[celli] += fN[i] + fT[i] + fP[i];
    }
}

void preciceAdapter::FSI::Force::calcForcesMoment()
{
    Force_[0] = Foam::vector::zero;
    Force_[1] = Foam::vector::zero;
    Force_[2] = Foam::vector::zero;

    if (mesh_.foundObject<transportModel>(nameTransportProperties_))
    {
        nameRho_ = nameRhoInf_;
        const dictionary &transportProperties = mesh_.lookupObject<IOdictionary>(nameTransportProperties_);

        dimensionedScalar rhoRef_(
            nameRho_,
            dimDensity,
            transportProperties.lookup(nameRho_));
    }
    
    if (directForceDensity_)
    {
        const volVectorField &fD = mesh_.lookupObject<volVectorField>(namefD_);

        const surfaceVectorField::Boundary &Sfb = mesh_.Sf().boundaryField();

        for (uint j = 0; j < patchIDs_.size(); j++)
        {
            int patchID = patchIDs_.at(j);

            scalarField sA(mag(Sfb[patchID]));

            // Normal force = surfaceUnitNormal*(surfaceNormal & forceDensity)
            vectorField fN(Sfb[patchID] / sA * (Sfb[patchID] & fD.boundaryField()[patchID]));

            // Tangential force (total force minus normal fN)
            vectorField fT(sA * fD.boundaryField()[patchID] - fN);

            //- Porous force
            vectorField fP(fT.size(), Zero);

            addToFields(patchID, fN, fT, fP);
        }
    }
    else
    {
        const volScalarField &p = mesh_.lookupObject<volScalarField>(nameP_);

        const surfaceVectorField::Boundary &Sfb = mesh_.Sf().boundaryField();

        tmp<volSymmTensorField> tdevRhoReff = devRhoReff();
        const volSymmTensorField::Boundary &devRhoReffb = tdevRhoReff().boundaryField();

        for (uint j = 0; j < patchIDs_.size(); j++)
        {
            int patchID = patchIDs_.at(j);

            vectorField fN(rho(p) * Sfb[patchID] * (p.boundaryField()[patchID]));

            vectorField fT(Sfb[patchID] & devRhoReffb[patchID]);

            vectorField fP(fT.size(), Zero);

            addToFields(patchID, fN, fT, fP);
        }
    }

    if (porosity_)
    {
        const volVectorField &U = mesh_.lookupObject<volVectorField>(nameU_);
        const volScalarField rho(this->rho());
        const volScalarField mu(this->mu());

        const HashTable<const porosityModel *> models = mesh_.lookupClass<porosityModel>();

        if (models.empty())
        {
            WarningInFunction
                << "Porosity effects requested, but no porosity models found "
                << "in the database"
                << endl;
        }

        forAllConstIter(HashTable<const porosityModel *>, models, iter)
        {
            // non-const access required if mesh is changing
            porosityModel &pm = const_cast<porosityModel &>(*iter());

            vectorField fPTot(pm.force(U, rho, mu));

            const labelList &cellZoneIDs = pm.cellZoneIDs();

            forAll(cellZoneIDs, i)
            {
                label zoneI = cellZoneIDs[i];
                const cellZone &cZone = mesh_.cellZones()[zoneI];

                const vectorField d(mesh_.C(), cZone);
                const vectorField fP(fPTot, cZone);

                const vectorField fDummy(fP.size(), Zero);

                addToFields(cZone, fDummy, fDummy, fP);
            }
        }
    }
}

void preciceAdapter::FSI::Force::write(double *buffer)
{
    calcForcesMoment();

    mesh_.lookupObject<volVectorField>("Force").write();
}

void preciceAdapter::FSI::Force::read(double *buffer)
{
       /* TODO: Implement
    * We need two nested for-loops for each patch,
    * the outer for the locations and the inner for the dimensions.
    * See the preCICE readBlockVectorData() implementation.
    */
    FatalErrorInFunction
        << "Reading forces is not supported."
        << exit(FatalError);
}

preciceAdapter::FSI::Force::~Force()
{
    // TODO: Is this enough?
    delete Force_;
}
