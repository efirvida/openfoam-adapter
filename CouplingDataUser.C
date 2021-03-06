#include "CouplingDataUser.H"

preciceAdapter::CouplingDataUser::CouplingDataUser()
{
}

bool preciceAdapter::CouplingDataUser::hasScalarData()
{
    return dataType_ == scalar;
}

bool preciceAdapter::CouplingDataUser::hasVectorData()
{
    return dataType_ == vector;
}

void preciceAdapter::CouplingDataUser::setDataID(int dataID)
{
    dataID_ = dataID;

    return;
}

int preciceAdapter::CouplingDataUser::dataID()
{
    return dataID_;
}

void preciceAdapter::CouplingDataUser::setPatchIDs(std::vector<int> patchIDs)
{
    patchIDs_ = patchIDs;
}
