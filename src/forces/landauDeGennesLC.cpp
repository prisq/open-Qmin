#include "landauDeGennesLC.h"
#include "landauDeGennesLC.cuh"
#include "qTensorFunctions.h"
#include "utilities.cuh"
/*! \file landauDeGennesLC.cpp */

landauDeGennesLC::landauDeGennesLC(bool _neverGPU)
    {
    neverGPU = _neverGPU;
    if(neverGPU)
        {
        energyDensity.noGPU =true;
        objectForceArray.noGPU = true;
        forceCalculationAssist.noGPU=true;
        energyPerParticle.noGPU = true;
        }
    baseInitialization();
    }

landauDeGennesLC::landauDeGennesLC(scalar _A, scalar _B, scalar _C, scalar _L1) :
    A(_A), B(_B), C(_C), L1(_L1)
    {
    baseInitialization();
    setNumberOfConstants(distortionEnergyType::oneConstant);
    };

landauDeGennesLC::landauDeGennesLC(scalar _A, scalar _B, scalar _C, scalar _L1, scalar _L2,scalar _L3orWavenumber,
                                   distortionEnergyType _type) :
                                   A(_A), B(_B), C(_C), L1(_L1), L2(_L2), L3(_L3orWavenumber), q0(_L3orWavenumber)
    {
    baseInitialization();
    setNumberOfConstants(_type);
    };

landauDeGennesLC::landauDeGennesLC(scalar _A, scalar _B, scalar _C, scalar _L1, scalar _L2,scalar _L3, scalar _L4, scalar _L6) :
                                   A(_A), B(_B), C(_C), L1(_L1), L2(_L2), L3(_L3), L4(_L4), L6(_L6)
    {
    baseInitialization();
    setNumberOfConstants(distortionEnergyType::multiConstant);
    };

void landauDeGennesLC::baseInitialization()
    {
    useNeighborList=false;
    computeEfieldContribution=false;
    computeHfieldContribution=false;
    forceTuner = make_shared<kernelTuner>(128,256,32,10,200000);
    boundaryForceTuner = make_shared<kernelTuner>(128,256,32,10,200000);
    l24ForceTuner = make_shared<kernelTuner>(128,256,32,10,200000);
    fieldForceTuner = make_shared<kernelTuner>(128,256,32,10,200000);
    forceAssistTuner = make_shared<kernelTuner>(128,256,32,10,200000);
    energyComponents.resize(5);
    }

void landauDeGennesLC::setModel(shared_ptr<cubicLattice> _model)
    {
    lattice=_model;
    model = _model;
    if(numberOfConstants == distortionEnergyType::multiConstant)
        {
        lattice->fillNeighborLists(1);//fill neighbor lists to allow computing mixed partials
        }
    else // one constant approx
        {
        lattice->fillNeighborLists(0);
        }
    int N = lattice->getNumberOfParticles();
    energyDensity.resize(N);
    /*
    forceCalculationAssist.resize(N);
    if(useGPU)
        {
        ArrayHandle<cubicLatticeDerivativeVector> fca(forceCalculationAssist,access_location::device,access_mode::overwrite);
        cubicLatticeDerivativeVector zero(0.0);
        gpu_set_array(fca.data,zero,N,512);
        }
    else
        {
        ArrayHandle<cubicLatticeDerivativeVector> fca(forceCalculationAssist);
        cubicLatticeDerivativeVector zero(0.0);
        for(int ii = 0; ii < N; ++ii)
            fca.data[ii] = zero;
        };
    */
    };

void landauDeGennesLC::computeForces(GPUArray<dVec> &forces,bool zeroOutForce, int type)
    {
    if(useGPU)
        computeForceGPU(forces,zeroOutForce);
    else
        computeForceCPU(forces,zeroOutForce,type);
    }

void landauDeGennesLC::computeForceGPU(GPUArray<dVec> &forces,bool zeroOutForce)
    {
    int N = lattice->getNumberOfParticles();
    ArrayHandle<dVec> d_force(forces,access_location::device,access_mode::readwrite);
    ArrayHandle<dVec> d_spins(lattice->returnPositions(),access_location::device,access_mode::read);
    ArrayHandle<int>  d_latticeTypes(lattice->returnTypes(),access_location::device,access_mode::read);
    ArrayHandle<int>  d_latticeNeighbors(lattice->neighboringSites,access_location::device,access_mode::read);
    switch (numberOfConstants)
        {
        case distortionEnergyType::oneConstant :
            {
            forceTuner->begin();
            gpu_qTensor_oneConstantForce(d_force.data, d_spins.data, d_latticeTypes.data, d_latticeNeighbors.data,
                                         lattice->neighborIndex,
                                         A,B,C,L1,N,
                                         zeroOutForce,forceTuner->getParameter());
            break;
            };
        case distortionEnergyType::multiConstant:
            {
            bool zeroForce = zeroOutForce;
            computeFirstDerivatives();
            break;
            };
        };
    forceTuner->end();
    if(lattice->boundaries.getNumElements() >0)
        {
        computeBoundaryForcesGPU(forces,false);
        };
    if(computeEfieldContribution)
        computeEorHFieldForcesGPU(forces,false,Efield,deltaEpsilon,epsilon0);
    if(computeHfieldContribution)
        computeEorHFieldForcesGPU(forces,false,Hfield,deltaChi,mu0);
    };

void landauDeGennesLC::computeForceCPU(GPUArray<dVec> &forces,bool zeroOutForce, int type)
    {
    switch (numberOfConstants)
        {
        case distortionEnergyType::oneConstant :
            {
            if(type ==0)
                computeL1BulkCPU(forces,zeroOutForce);
            if(type ==1)
                computeL1BoundaryCPU(forces,zeroOutForce);
            break;
            };
        case distortionEnergyType::multiConstant :
            {
            bool zeroForce = zeroOutForce;
            computeFirstDerivatives();
            if(L1 != 0)
                {
                if(type ==0)
                    computeL1BulkCPU(forces,zeroForce);
                if(type ==1)
                    computeL1BoundaryCPU(forces,zeroForce);
                zeroForce = false;
                }
            if(L2 != 0)
                {
                if(type ==0)
                    computeL2BulkCPU(forces,zeroOutForce);
                if(type ==1)
                    computeL2BoundaryCPU(forces,zeroOutForce);
                zeroForce = false;
                }
            if(L3 != 0)
                {
                if(type ==0)
                    computeL3BulkCPU(forces,zeroOutForce);
                if(type ==1)
                    computeL3BoundaryCPU(forces,zeroOutForce);
                zeroForce = false;
                }
            if(L4 != 0)
                {
                if(type ==0)
                    computeL4BulkCPU(forces,zeroOutForce);
                if(type ==1)
                    computeL4BoundaryCPU(forces,zeroOutForce);
                zeroForce = false;
                }
            if(L6 != 0)
                {
                if(type ==0)
                    computeL6BulkCPU(forces,zeroOutForce);
                if(type ==1)
                    computeL6BoundaryCPU(forces,zeroOutForce);
                zeroForce = false;
                }
            break;
            };
        };
    if(type != 0 )
        {
        if(lattice->boundaries.getNumElements() >0)
            {
            computeBoundaryForcesCPU(forces,false);
            };
        if(computeEfieldContribution)
            {
            computeEorHFieldForcesCPU(forces,false, Efield,deltaEpsilon,epsilon0);
            };
        if(computeHfieldContribution)
            {
            computeEorHFieldForcesCPU(forces,false,Hfield,deltaChi,mu0);
            };
        };
    };

void landauDeGennesLC::computeEnergyGPU(bool verbose)
    {
    computeEnergyCPU(verbose);
    /*
    int N = lattice->getNumberOfParticles();
    ArrayHandle<dVec> d_force(forces,access_location::device,access_mode::readwrite);
    ArrayHandle<dVec> d_spins(lattice->returnPositions(),access_location::device,access_mode::read);
    ArrayHandle<int>  d_latticeTypes(lattice->returnTypes(),access_location::device,access_mode::read);
    ArrayHandle<int>  d_latticeNeighbors(lattice->neighboringSites,access_location::device,access_mode::read);
            forceTuner->begin();
            gpu_qTensor_oneConstantForce(d_force.data, d_spins.data, d_latticeTypes.data, d_latticeNeighbors.data,
                                         lattice->neighborIndex,
                                         A,B,C,L1,N,
                                         zeroOutForce,forceTuner->getParameter());
    */
    }


void landauDeGennesLC::computeEnergyCPU(bool verbose)
    {
    scalar phaseEnergy = 0.0;
    scalar distortionEnergy = 0.0;
    scalar anchoringEnergy = 0.0;
    scalar eFieldEnergy = 0.0;
    scalar hFieldEnergy = 0.0;
    energy=0.0;
    ArrayHandle<dVec> Qtensors(lattice->returnPositions());
    ArrayHandle<int> latticeTypes(lattice->returnTypes());
    ArrayHandle<boundaryObject> bounds(lattice->boundaries);
    ArrayHandle<scalar> energyPerSite(energyDensity);
    scalar a = 0.5*A;
    scalar b = B/3.0;
    scalar c = 0.25*C;
    int LCSites = 0;
    for (int i = 0; i < lattice->getNumberOfParticles(); ++i)
        {
        energyPerSite.data[i] = 0.0;
        //the current scheme for getting the six nearest neighbors
        int neighNum;
        vector<int> neighbors(6);
        int currentIndex;
        dVec qCurrent, xDown, xUp, yDown,yUp,zDown,zUp;
        currentIndex = lattice->getNeighbors(i,neighbors,neighNum);
        qCurrent = Qtensors.data[currentIndex];
        if(latticeTypes.data[currentIndex] <=0)
            {
            LCSites +=1;
            scalar phaseAtSite = a*TrQ2(qCurrent) + b*TrQ3(qCurrent) + c* TrQ2Squared(qCurrent);
            energyPerSite.data[i] += phaseAtSite;
            phaseEnergy += phaseAtSite;

            if(computeEfieldContribution)
                {
                    scalar eFieldAtSite = epsilon0*(-0.5*Efield.x*Efield.x*(epsilon + deltaEpsilon*qCurrent[0]) -
                              deltaEpsilon*Efield.x*Efield.y*qCurrent[1] - deltaEpsilon*Efield.x*Efield.z*qCurrent[2] -
                              0.5*Efield.z*Efield.z*(epsilon - deltaEpsilon*qCurrent[0] - deltaEpsilon*qCurrent[3]) -
                              0.5*Efield.y*Efield.y*(epsilon + deltaEpsilon*qCurrent[3]) - deltaEpsilon*Efield.y*Efield.z*qCurrent[4]);
                    eFieldEnergy+=eFieldAtSite;
                    energyPerSite.data[i] +=eFieldAtSite;
                }
            if(computeHfieldContribution)
                {
                    scalar hFieldAtSite=mu0*(-0.5*Hfield.x*Hfield.x*(Chi + deltaChi*qCurrent[0]) -
                              deltaChi*Hfield.x*Hfield.y*qCurrent[1] - deltaChi*Hfield.x*Hfield.z*qCurrent[2] -
                              0.5*Hfield.z*Hfield.z*(Chi - deltaChi*qCurrent[0] - deltaChi*qCurrent[3]) -
                              0.5*Hfield.y*Hfield.y*(Chi + deltaChi*qCurrent[3]) - deltaChi*Hfield.y*Hfield.z*qCurrent[4]);
                    hFieldEnergy+=hFieldAtSite;
                    energyPerSite.data[i] +=hFieldAtSite;
                }

            xDown = Qtensors.data[neighbors[0]];
            xUp = Qtensors.data[neighbors[1]];
            yDown = Qtensors.data[neighbors[2]];
            yUp = Qtensors.data[neighbors[3]];
            zDown = Qtensors.data[neighbors[4]];
            zUp = Qtensors.data[neighbors[5]];

            dVec firstDerivativeX = 0.5*(xUp - xDown);
            dVec firstDerivativeY = 0.5*(yUp - yDown);
            dVec firstDerivativeZ = 0.5*(zUp - zDown);
            scalar anchoringEnergyAtSite = 0.0;
            if(latticeTypes.data[currentIndex] <0)
                {
                if(latticeTypes.data[neighbors[0]]>0)
                    {
                    anchoringEnergyAtSite+= computeBoundaryEnergy(qCurrent, xDown, bounds.data[latticeTypes.data[neighbors[0]]-1]);
                    firstDerivativeX = xUp - qCurrent;
                    }
                if(latticeTypes.data[neighbors[1]]>0)
                    {
                    anchoringEnergyAtSite += computeBoundaryEnergy(qCurrent, xUp, bounds.data[latticeTypes.data[neighbors[1]]-1]);
                    firstDerivativeX = qCurrent - xDown;
                    }
                if(latticeTypes.data[neighbors[2]]>0)
                    {
                    anchoringEnergyAtSite += computeBoundaryEnergy(qCurrent, yDown, bounds.data[latticeTypes.data[neighbors[2]]-1]);
                    firstDerivativeY = yUp - qCurrent;
                    }
                if(latticeTypes.data[neighbors[3]]>0)
                    {
                    anchoringEnergyAtSite += computeBoundaryEnergy(qCurrent, yUp, bounds.data[latticeTypes.data[neighbors[3]]-1]);
                    firstDerivativeY = qCurrent - yDown;
                    }
                if(latticeTypes.data[neighbors[4]]>0)
                    {
                    anchoringEnergyAtSite += computeBoundaryEnergy(qCurrent, zDown, bounds.data[latticeTypes.data[neighbors[4]]-1]);
                    firstDerivativeZ = zUp - qCurrent;
                    }
                if(latticeTypes.data[neighbors[5]]>0)
                    {
                    anchoringEnergyAtSite += computeBoundaryEnergy(qCurrent, zUp, bounds.data[latticeTypes.data[neighbors[5]]-1]);
                    firstDerivativeZ = qCurrent - zDown;
                    }
                anchoringEnergy += anchoringEnergyAtSite;
                energyPerSite.data[i] +=anchoringEnergyAtSite;
                }
            scalar distortionEnergyAtSite=0.0;
            if(L1 != 0)
                {
                //L1
                distortionEnergyAtSite += L1*(dot(firstDerivativeX,firstDerivativeX) + firstDerivativeX[0]*firstDerivativeX[3]);
                distortionEnergyAtSite += L1*(dot(firstDerivativeY,firstDerivativeY) + firstDerivativeY[0]*firstDerivativeY[3]);
                distortionEnergyAtSite += L1*(dot(firstDerivativeZ,firstDerivativeZ) + firstDerivativeZ[0]*firstDerivativeZ[3]);
                };

            /*
            switch (numberOfConstants)
                {
                case distortionEnergyType::oneConstant :
                    {
                    break;
                    }
                case distortionEnergyType::twoConstant :
                    {
                    //L2
                    distortionEnergyAtSite += (L2*(2*firstDerivativeX[2]*firstDerivativeY[4] - 2*firstDerivativeX[2]*firstDerivativeZ[0] - 2*firstDerivativeY[4]*firstDerivativeZ[0] + 2*firstDerivativeY[1]*firstDerivativeZ[2] + 2*firstDerivativeX[0]*(firstDerivativeY[1] + firstDerivativeZ[2]) - 2*firstDerivativeX[2]*firstDerivativeZ[3] - 2*firstDerivativeY[4]*firstDerivativeZ[3] + 2*firstDerivativeZ[0]*firstDerivativeZ[3] + 2*firstDerivativeY[3]*firstDerivativeZ[4] + 2*firstDerivativeX[1]*(firstDerivativeY[3] + firstDerivativeZ[4]) + firstDerivativeX[0]*firstDerivativeX[0] + firstDerivativeX[1]*firstDerivativeX[1] + firstDerivativeX[2]*firstDerivativeX[2] + firstDerivativeY[1]*firstDerivativeY[1] + firstDerivativeY[3]*firstDerivativeY[3] + firstDerivativeY[4]*firstDerivativeY[4] + firstDerivativeZ[0]*firstDerivativeZ[0] + firstDerivativeZ[2]*firstDerivativeZ[2] + firstDerivativeZ[3]*firstDerivativeZ[3] + firstDerivativeZ[4]*firstDerivativeZ[4]))/2.;
                    //L1 with chiral
                    distortionEnergyAtSite +=(L1*(-2*firstDerivativeX[3]*firstDerivativeY[1] - 2*firstDerivativeX[4]*firstDerivativeY[2] + 2*firstDerivativeY[0]*firstDerivativeY[3] - 2*firstDerivativeX[2]*firstDerivativeZ[0] - 2*firstDerivativeX[4]*firstDerivativeZ[1] - 2*firstDerivativeY[2]*firstDerivativeZ[1] + 2*firstDerivativeX[3]*firstDerivativeZ[2] - 2*firstDerivativeY[4]*firstDerivativeZ[3] + 2*firstDerivativeY[0]*firstDerivativeZ[4] + 2*firstDerivativeY[3]*firstDerivativeZ[4] + firstDerivativeX[0]*firstDerivativeX[0] + firstDerivativeX[1]*firstDerivativeX[1] + firstDerivativeX[2]*firstDerivativeX[2] + 2*(firstDerivativeX[3]*firstDerivativeX[3]) + 2*(firstDerivativeX[4]*firstDerivativeX[4]) + 2*(firstDerivativeY[0]*firstDerivativeY[0]) + firstDerivativeY[1]*firstDerivativeY[1] + 2*(firstDerivativeY[2]*firstDerivativeY[2]) + firstDerivativeY[3]*firstDerivativeY[3] + firstDerivativeY[4]*firstDerivativeY[4] + firstDerivativeZ[0]*firstDerivativeZ[0] + 2*(firstDerivativeZ[1]*firstDerivativeZ[1]) + firstDerivativeZ[2]*firstDerivativeZ[2] + firstDerivativeZ[3]*firstDerivativeZ[3] + firstDerivativeZ[4]*firstDerivativeZ[4] + 648*(q0*q0)*(qCurrent[0]*qCurrent[0]) + 648*(q0*q0)*(qCurrent[1]*qCurrent[1]) + 648*(q0*q0)*(qCurrent[2]*qCurrent[2]) + 648*(q0*q0)*(qCurrent[3]*qCurrent[3]) + 648*(q0*q0)*(qCurrent[4]*qCurrent[4]) - 36*q0*firstDerivativeX[4]*qCurrent[0] + 72*q0*firstDerivativeY[2]*qCurrent[0] - 36*q0*firstDerivativeZ[1]*qCurrent[0] - 36*q0*firstDerivativeX[2]*qCurrent[1] + 36*q0*firstDerivativeY[4]*qCurrent[1] + 36*q0*firstDerivativeZ[0]*qCurrent[1] - 36*q0*firstDerivativeZ[3]*qCurrent[1] - 72*q0*firstDerivativeY[0]*qCurrent[2] - 36*q0*firstDerivativeY[3]*qCurrent[2] - 36*q0*firstDerivativeZ[4]*qCurrent[2] - 2*firstDerivativeX[1]*(firstDerivativeY[0] -18*q0*qCurrent[2]) - 72*q0*firstDerivativeX[4]*qCurrent[3] + 36*q0*firstDerivativeY[2]*qCurrent[3] + 36*q0*firstDerivativeZ[1]*qCurrent[3] + 648*(q0*q0)*qCurrent[0]*qCurrent[3] + 72*q0*firstDerivativeX[3]*qCurrent[4] - 36*q0*firstDerivativeY[1]*qCurrent[4] + 36*q0*firstDerivativeZ[2]*qCurrent[4] + 2*firstDerivativeX[0]*(firstDerivativeX[3] + firstDerivativeZ[2] + 18*q0*qCurrent[4])))/2.;
                    break;
                    }
                case distortionEnergyType::threeConstant :
                    {
                    //L1
                    distortionEnergyAtSite += l*(dot(firstDerivativeX,firstDerivativeX) + firstDerivativeX[0]*firstDerivativeX[3]);
                    distortionEnergyAtSite += l*(dot(firstDerivativeY,firstDerivativeY) + firstDerivativeY[0]*firstDerivativeY[3]);
                    distortionEnergyAtSite += l*(dot(firstDerivativeZ,firstDerivativeZ) + firstDerivativeZ[0]*firstDerivativeZ[3]);
                    //L2
                    distortionEnergyAtSite += (L2*(2*firstDerivativeX[2]*firstDerivativeY[4] - 2*firstDerivativeX[2]*firstDerivativeZ[0] - 2*firstDerivativeY[4]*firstDerivativeZ[0] + 2*firstDerivativeY[1]*firstDerivativeZ[2] + 2*firstDerivativeX[0]*(firstDerivativeY[1] + firstDerivativeZ[2]) - 2*firstDerivativeX[2]*firstDerivativeZ[3] - 2*firstDerivativeY[4]*firstDerivativeZ[3] + 2*firstDerivativeZ[0]*firstDerivativeZ[3] + 2*firstDerivativeY[3]*firstDerivativeZ[4] + 2*firstDerivativeX[1]*(firstDerivativeY[3] + firstDerivativeZ[4]) + firstDerivativeX[0]*firstDerivativeX[0] + firstDerivativeX[1]*firstDerivativeX[1] + firstDerivativeX[2]*firstDerivativeX[2] + firstDerivativeY[1]*firstDerivativeY[1] + firstDerivativeY[3]*firstDerivativeY[3] + firstDerivativeY[4]*firstDerivativeY[4] + firstDerivativeZ[0]*firstDerivativeZ[0] + firstDerivativeZ[2]*firstDerivativeZ[2] + firstDerivativeZ[3]*firstDerivativeZ[3] + firstDerivativeZ[4]*firstDerivativeZ[4]))/2.;
                    //L3
                    distortionEnergyAtSite += L3*((firstDerivativeX[0]*firstDerivativeX[3] - firstDerivativeZ[0]*firstDerivativeZ[3] + firstDerivativeX[0]*firstDerivativeX[0] + firstDerivativeX[1]*firstDerivativeX[1] + firstDerivativeX[2]*firstDerivativeX[2] + firstDerivativeX[3]*firstDerivativeX[3] + firstDerivativeX[4]*firstDerivativeX[4] - firstDerivativeZ[0]*firstDerivativeZ[0] - firstDerivativeZ[1]*firstDerivativeZ[1] - firstDerivativeZ[2]*firstDerivativeZ[2] - firstDerivativeZ[3]*firstDerivativeZ[3] - firstDerivativeZ[4]*firstDerivativeZ[4])*qCurrent[0] + 2*firstDerivativeX[0]*firstDerivativeY[0]*qCurrent[1] + firstDerivativeX[3]*firstDerivativeY[0]*qCurrent[1] + 2*firstDerivativeX[1]*firstDerivativeY[1]*qCurrent[1] + 2*firstDerivativeX[2]*firstDerivativeY[2]*qCurrent[1] + firstDerivativeX[0]*firstDerivativeY[3]*qCurrent[1] + 2*firstDerivativeX[3]*firstDerivativeY[3]*qCurrent[1] + 2*firstDerivativeX[4]*firstDerivativeY[4]*qCurrent[1] + 2*firstDerivativeX[0]*firstDerivativeZ[0]*qCurrent[2] + firstDerivativeX[3]*firstDerivativeZ[0]*qCurrent[2] + 2*firstDerivativeX[1]*firstDerivativeZ[1]*qCurrent[2] + 2*firstDerivativeX[2]*firstDerivativeZ[2]*qCurrent[2] + firstDerivativeX[0]*firstDerivativeZ[3]*qCurrent[2] + 2*firstDerivativeX[3]*firstDerivativeZ[3]*qCurrent[2] + 2*firstDerivativeX[4]*firstDerivativeZ[4]*qCurrent[2] + (firstDerivativeY[0]*firstDerivativeY[3] - firstDerivativeZ[0]*firstDerivativeZ[3] + firstDerivativeY[0]*firstDerivativeY[0] + firstDerivativeY[1]*firstDerivativeY[1] + firstDerivativeY[2]*firstDerivativeY[2] + firstDerivativeY[3]*firstDerivativeY[3] + firstDerivativeY[4]*firstDerivativeY[4] - firstDerivativeZ[0]*firstDerivativeZ[0] - firstDerivativeZ[1]*firstDerivativeZ[1] - firstDerivativeZ[2]*firstDerivativeZ[2] - firstDerivativeZ[3]*firstDerivativeZ[3] - firstDerivativeZ[4]*firstDerivativeZ[4])*qCurrent[3] + 2*firstDerivativeY[0]*firstDerivativeZ[0]*qCurrent[4] + firstDerivativeY[3]*firstDerivativeZ[0]*qCurrent[4] + 2*firstDerivativeY[1]*firstDerivativeZ[1]*qCurrent[4] + 2*firstDerivativeY[2]*firstDerivativeZ[2]*qCurrent[4] + firstDerivativeY[0]*firstDerivativeZ[3]*qCurrent[4] + 2*firstDerivativeY[3]*firstDerivativeZ[3]*qCurrent[4] + 2*firstDerivativeY[4]*firstDerivativeZ[4]*qCurrent[4]);
                    break;
                    }

                };
                */
                /*
            if(useL24)
                {
                distortionEnergyAtSite += 3*L24*(firstDerivativeX[1]*firstDerivativeY[0] - firstDerivativeX[0]*firstDerivativeY[1] + firstDerivativeX[3]*firstDerivativeY[1] + firstDerivativeX[4]*firstDerivativeY[2] - firstDerivativeX[1]*firstDerivativeY[3] - firstDerivativeX[2]*firstDerivativeY[4] + 2*firstDerivativeX[2]*firstDerivativeZ[0] + firstDerivativeY[4]*firstDerivativeZ[0] + firstDerivativeX[4]*firstDerivativeZ[1] + firstDerivativeY[2]*firstDerivativeZ[1] - (2*firstDerivativeX[0] + firstDerivativeX[3] + firstDerivativeY[1])*firstDerivativeZ[2] + firstDerivativeX[2]*firstDerivativeZ[3] + 2*firstDerivativeY[4]*firstDerivativeZ[3] - (firstDerivativeX[1] + firstDerivativeY[0] + 2*firstDerivativeY[3])*firstDerivativeZ[4]);
                }
                */
            distortionEnergy+=distortionEnergyAtSite;
            energyPerSite.data[i] +=distortionEnergyAtSite;
            }

        };
    energy = (phaseEnergy + distortionEnergy + anchoringEnergy + eFieldEnergy + hFieldEnergy);
    energyComponents[0] = phaseEnergy;
    energyComponents[1] = distortionEnergy;
    energyComponents[2] = anchoringEnergy;
    energyComponents[3] = eFieldEnergy;
    energyComponents[4] = hFieldEnergy;

//    if(verbose)
//        printf("%f %f %f %f %f\n",phaseEnergy , distortionEnergy , anchoringEnergy , eFieldEnergy , hFieldEnergy);
    };
