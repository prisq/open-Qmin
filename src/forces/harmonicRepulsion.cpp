#include "harmonicRepulsion.h"
#include "harmonicRepulsion.cuh"
#include "utilities.cuh"
/*! \file harmonicRepulsion.cpp */
/*!
Force defined by....
*/
/*!
\pre the vector has n^2 elements, where n is the number of types, and the values of type in the system
must be of the form 0, 1, 2, ...n. The input vector should be laid out as:
params[0] = k_{0,0}
params[1] = k_{0,1}
params[n] = k_{0,n}
params[n+1] = k_{1,0} (physically, this better be the same as g_{0,1})
params[n+2] = k_{1,1}
...
params[n^2-1] = k_{n,n}
*/
void harmonicRepulsion::setForceParameters(vector<scalar> &params)
    {
    parameters.resize(params.size());
    nTypes = sqrt(params.size());
    particleTypeIndexer = Index2D(nTypes);
    ArrayHandle<scalar> h_p(parameters);
    for(int ii = 0; ii < params.size(); ++ii)
        {
        int typeI = ii/nTypes;
        int typeJ = ii - typeI*nTypes;
        h_p.data[particleTypeIndexer(typeJ,typeI)] = params[ii];
        };
    };

/*
Need to get the type-type stiffness and the sum of radii
*/
void harmonicRepulsion::getParametersForParticlePair(int index1, int index2, vector<scalar> &params)
    {
    ArrayHandle<int> particleType(model->returnTypes());
    ArrayHandle<scalar> h_p(parameters);
    ArrayHandle<scalar> h_r(model->returnRadii());
    params[0] = h_p.data[particleTypeIndexer(particleType.data[index2],particleType.data[index1])];
    params[1] = (h_r.data[index1]+h_r.data[index2]);
    };

void harmonicRepulsion::computePairwiseForce(dVec &relativeDistance, vector<scalar> &params, dVec &f)
    {
    scalar dnorm = norm(relativeDistance);
    //scalar delta = (1.0 - dnorm/params[1]);
    if(dnorm <= params[1])
        f=  params[0]*(1.0/params[1])*(1.0 - dnorm/params[1])*(1.0/dnorm)*relativeDistance;
    else
        f = make_dVec(0.0);
    };

void harmonicRepulsion::computeForcesGPU(GPUArray<dVec> &forces, bool zeroOutForce)
    {
    if(!useNeighborList)
        UNWRITTENCODE("pairwise forces without neighbor list not written yet");

    ArrayHandle<dVec> d_force(forces,access_location::device,access_mode::readwrite);
    int N = model->getNumberOfParticles();
    if(zeroOutForce)
        gpu_zero_array(d_force.data,N);

    neighbors->computeNeighborLists(model->returnPositions());
    ArrayHandle<unsigned int> d_npp(neighbors->neighborsPerParticle,access_location::device,access_mode::read);
    ArrayHandle<int> d_n(neighbors->particleIndices,access_location::device,access_mode::read);
    ArrayHandle<dVec> d_nv(neighbors->neighborVectors,access_location::device,access_mode::read);

    ArrayHandle<int> particleType(model->returnTypes(),access_location::device,access_mode::read);
    ArrayHandle<scalar> d_params(parameters,access_location::device,access_mode::read);
    ArrayHandle<scalar> d_radii(model->returnRadii(),access_location::device,access_mode::read);
    gpu_harmonic_repulsion_calculation(d_force.data,
                                       d_npp.data,
                                       d_n.data,
                                       d_nv.data,
                                       particleType.data,
                                       d_radii.data,
                                       d_params.data,
                                       neighbors->neighborIndexer,
                                       particleTypeIndexer,
                                       N);
    };
