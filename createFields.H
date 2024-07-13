    Info<< "Reading transportProperties\n" << endl;
      IOdictionary fluidTransportProperties
    (
        IOobject
        (
            "fluidTransportProperties",
            runTime.constant(),
            fluidmesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
      IOdictionary solidTransportProperties
    (
        IOobject
        (
            "solidTransportProperties",
            runTime.constant(),
            solidmesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    ////fluid part  
    // Density [kg/m^3]
    //dimensionedScalar rhof(fluidTransportProperties.lookup("rhof"));
    dimensionedScalar rhof
    (
    "rhof",
    fluidTransportProperties.lookup("rhof")
    );
    // Heat capacity [J/(kg*K)]
    //dimensionedScalar Cpf(fluidTransportProperties.lookup("Cpf"));
    dimensionedScalar Cpf
    (
    "Cpf",
    fluidTransportProperties.lookup("Cpf")
    );
    // Thermal conduction [W/(m*K)]
    //dimensionedScalar Kf(fluidTransportProperties.lookup("Kf"));    
    dimensionedScalar Kf
    (
    "Kf",
    fluidTransportProperties.lookup("Kf")
    );
    //kenematic viscosity [m^2/s]
   // dimensionedScalar nu(fluidTransportProperties.lookup("nu"));
    dimensionedScalar nu
    (
    "nu",
    fluidTransportProperties.lookup("nu")
    );
    
    ////solid part   
    // Density [kg/m^3]
    dimensionedScalar rhos(solidTransportProperties.lookup("rhos"));

    // Heat capacity [J/(kg*K)]
    dimensionedScalar Cps(solidTransportProperties.lookup("Cps"));
    
    // Thermal conduction [W/(m*K)]
    dimensionedScalar Ks(solidTransportProperties.lookup("Ks"));  
    
    //for coupled interface
    dictionary  coupledDict = 
    IOdictionary
    (
        IOobject
        (
            "coupledInterface",
            runTime.constant(),
            fluidmesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );// located in file "constant"

    label fluidPatchIndex = -1;
    label solidPatchIndex = -1;
    {
      word fluidPatchName(coupledDict.lookup("fluidSideCoupledPatchName"));
      fluidPatchIndex = fluidmesh.boundaryMesh().findPatchID(fluidPatchName);

      word solidPatchName(coupledDict.lookup("solidSideCoupledPatchName"));
      solidPatchIndex = solidmesh.boundaryMesh().findPatchID(solidPatchName);
    }

    scalar tolerance = coupledDict.lookupOrDefault<scalar>("tolerance", 1e-4);
    label iterationN = coupledDict.lookupOrDefault<label>("numberOfIteration", 50);
    
    #include "createFluidFields.H"
    #include "createSolidFields.H"   
