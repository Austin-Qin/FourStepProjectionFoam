    Info<< "Reading fluid field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            fluidmesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fluidmesh
    );

    Info<< "Reading fluid field T\n" << endl;    
      volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            fluidmesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fluidmesh
    );


    Info<< "Reading fluid field U\n" << endl;
      volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            fluidmesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fluidmesh
    );


     surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        fluidmesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    linearInterpolate(U) & fluidmesh.Sf()
);


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, fluidmesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
    
