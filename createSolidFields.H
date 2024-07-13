Info << "Reading field solid T\n" << endl;
volScalarField TS
(
    IOobject
    (
        "TS",
        runTime.timeName(),
        solidmesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    solidmesh
);
