/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    FourStepProjection

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "pisoControl.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    //pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    #include "preparation.H"
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"
        #include "continuityErrs.H"



        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + 1.5 * fvc::div(phi,U) - 0.5 * fvc::div(phi.oldTime(),U.oldTime())
          - 0.5 * fvm::laplacian(nu, U) - 0.5 * fvc::laplacian(nu,U) 
        );

            solve(UEqn == -fvc::grad(p));
        
        
            U = U + fluidmesh.time().deltaT() * fvc::grad(p); 
            U.correctBoundaryConditions();

        surfaceScalarField Uphi
        (
            "Uphi",
            fvc::interpolate(U)& fluidmesh.Sf()
        );

            fvScalarMatrix pEqn
            (
            fvm::laplacian(fluidmesh.time().deltaT(),p) == fvc::div(Uphi)
            );
            pEqn.setReference(pRefCell, pRefValue);
            solve(pEqn);

            phi = Uphi - pEqn.flux();




            U = U - fluidmesh.time().deltaT() * fvc::grad(p);
            U.correctBoundaryConditions();


        
        #include "Temperature.H"

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;
return 0;
}