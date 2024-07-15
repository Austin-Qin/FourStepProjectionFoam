/*---------------------------------------------------------------------------------------------*\                        
Using so-called "D-D" iteration to solve temperature field of all region, wall temperature(Tw) of one region is calculated by the equality of heat flux between two regions.

\*---------------------------------------------------------------------------------------------*/
//// "D-D" iteration ////

 label counter = 0;
 scalar errorT = 0;

        scalarField boundaryFluidTemp(T, fluidFaceCells);
        scalarField boundarySolidTemp(TS, solidFaceCells);  //these are the temperature of cells which next to the coupled face 
        
        fvScalarMatrix TEqn
            (
                 fvm::laplacian(Kf, T) == rhof*Cpf*(fvm::div(phi, T) + fvc::ddt(T) ) 
                              
            );

       fvScalarMatrix TSEqn
            (
                fvm::laplacian(Ks, TS) == rhos*Cps*fvc::ddt(TS)
            );
            
       const fvPatchField<scalar>& psfFluid = T.boundaryField()[fluidPatchIndex];
       const fvPatchField<scalar>& psfSolid = TS.boundaryField()[solidPatchIndex];
    do
    {      
       // coupled boundary in fluid  
       scalarField Tw(nFace,0.0);
       forAll(fluidToSolid,fI)
       {
        Tw[fI] =  (fluidKD[fI]*boundaryFluidTemp[fI] + solidKD[fluidToSolid[fI]]*boundarySolidTemp[fluidToSolid[fI]])/(fluidKD[fI] + solidKD[fluidToSolid[fI]]);
       }// boundary temperature feild from Nuemann B.C.
       scalarField& TB = T.boundaryFieldRef()[fluidPatchIndex];
       TB = Tw;
       T.correctBoundaryConditions();
       TEqn.boundaryCoeffs()[fluidPatchIndex] = -patchKFluid * psfFluid.gradientBoundaryCoeffs(); 

       TEqn.solve();
       
         forAll(boundaryFluidTemp, cI)
         {
          boundaryFluidTemp[cI] = T.internalField()[fluidFaceCells[cI]];
         } //renew fluid temperature of the boundary cells
       

       // coupled boundary in solid
       scalarField TSw(nFace,0.0);
       forAll(solidToFluid,fI)
       {
        TSw[fI] =  (fluidKD[solidToFluid[fI]]*boundaryFluidTemp[solidToFluid[fI]] + solidKD[fI]*boundarySolidTemp[fI])/(fluidKD[solidToFluid[fI]] + solidKD[fI]);
       }// boundary temperature feild from Nuemann B.C.
       scalarField& TSB = TS.boundaryFieldRef()[solidPatchIndex];
       TSB = TSw;
       TS.correctBoundaryConditions();
       
       TSEqn.boundaryCoeffs()[solidPatchIndex] = -patchKSolid * psfSolid.gradientBoundaryCoeffs();
       TSEqn.solve();
       
        forAll(boundarySolidTemp, cI)
         {
          boundarySolidTemp[cI] = TS.internalField()[solidFaceCells[cI]];
         }//renew solid temperature of the boundary cells
         
       // iteration creterion
       errorT = max
                (
                    mag(Tw-TSw)/max(max(mag(boundaryFluidTemp-boundarySolidTemp)),SMALL)
                );
        counter++;


    }
    while(errorT > tolerance && counter < iterationN);
    
    Info<< "\nloop T counter = " << counter << endl;
    Info<< "\nloop T errorT = " << errorT << endl;
    
    ////finish////
        




