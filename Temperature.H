

 label counter = 0;
 scalar errorT = 0;

 scalarField boundaryFluidTemp(T, fluidFaceCells);
 scalarField boundarySolidTemp(TS, solidFaceCells);  
        fvScalarMatrix TEqn
            (
                rhof*Cpf*(fvm::ddt(T) + fvm::div(phi, T)) == fvm::laplacian(Kf, T)
                              
            );

       fvScalarMatrix TSEqn
            (
                rhos*Cps*fvm::ddt(TS) == fvm::laplacian(Ks, TS)


            );
            

    do
    {      
       scalarField Tw(nFace,0.0);
       forAll(fluidToSolid,fI)
       {
        Tw[fI] =  (fluidKD[fI]*boundaryFluidTemp[fI] + solidKD[fluidToSolid[fI]]*boundarySolidTemp[fluidToSolid[fI]])/(fluidKD[fI] + solidKD[fluidToSolid[fI]]);
       }
       scalarField& TB = T.boundaryFieldRef()[fluidPatchIndex];
       TB = Tw;
       T.correctBoundaryConditions();
       
       
       TEqn.solve();

       
       scalarField TSw(nFace,0.0);
       forAll(solidToFluid,fI)
       {
        TSw[fI] =  (fluidKD[solidToFluid[fI]]*boundaryFluidTemp[solidToFluid[fI]] + solidKD[fI]*boundarySolidTemp[fI])/(fluidKD[solidToFluid[fI]] + solidKD[fI]);
       }
       scalarField& TSB = TS.boundaryFieldRef()[solidPatchIndex];
       TSB = TSw;
       TS.correctBoundaryConditions();
       
       
       TSEqn.solve();
       
       errorT = max
                (
                    mag(Tw-TSw)/max(max(mag(boundaryFluidTemp-boundarySolidTemp)),SMALL)
                );
        counter++;


    }
    while(errorT > tolerance && counter < iterationN);
    
    Info<< "\nloop T counter = " << counter << endl;
    Info<< "\nloop T errorT = " << errorT << endl;
    

        




