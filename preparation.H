

const unallocLabelList& fluidFaceCells = fluidmesh.boundaryMesh()[fluidPatchIndex].faceCells();
const unallocLabelList& solidFaceCells = solidmesh.boundaryMesh()[solidPatchIndex].faceCells();

label nFace = fluidFaceCells.size();

const scalarField patchKFluid = Kf.value()*fluidmesh.magSf().boundaryField()[fluidPatchIndex];
const scalarField patchKSolid = Ks.value()*solidmesh.magSf().boundaryField()[solidPatchIndex];


const vectorField& faceCenter_F = fluidmesh.Cf().boundaryField()[fluidPatchIndex];
const vectorField& faceCenter_S = solidmesh.Cf().boundaryField()[solidPatchIndex];

labelList fluidToSolid(nFace,-1);
labelList solidToFluid(nFace,-1);

forAll(fluidToSolid,fI)
{
 vector fCenter  = faceCenter_F[fI];
 scalar sup      = GREAT;
 label fNumber   = -1;
  forAll(faceCenter_S,cI)
  {
   scalar distance = mag(faceCenter_S[cI]-fCenter);
   if(distance < SMALL)
   {
     fNumber = cI;
     break;
   }
   else if (distance < sup )
   {
     sup  = distance;
     fNumber =   cI; 
   }
  }
  fluidToSolid[fI] = fNumber; 
}

forAll(solidToFluid,fI)
{
  solidToFluid[fluidToSolid[fI]] = fI;
}



  const scalarField& fluidDelta = fluidmesh.deltaCoeffs().boundaryField()[fluidPatchIndex];
  const scalarField& solidDelta = solidmesh.deltaCoeffs().boundaryField()[solidPatchIndex];
  
  scalarField fluidKD = Kf.value() * fluidDelta ;
  scalarField solidKD = Ks.value() * solidDelta ;


