function [PHTelem,controlPts,meshInfo] = initMesh(geometry)

meshInfo = struct;
L = geometry.L;
W_left = geometry.W_left;
W_right = geometry.W_right;
p = geometry.p;
q = geometry.q;
numberElemU = geometry.numElemU;
numberElemV = geometry.numElemV;
numPatches = geometry.numPatches;
numElements_patch = zeros(1,geometry.numPatches);
numberElements = 0;
% Initialize the PHT geometry on coarse mesh
controlPts = cell(numPatches,1);
PHTelem = cell(numPatches, 1);
dimBasis = zeros(1, numPatches);
quadList = cell(numPatches,1);

patchIndex = 1;
% Initialize geometry on coarsest mesh
coefs(1:3,1,1) = [0; -W_left/2; 0];
coefs(1:3,1,2) = [0; 0; 0];
coefs(1:3,2,1) = [L; -W_right/2; 0];
coefs(1:3,2,2) = [L; 0;0];
coefs(4,1,1) = 1;
coefs(4,1,2) = 1;
coefs(4,2,1) = 1;
coefs(4,2,2) = 1;

knotU = [0 0 1 1];
knotV = [0 0 1 1];

nurbs = nrbmak(coefs,{knotU,knotV});
p_init = nurbs.order(1)-1;
q_init = nurbs.order(2)-1;

% Refine into numberElemU by numberElemV knotspans
knotU = linspace(0,1,numberElemU+1);
knotV = linspace(0,1,numberElemV+1);

numberElementsU = length(unique(knotU))-1;
numberElementsV = length(unique(knotV))-1;

% Increase polynomial order
nurbs = nrbdegelev(nurbs,[p-p_init,q-q_init]);
nurbs = nrbkntins(nurbs,{knotU(2:end-1) knotV(2:end-1)});
% Repeat the knots to get C1 continuity
nurbs = nrbkntins(nurbs,{knotU(2:end-1) knotV(2:end-1)});

[controlPts{patchIndex},PHTelem{patchIndex},dimBasis(patchIndex),quadList{patchIndex}] = genControlPtsNoRep(nurbs,p,q,numberElementsU,numberElementsV);

% Refine the area near the support
[quadRef] = refineMeshCoordinates(PHTelem{patchIndex},controlPts{patchIndex},quadList{patchIndex},p,q,0.0,1.0-0.01,(-W_left/2),0.0);
[quadList{patchIndex},PHTelem{patchIndex}, controlPts{patchIndex},dimBasis(patchIndex)] = refineMeshIso(quadRef,quadList{patchIndex},PHTelem{patchIndex},controlPts{patchIndex},p,q,dimBasis(patchIndex));

% Refine the area near the support
[quadRef] = refineMeshCoordinates(PHTelem{patchIndex},controlPts{patchIndex},quadList{patchIndex},p,q,0.0,0.5-0.01,(-W_left/2),0.0);
[quadList{patchIndex},PHTelem{patchIndex}, controlPts{patchIndex},dimBasis(patchIndex)] = refineMeshIso(quadRef,quadList{patchIndex},PHTelem{patchIndex},controlPts{patchIndex},p,q,dimBasis(patchIndex));

numberElements = numberElements + 4*size(quadList{patchIndex},1);
numElements_patch(patchIndex) = 4*size(quadList{patchIndex},1);

patchIndex = 2;
% Initialize geometry on coarsest mesh
coefs(1:3,1,1) = [0; 0; 0];
coefs(1:3,1,2) = [0; W_left/2; 0];
coefs(1:3,2,1) = [L; 0; 0];
coefs(1:3,2,2) = [L; W_right/2;0];
coefs(4,1,1) = 1;
coefs(4,1,2) = 1;
coefs(4,2,1) = 1;
coefs(4,2,2) = 1;

knotU = [0 0 1 1];
knotV = [0 0 1 1];

nurbs = nrbmak(coefs,{knotU,knotV});
p_init = nurbs.order(1)-1;
q_init = nurbs.order(2)-1;

% Refine into numberElemU by numberElemV knotspans
knotU = linspace(0,1,numberElemU+1);
knotV = linspace(0,1,numberElemV+1);

numberElementsU = length(unique(knotU))-1;
numberElementsV = length(unique(knotV))-1;

% Increase polynomial order
nurbs = nrbdegelev(nurbs,[p-p_init,q-q_init]);
nurbs = nrbkntins(nurbs,{knotU(2:end-1) knotV(2:end-1)});
% Repeat the knots to get C1 continuity
nurbs = nrbkntins(nurbs,{knotU(2:end-1) knotV(2:end-1)});

[controlPts{patchIndex},PHTelem{patchIndex},dimBasis(patchIndex),quadList{patchIndex}] = genControlPtsNoRep(nurbs,p,q,numberElementsU,numberElementsV);

% Refine the area near the support
[quadRef] = refineMeshCoordinates(PHTelem{patchIndex},controlPts{patchIndex},quadList{patchIndex},p,q,0.0,1.0-0.01,0.0,(W_left/2));
[quadList{patchIndex},PHTelem{patchIndex}, controlPts{patchIndex},dimBasis(patchIndex)] = refineMeshIso(quadRef,quadList{patchIndex},PHTelem{patchIndex},controlPts{patchIndex},p,q,dimBasis(patchIndex));

% Refine the area near the support
[quadRef] = refineMeshCoordinates(PHTelem{patchIndex},controlPts{patchIndex},quadList{patchIndex},p,q,0.0,0.5-0.01,0.0,(W_left/2));
[quadList{patchIndex},PHTelem{patchIndex}, controlPts{patchIndex},dimBasis(patchIndex)] = refineMeshIso(quadRef,quadList{patchIndex},PHTelem{patchIndex},controlPts{patchIndex},p,q,dimBasis(patchIndex));

numberElements = numberElements + 4*size(quadList{patchIndex},1);
numElements_patch(patchIndex) = 4*size(quadList{patchIndex},1);

meshInfo.dimBasis = dimBasis;
meshInfo.quadList = quadList;
meshInfo.numElements = numberElements;
meshInfo.numElement_patch = numElements_patch;

end