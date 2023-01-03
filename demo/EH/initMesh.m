function [PHTelem,controlPts,meshInfo] = initMesh(geometry)

meshInfo = struct;
numPatches = geometry.numPatches;
p = geometry.p;
q = geometry.q;
nurbsList = cell(1,numPatches);
controlPts = cell(numPatches,1);
PHTelem = cell(numPatches,1);
quadList = cell(numPatches,1);
numElemU = geometry.numElemU;
numElemV = geometry.numElemV;
dimBasis = zeros(1, numPatches);
numElements_patch = zeros(1,numPatches);

% patchCordsTri = [PatchIndex for the 1st Quad in the square, PatchXMin, PatchXMax, PatchYMin, PatchYMax, X-Centre of the circle, Y-Centre of the circle]
patchCordsTri = [1,0.0,1.0,0.0,1.0,0.6,0.5];
numberElements = 0;
for iPatch = 1:size(patchCordsTri,1)
    
    patchMinX = patchCordsTri(iPatch,2);
    patchMaxX = patchCordsTri(iPatch,3);
    patchMinY = patchCordsTri(iPatch,4);
    patchMaxY = patchCordsTri(iPatch,5);
    centreX = patchCordsTri(iPatch,6);
    centreY = patchCordsTri(iPatch,7);
    rad = 0.2;
    
    patchIndex = patchCordsTri(iPatch,1);
    startAngle = pi/4;
    endAngle = 3*pi/4;
    arcSide = nrbcirc(rad, [centreX, centreY, 0], startAngle, endAngle);
    arcSide.coefs = fliplr(arcSide.coefs);
    startPtX = arcSide.coefs(1,1);
    startPtY = arcSide.coefs(2,1);
    side1 = nrbline([startPtX, startPtY],[patchMinX,patchMaxY]);
    side2 = nrbline([patchMinX, patchMaxY],[patchMaxX, patchMaxY]);
    endPtX = arcSide.coefs(1,end);
    endPtY = arcSide.coefs(2,end);
    side3 = nrbline([endPtX, endPtY],[patchMaxX, patchMaxY]);
    nurbs = nrbcoons(arcSide,side2,side1,side3);
    nurbsList{patchIndex} = nurbs;
    [controlPts{patchIndex},PHTelem{patchIndex},dimBasis(patchIndex),quadList{patchIndex}] = initElements(nurbs,p,q,numElemU,numElemV);
    numberElements = numberElements + 4*size(quadList{patchIndex},1);
    
    patchIndex = patchCordsTri(iPatch,1)+ 1;
    startAngle = 3*pi/4;
    endAngle = 5*pi/4;
    arcSide = nrbcirc(rad, [centreX, centreY, 0], startAngle, endAngle);
    arcSide.coefs = fliplr(arcSide.coefs);
    startPtX = arcSide.coefs(1,1);
    startPtY = arcSide.coefs(2,1);
    side1 = nrbline([startPtX, startPtY],[patchMinX,patchMinY]);
    side2 = nrbline([patchMinX,patchMinY], [patchMinX, patchMaxY]);
    endPtX = arcSide.coefs(1,end);
    endPtY = arcSide.coefs(2,end);
    side3 = nrbline([endPtX, endPtY],[patchMinX, patchMaxY]);
    nurbs = nrbcoons(arcSide,side2,side1,side3);
    nurbsList{patchIndex} = nurbs;
    [controlPts{patchIndex},PHTelem{patchIndex},dimBasis(patchIndex),quadList{patchIndex}] = initElements(nurbs,p,q,numElemU,numElemV);
    numberElements = numberElements + 4*size(quadList{patchIndex},1);
    
    patchIndex = patchCordsTri(iPatch,1)+ 2;
    startAngle = 5*pi/4;
    endAngle = 7*pi/4;
    arcSide = nrbcirc(rad, [centreX, centreY, 0], startAngle, endAngle);
    arcSide.coefs = fliplr(arcSide.coefs);
    startPtX = arcSide.coefs(1,1);
    startPtY = arcSide.coefs(2,1);
    side1 = nrbline([startPtX, startPtY],[patchMaxX,patchMinY]);
    side2 = nrbline([patchMaxX,patchMinY] , [patchMinX, patchMinY]);
    endPtX = arcSide.coefs(1,end);
    endPtY = arcSide.coefs(2,end);
    side3 = nrbline([endPtX, endPtY],[patchMinX, patchMinY]);
    nurbs = nrbcoons(arcSide,side2,side1,side3);
    nurbsList{patchIndex} = nurbs;
    [controlPts{patchIndex},PHTelem{patchIndex},dimBasis(patchIndex),quadList{patchIndex}] = initElements(nurbs,p,q,numElemU,numElemV);
    numberElements = numberElements + 4*size(quadList{patchIndex},1);
    
    patchIndex = patchCordsTri(iPatch,1)+ 3;
    startAngle = 7*pi/4;
    endAngle = pi/4;
    arcSide = nrbcirc(rad, [centreX, centreY, 0], startAngle, endAngle);
    arcSide.coefs = fliplr(arcSide.coefs);
    startPtX = arcSide.coefs(1,1);
    startPtY = arcSide.coefs(2,1);
    side1 = nrbline([startPtX, startPtY],[patchMaxX,patchMaxY]);
    side2 = nrbline([patchMaxX,patchMaxY], [patchMaxX, patchMinY]);
    endPtX = arcSide.coefs(1,end);
    endPtY = arcSide.coefs(2,end);
    side3 = nrbline([endPtX, endPtY],[patchMaxX, patchMinY]);
    nurbs = nrbcoons(arcSide,side2,side1,side3);
    nurbs = nrbreverse(nurbs,[1]);
    nurbsList{patchIndex} = nurbs;
    [controlPts{patchIndex},PHTelem{patchIndex},dimBasis(patchIndex),quadList{patchIndex}] = initElements(nurbs,p,q,numElemU,numElemV);
    numberElements = numberElements + 4*size(quadList{patchIndex},1);
end

meshInfo.dimBasis = dimBasis;
meshInfo.quadList = quadList;
meshInfo.numElements = numberElements;
meshInfo.numElement_patch = numElements_patch;
end