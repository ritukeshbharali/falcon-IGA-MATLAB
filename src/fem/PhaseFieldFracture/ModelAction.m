%-------------------------------------------------------------------------%
% ModelAction is a function that carries out an action requested by the 
% solver. Here is a list of actions for the Phase-Field Fracture FE Model:
%    - INIT
%    - GET_MATRIX0
%    - GET_MATRIX0_EXT
%    - GET_INT_FORCE
%    - GET_MATRIX0_1
%    - GET_MATRIX0_2
%    - GET_DISSIPATION
%    - GET_AT2_DISSIPATION
%    - GET_MATRIX2
%    - COMMIT
%    - REVERT
%    - REFINE_N_TRANSFER
%-------------------------------------------------------------------------%

function globdat = ModelAction(props,globdat,action)

switch action

    case 'INIT'

        globdat = init(props,globdat);

    case 'GET_MATRIX0'
        
        [globdat.K, globdat.fint, globdat.markRef] = ...
            getMatrix0(props,globdat);

    case 'GET_MATRIX0_EXT'

        [globdat.K, globdat.fint, globdat.markRef] = ...
            getMatrix0_Ext(props,globdat);     

    case 'GET_INT_FORCE'
        
        [globdat.fint, globdat.markRef] = ...
            getIntForce(props,globdat);     
    
    case 'GET_MATRIX0_1'
        
        [globdat.K1, globdat.fint1, globdat.markRef] = ...
            getMatrix0_1(props,globdat);
        
    case 'GET_MATRIX0_2'
        
        [globdat.K2, globdat.fint2] = ...
            getMatrix0_2(props,globdat);

    case 'GET_DISSIPATION'
        
        [globdat.h, globdat.g] = ...
            getDissipation(props,globdat);

    case 'GET_AT2_DISSIPATION'
        
        [globdat.h, globdat.g] = ...
            getAT2Dissipation(props,globdat);    
        
    case 'GET_MATRIX20'
        
        globdat.M = ...
            getMatrix20(props,globdat);
        
    case 'GET_MATRIX21'
        
        globdat.M = ...
            getMatrix21(props,globdat);    
        
    case 'COMMIT'    
        
        globdat.state00 = globdat.state0;
        globdat.state0  = globdat.state;

    case 'REVERT'    
        
        globdat.state   = globdat.state0;
        globdat.state0  = globdat.state00;    

    case 'REFINE_N_TRANSFER'

        % Refine and solution transfer
        globdat = refineNTransferSolution(props,globdat);

        % Recompute constraint matrix
        globdat = getConstraints(props,globdat);        

end

end


%=========================================================================%
% INIT  
%=========================================================================%

function globdat = init(props,globdat)

% Compute problem dimension and total number of system dofs
globdat.dim   = props.geom.dim;
globdat.ndofs = (globdat.dim+1) * globdat.mesh.sizeBasis; 

% Check phase-field model type
switch props.mat.eta
    case 1
        disp('   + Fracture Type: Brittle AT1')
    case 0
        disp('   + Fracture Type: Brittle AT2')
    case 2
        disp('   + Fracture Type: PFCZM')
    otherwise
        error('Choose eta: 0, 1 or 2 for AT2, AT1 and PFCZM')
end

% Pvol 
Pvol               = zeros(6);
Pidx               = [1 2 3];
Pvol(Pidx,Pidx)    = 1;

% Pdev
Pdev               = zeros(6);
Pdev(1:3,1:3)      = -1/3;
for i = 1:3
    Pdev(i,i)      = 2/3;
end
for i = 4:6
    Pdev(i,i)      = 1/2;
end

if globdat.dim == 2
    Pidx = [3 5 6];
    Pvol(Pidx,:) = [];
    Pvol(:,Pidx) = [];
    Pdev(Pidx,:) = [];
    Pdev(:,Pidx) = [];
    globdat.Pvol = Pvol;
    globdat.Pdev = Pdev;
else
    globdat.Pvol = Pvol;
    globdat.Pdev = Pdev;
end

end


%=========================================================================%
% GET_MATRIX0
%=========================================================================%

function [K, fint, markRef] = getMatrix0(props,globdat)

% Some useful variables
material   = props.mat;
dim        = props.geom.dim;
nstress    = props.geom.nstress;
numPatches = props.geom.numPatches;
ngaussX    = props.geom.ngaussX;
ngaussY    = props.geom.ngaussY;
ndofs      = globdat.ndofs;
udofs      = dim*globdat.mesh.sizeBasis;

% Initialise cell for mesh refinement marker
markRef    = cell(numPatches,1);

% Initialise vector for global residual
fint        = zeros(ndofs,1);

% Initialise cells for local stiffness matrices
nelem      = globdat.mesh.numElements;
kUU        = cell(nelem,1);
kUPhi      = cell(nelem,1);
kPhiU      = cell(nelem,1);
kPhiPhi    = cell(nelem,1);

% Set index counter to zero
elementCounter = 0;
indexCounter   = 0;

% Loop over patches
for indexPatch = 1:numPatches
    
    % Initialise vectors for elements in a patch
    markRef(indexPatch) = {zeros(1,length(globdat.PHTelem{indexPatch}))};
    
    % Loop over elements
    for i=1:length(globdat.PHTelem{indexPatch})
        if isempty(globdat.PHTelem{indexPatch}(i).children)
            
            % Extract element displacement and phase-field
            nument       = size(globdat.PHTelem{indexPatch}(i).C,1);
            elemNodes    = globdat.PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            elemUDofs    = reshape([2*elemNodes-1;2*elemNodes],1,dim*nument);
            elemPhiDofs  = elemNodes+udofs;
            elemDisp     = globdat.state(elemUDofs);
            elemPhi      = globdat.state(elemPhiDofs);
            elemPhi0     = globdat.state0(elemPhiDofs);
            
            % Increment counter
            elementCounter = elementCounter + 1;
            indexCounter   = indexCounter   + dim^2*nument^2;
            
            % Set local stiffness matrix, residual, displacement to zero
            localkUU      = zeros(dim*nument,dim*nument);
            localkUPhi    = zeros(dim*nument,nument);
            localkPhiPhi  = zeros(nument,nument);
            localfintU    = zeros(dim*nument,1);
            localfintPhi  = zeros(nument,1);
            
            % Loop over Gauss points
            kgauss = 0;
            for ii=1:ngaussX
                for jj=1:ngaussY
                    kgauss = kgauss+1;
                    
                    % Get phase-field and compute degradation function
                    phigp   = globdat.Basis.shape{indexPatch}{i}(kgauss,:) * elemPhi;
                    phigp0  = globdat.Basis.shape{indexPatch}{i}(kgauss,:) * elemPhi0;
                    
                    fac1   = (1-phigp)^material.p;
                    dfac1  = -material.p*(1-phigp)^(material.p-1); 
                    ddfac1 = material.p*(material.p-1)*(1-phigp)^(material.p-2); 
                    
                    fac2   = fac1 + material.a1*phigp + material.a1*material.a2*phigp*phigp + ...
                               material.a1*material.a2*material.a3*phigp*phigp*phigp;
                    dfac2  = dfac1 + material.a1 + 2*material.a1*material.a2*phigp + ...
                               3*material.a1*material.a2*material.a3*phigp*phigp;
                    ddfac2 = ddfac1 + 2*material.a1*material.a2 + ...
                               6*material.a1*material.a2*material.a3*phigp;       
                           
                    gphi     = fac1/fac2;
                    dgphi    = (dfac1*fac2-fac1*dfac2)/(fac2*fac2);
                    ddgphi   = ((ddfac1*fac2-fac1*ddfac2)*fac2-2*...
                               (dfac1*fac2-fac1*dfac2)*dfac2)/(fac2*fac2*fac2);
                           
                    % alpha  = material.eta*phigp + (1-material.eta)*phigp*phigp                           
                    dalpha   = material.eta + 2*(1-material.eta)*phigp;
                    ddalpha  = 2*(1-material.eta);               
                    
                    % Set refinement based on phase-field
                    if (phigp > props.geom.threshPhi) && (globdat.PHTelem{indexPatch}(i).level < props.geom.maxRefLevel)
                        markRef{indexPatch}(i)=1;
                    end
                    
                    % Get shape functions derivatives
                    [Bu,BPhi] = getBMatrices(globdat.Basis.dgdx{indexPatch}{i},...
                                             nument,nstress,dim,kgauss);
                    
                    % Compute gradients (strain and phase-field gradient)
                    strain  = Bu*elemDisp;
                    gradPhi = BPhi*elemPhi;

                    % Strain: Voigt to Tensor
                    trueStrain    = Voigt2Tensor(strain,dim,true);
                    trStrainPos   = max(0,trace(trueStrain));
                    volStrain     = trStrainPos/dim;
                    devStrain     = trueStrain-eye(dim)*volStrain;

                    stressPos     = material.K * MyHeaviside(volStrain,-1e-15) * trace(trueStrain) * eye(dim) + ...
                                    2 * material.mu * devStrain;
                    stressPos     = Tensor2Voigt(stressPos,dim);

                    devStrain     = Tensor2Voigt(devStrain,dim);
                    Psi           = 0.5 * material.K * trStrainPos^2 + ...
                                    material.mu * dot(devStrain,devStrain);
                    volStrainP    = MyHeaviside(volStrain,-1e-15);
                    volStrainN    = 1-volStrainP;

                    % Compute D tensor and stress
                    D = gphi * (material.K * volStrainP * globdat.Pvol + ...
                        2 * material.mu * globdat.Pdev ) + ...
                        material.K * volStrainN * globdat.Pvol;
                    stress = D * strain;

                    % Get initial Psi
                    Psi = max(Psi,material.Psi0); 
                        
                    % Compute local stiffness matrix
                    localkUU     = localkUU + (Bu'*D*Bu).* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    localkUPhi   = localkUPhi + dgphi .* (Bu' * stressPos * globdat.Basis.shape{indexPatch}{i}(kgauss,:)) .* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    localkPhiPhi = localkPhiPhi + 2*material.gc*material.l0/material.calpha* (BPhi'* BPhi).* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    localkPhiPhi = localkPhiPhi + (material.gc/(material.l0*material.calpha)*ddalpha + ddgphi * Psi + material.visc).*(globdat.Basis.shape{indexPatch}{i}(kgauss,:)'* globdat.Basis.shape{indexPatch}{i}(kgauss,:)) .* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    

                    % Compute local residuals
                    localfintU    = localfintU + Bu' * stress .* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    localfintPhi  = localfintPhi + 2*material.gc*material.l0/material.calpha* (BPhi'* gradPhi).* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    localfintPhi  = localfintPhi + (material.gc/(material.l0*material.calpha)*dalpha + dgphi * Psi + material.visc * (phigp-phigp0) ) .* globdat.Basis.shape{indexPatch}{i}(kgauss,:)' .* globdat.Basis.volume{indexPatch}{i}(kgauss);

                end % j Gauss
            end % i Gauss
            
            % Assemble global residual
            fint(elemUDofs)    = fint(elemUDofs) + localfintU;
            fint(elemPhiDofs)  = fint(elemPhiDofs) + localfintPhi;
            
            % Store element matrices in pre-allocated cells
            kUU{elementCounter}     = localkUU;
            kUPhi{elementCounter}   = localkUPhi;
            kPhiPhi{elementCounter} = localkPhiPhi;
            kPhiU{elementCounter}   = localkUPhi';
        end
    end
end
            
clear dsctrx

% COMMENT: Can we include the indexing in the previous loop?
II = zeros(1,indexCounter);
JJ = zeros(1,indexCounter);
S  = zeros(1,indexCounter);

% Set element and index counter to zero
elementCounter = 0;
indexCounter   = 0;

% TO-DO: Can we move this inside the previous loop? We will have to 
% pre-allocate II,JJ,S.
for indexPatch = 1:numPatches
    for i=1:length(globdat.PHTelem{indexPatch})
        if isempty(globdat.PHTelem{indexPatch}(i).children)
            elementCounter = elementCounter+1;
            nument    = size(globdat.PHTelem{indexPatch}(i).C,1);
            sctrx     = globdat.PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            dsctrx    = [reshape([2*sctrx-1;2*sctrx],1,(dim)*nument), udofs + sctrx];
            elmtStiff = [kUU{elementCounter} kUPhi{elementCounter};kPhiU{elementCounter} kPhiPhi{elementCounter}];
            II(indexCounter+1:indexCounter+(dim+1)^2*nument^2) = repmat(dsctrx,1,(dim+1)*nument);
            JJ(indexCounter+1:indexCounter+(dim+1)^2*nument^2) = reshape(repmat(dsctrx,(dim+1)*nument,1),1,(dim+1)^2*nument^2);
            S(indexCounter+1:indexCounter+(dim+1)^2*nument^2)  = reshape(elmtStiff,1,(dim+1)^2*nument^2);
            indexCounter = indexCounter +(dim+1)^2*nument^2;
        end
    end
end

% Assemble (sparse) stiffness matrix
K = sparse(II,JJ,S,ndofs,ndofs);

end



%=========================================================================%
%  GET_MATRIX0_EXT
%=========================================================================%

function [K, fint, markRef] = getMatrix0_Ext(props,globdat)

% Some useful variables
material   = props.mat;
dim        = props.geom.dim;
nstress    = props.geom.nstress;
numPatches = props.geom.numPatches;
ngaussX    = props.geom.ngaussX;
ngaussY    = props.geom.ngaussY;
NN         = globdat.ndofs;
UEnd       = dim*globdat.mesh.sizeBasis;

% Initialise cell for mesh refinement marker
markRef    = cell(numPatches,1);

% Initialise vector for global residual
fint        = zeros(NN,1);

% Initialise cells for local stiffness matrices
nelem      = globdat.mesh.numElements;
kUU        = cell(nelem,1);
kUPhi      = cell(nelem,1);
kPhiU      = cell(nelem,1);
kPhiPhi    = cell(nelem,1);

% Set index counter to zero
elementCounter = 0;
indexCounter   = 0;

% Loop over patches
for indexPatch = 1:numPatches
    
    % Initialise vectors for elements in a patch
    markRef(indexPatch) = {zeros(1,length(globdat.PHTelem{indexPatch}))};
    
    % Loop over elements
    for i=1:length(globdat.PHTelem{indexPatch})
        if isempty(globdat.PHTelem{indexPatch}(i).children)
            
            % Extract element displacement and phase-field
            nument       = size(globdat.PHTelem{indexPatch}(i).C,1);
            elemNodes    = globdat.PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            elemUDofs    = reshape([2*elemNodes-1;2*elemNodes],1,dim*nument);
            elemPhiDofs  = elemNodes+UEnd;
            elemDisp     = globdat.state(elemUDofs);
            elemPhi      = globdat.state(elemPhiDofs);
            elemPhi0     = globdat.state0(elemPhiDofs);
            elemPhiEx    = globdat.Exstate(elemPhiDofs);
            
            % Increment counter
            elementCounter = elementCounter + 1;
            indexCounter   = indexCounter   + dim^2*nument^2;
            
            % Set local stiffness matrix, residual, displacement to zero
            localkUU      = zeros(dim*nument,dim*nument);
            localkUPhi    = zeros(dim*nument,nument);
            localkPhiPhi  = zeros(nument,nument);
            localfintU    = zeros(dim*nument,1);
            localfintPhi  = zeros(nument,1);
            
            % Loop over Gauss points
            kgauss = 0;
            for ii=1:ngaussX
                for jj=1:ngaussY
                    kgauss = kgauss+1;

                    % Get extrapolated phase-field and degradation function
                    phigpEx   = globdat.Basis.shape{indexPatch}{i}(kgauss,:) * elemPhiEx;

                    fac1   = (1-phigpEx)^material.p;
                    fac2   = fac1 + material.a1*phigpEx + material.a1*material.a2*phigpEx*phigpEx + ...
                               material.a1*material.a2*material.a3*phigpEx*phigpEx*phigpEx;      
                           
                    gphiEx     = fac1/fac2;
                    
                    % Get phase-field and compute degradation function
                    phigp   = globdat.Basis.shape{indexPatch}{i}(kgauss,:) * elemPhi;
                    phigp0  = globdat.Basis.shape{indexPatch}{i}(kgauss,:) * elemPhi0;
                    
                    fac1   = (1-phigp)^material.p;
                    dfac1  = -material.p*(1-phigp)^(material.p-1); 
                    ddfac1 = material.p*(material.p-1)*(1-phigp)^(material.p-2); 
                    
                    fac2   = fac1 + material.a1*phigp + material.a1*material.a2*phigp*phigp + ...
                               material.a1*material.a2*material.a3*phigp*phigp*phigp;
                    dfac2  = dfac1 + material.a1 + 2*material.a1*material.a2*phigp + ...
                               3*material.a1*material.a2*material.a3*phigp*phigp;
                    ddfac2 = ddfac1 + 2*material.a1*material.a2 + ...
                               6*material.a1*material.a2*material.a3*phigp;       
                           
                    dgphi    = (dfac1*fac2-fac1*dfac2)/(fac2*fac2);
                    ddgphi   = ((ddfac1*fac2-fac1*ddfac2)*fac2-2*...
                               (dfac1*fac2-fac1*dfac2)*dfac2)/(fac2*fac2*fac2);
                           
                    % alpha  = material.eta*phigp + (1-material.eta)*phigp*phigp                           
                    dalpha   = material.eta + 2*(1-material.eta)*phigp;
                    ddalpha  = 2*(1-material.eta);

                    % Penalisation stuff
                    if phigp - phigp0 > 1e-10
                        delphigp  = 0;
                    else
                        delphigp  = phigp0-phigp;
                    end                
                    
                    % Set refinement based on phase-field
                    if (phigp > props.geom.threshPhi) && (globdat.PHTelem{indexPatch}(i).level < props.geom.maxRefLevel)
                        markRef{indexPatch}(i)=1;
                    end
                    
                    % Get shape functions derivatives
                    [Bu,BPhi] = getBMatrices(globdat.Basis.dgdx{indexPatch}{i},...
                                             nument,nstress,dim,kgauss);
                    
                    % Compute gradients (strain and phase-field gradient)
                    strain  = Bu*elemDisp;
                    gradPhi = BPhi*elemPhi;
                    
                    % Strain: Voigt to Tensor
                    trueStrain    = Voigt2Tensor(strain,dim,true);
                    trStrainPos   = max(0,trace(trueStrain));
                    volStrain     = trStrainPos/dim;
                    devStrain     = trueStrain-eye(dim)*volStrain;

                    stressPos     = material.K * MyHeaviside(volStrain,-1e-15) * trace(trueStrain) * eye(dim) + ...
                                    2 * material.mu * devStrain;
                    stressPos     = Tensor2Voigt(stressPos,dim);

                    devStrain     = Tensor2Voigt(devStrain,dim);
                    Psi           = 0.5 * material.K * trStrainPos^2 + ...
                                    material.mu * dot(devStrain,devStrain);
                    volStrainP    = MyHeaviside(volStrain,-1e-15);
                    volStrainN    = 1-volStrainP;

                    % Compute D tensor and stress
                    D = gphiEx * (material.K * volStrainP * globdat.Pvol + ...
                        2 * material.mu * globdat.Pdev ) + ...
                        material.K * volStrainN * globdat.Pvol;
                    stress = D * strain;
                    
                    % Get initial Psi
                    Psi = max(Psi,material.Psi0); 
                        
                    % Compute local stiffness matrix
                    localkUU     = localkUU + (Bu'*D*Bu).* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    localkUPhi   = localkUPhi + dgphi .* (Bu' * stressPos * globdat.Basis.shape{indexPatch}{i}(kgauss,:)) .* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    localkPhiPhi = localkPhiPhi + 2*material.gc*material.l0/material.calpha* (BPhi'* BPhi).* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    localkPhiPhi = localkPhiPhi + (material.gc/(material.l0*material.calpha)*ddalpha + ddgphi * Psi + material.visc ).*(globdat.Basis.shape{indexPatch}{i}(kgauss,:)'* globdat.Basis.shape{indexPatch}{i}(kgauss,:)) .* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    
                    % Add penalty term
                    
                    if phigp0 > phigp
                        localkPhiPhi = localkPhiPhi + material.pen*material.gc/(material.l0*material.calpha).*(globdat.Basis.shape{indexPatch}{i}(kgauss,:)'* globdat.Basis.shape{indexPatch}{i}(kgauss,:)) .* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    end

                    % Compute local residuals
                    localfintU    = localfintU + Bu' * stress .* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    localfintPhi  = localfintPhi + 2*material.gc*material.l0/material.calpha* (BPhi'* gradPhi).* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    localfintPhi  = localfintPhi + (material.gc/(material.l0*material.calpha)*dalpha + dgphi * Psi - material.pen*material.gc/(material.l0*material.calpha)*delphigp + material.visc * (phigp-phigp0) ) .* globdat.Basis.shape{indexPatch}{i}(kgauss,:)' .* globdat.Basis.volume{indexPatch}{i}(kgauss);

                end % j Gauss
            end % i Gauss
            
            % Assemble global residual
            fint(elemUDofs)    = fint(elemUDofs) + localfintU;
            fint(elemPhiDofs)  = fint(elemPhiDofs) + localfintPhi;
            
            % Store element matrices in pre-allocated cells
            kUU{elementCounter}     = localkUU;
            kUPhi{elementCounter}   = zeros(dim*nument,nument);
            kPhiPhi{elementCounter} = localkPhiPhi;
            kPhiU{elementCounter}   = localkUPhi';
        end
    end
end
            
clear dsctrx

% COMMENT: Can we include the indexing in the previous loop?
II = zeros(1,indexCounter);
JJ = zeros(1,indexCounter);
S  = zeros(1,indexCounter);

% Set element and index counter to zero
elementCounter = 0;
indexCounter   = 0;

% TO-DO: Can we move this inside the previous loop? We will have to 
% pre-allocate II,JJ,S.
for indexPatch = 1:numPatches
    for i=1:length(globdat.PHTelem{indexPatch})
        if isempty(globdat.PHTelem{indexPatch}(i).children)
            elementCounter = elementCounter+1;
            nument    = size(globdat.PHTelem{indexPatch}(i).C,1);
            sctrx     = globdat.PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            dsctrx    = [reshape([2*sctrx-1;2*sctrx],1,(dim)*nument), UEnd + sctrx];
            elmtStiff = [kUU{elementCounter} kUPhi{elementCounter};kPhiU{elementCounter} kPhiPhi{elementCounter}];
            II(indexCounter+1:indexCounter+(dim+1)^2*nument^2) = repmat(dsctrx,1,(dim+1)*nument);
            JJ(indexCounter+1:indexCounter+(dim+1)^2*nument^2) = reshape(repmat(dsctrx,(dim+1)*nument,1),1,(dim+1)^2*nument^2);
            S(indexCounter+1:indexCounter+(dim+1)^2*nument^2)  = reshape(elmtStiff,1,(dim+1)^2*nument^2);
            indexCounter = indexCounter +(dim+1)^2*nument^2;
        end
    end
end

% Assemble (sparse) stiffness matrix
K = sparse(II,JJ,S,NN,NN);

end



%=========================================================================%
%  GET_INT_FORCE
%=========================================================================%

function [fint, markRef] = getIntForce(props,globdat)

% Some useful variables
material   = props.mat;
dim        = props.geom.dim;
nstress    = props.geom.nstress;
numPatches = props.geom.numPatches;
ngaussX    = props.geom.ngaussX;
ngaussY    = props.geom.ngaussY;
ndofs      = globdat.ndofs;
UEnd       = dim*globdat.mesh.sizeBasis;

% Initialise cell for mesh refinement marker
markRef    = cell(numPatches,1);

% Initialise vector for global residual
fint        = zeros(ndofs,1);

% Set index counter to zero
elementCounter = 0;
indexCounter   = 0;

% Loop over patches
for indexPatch = 1:numPatches
    
    % Initialise vectors for elements in a patch
    markRef(indexPatch) = {zeros(1,length(globdat.PHTelem{indexPatch}))};
    
    % Loop over elements
    for i=1:length(globdat.PHTelem{indexPatch})
        if isempty(globdat.PHTelem{indexPatch}(i).children)
            
            % Extract element displacement and phase-field
            nument       = size(globdat.PHTelem{indexPatch}(i).C,1);
            elemNodes    = globdat.PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            elemUDofs    = reshape([2*elemNodes-1;2*elemNodes],1,dim*nument);
            elemPhiDofs  = elemNodes+UEnd;
            elemDisp     = globdat.state(elemUDofs);
            elemPhi      = globdat.state(elemPhiDofs);
            elemPhi0     = globdat.state0(elemPhiDofs);
            
            % Increment counter
            elementCounter = elementCounter + 1;
            indexCounter   = indexCounter   + dim^2*nument^2;
            
            % Set local stiffness matrix, residual, displacement to zero
            localfintU    = zeros(dim*nument,1);
            localfintPhi  = zeros(nument,1);
            
            % Loop over Gauss points
            kgauss = 0;
            for ii=1:ngaussX
                for jj=1:ngaussY
                    kgauss = kgauss+1;
                    
                    % Get phase-field and compute degradation function
                    phigp   = globdat.Basis.shape{indexPatch}{i}(kgauss,:) * elemPhi;
                    phigp0  = globdat.Basis.shape{indexPatch}{i}(kgauss,:) * elemPhi0;
                    
                    fac1   = (1-phigp)^material.p;
                    dfac1  = -material.p*(1-phigp)^(material.p-1); 
                    
                    fac2   = fac1 + material.a1*phigp + material.a1*material.a2*phigp*phigp + ...
                               material.a1*material.a2*material.a3*phigp*phigp*phigp;
                    dfac2  = dfac1 + material.a1 + 2*material.a1*material.a2*phigp + ...
                               3*material.a1*material.a2*material.a3*phigp*phigp;      
                           
                    dgphi    = (dfac1*fac2-fac1*dfac2)/(fac2*fac2);
                           
                    % alpha  = material.eta*phigp + (1-material.eta)*phigp*phigp                           
                    dalpha   = material.eta + 2*(1-material.eta)*phigp;

                    % Penalisation stuff
                    if phigp - phigp0 > 1e-10
                        delphigp  = 0;
                    else
                        delphigp  = phigp0-phigp;
                    end                
                    
                    % Set refinement based on phase-field
                    if (phigp > props.geom.threshPhi) && (globdat.PHTelem{indexPatch}(i).level < props.geom.maxRefLevel)
                        markRef{indexPatch}(i)=1;
                    end
                    
                    % Get shape functions derivatives
                    [Bu,BPhi] = getBMatrices(globdat.Basis.dgdx{indexPatch}{i},...
                                             nument,nstress,dim,kgauss);
                    
                    % Compute gradients (strain and phase-field gradient)
                    strain  = Bu*elemDisp;
                    gradPhi = BPhi*elemPhi;
                                        
                    % Strain: Voigt to Tensor
                    trueStrain    = Voigt2Tensor(strain,dim,true);
                    trStrainPos   = max(0,trace(trueStrain));
                    volStrain     = trStrainPos/dim;
                    devStrain     = trueStrain-eye(dim)*volStrain;

                    devStrain     = Tensor2Voigt(devStrain,dim);
                    Psi           = 0.5 * material.K * trStrainPos^2 + ...
                                    material.mu * dot(devStrain,devStrain);
                    volStrainP    = MyHeaviside(volStrain,-1e-15);
                    volStrainN    = 1-volStrainP;

                    % Compute D tensor and stress
                    D = gphiEx * (material.K * volStrainP * globdat.Pvol + ...
                        2 * material.mu * globdat.Pdev ) + ...
                        material.K * volStrainN * globdat.Pvol;
                    stress = D * strain;
                    
                    % Get initial Psi
                    Psi = max(Psi,material.Psi0); 
                        
                    % Compute local residuals
                    localfintU    = localfintU + Bu' * stress .* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    localfintPhi  = localfintPhi + 2*material.gc*material.l0/material.calpha* (BPhi'* gradPhi).* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    localfintPhi  = localfintPhi + (material.gc/(material.l0*material.calpha)*dalpha + dgphi * Psi - material.pen*material.gc/(material.l0*material.calpha)*delphigp + material.visc * (phigp-phigp0) ) .* globdat.Basis.shape{indexPatch}{i}(kgauss,:)' .* globdat.Basis.volume{indexPatch}{i}(kgauss);

                end % j Gauss
            end % i Gauss
            
            % Assemble global residual
            fint(elemUDofs)    = fint(elemUDofs) + localfintU;
            fint(elemPhiDofs)  = fint(elemPhiDofs) + localfintPhi;
            
        end
    end
end

end


%=========================================================================%
% GET_MATRIX0_1
%=========================================================================%

function [K, fint, markRef] = getMatrix0_1(props,globdat)

% Some useful variables
material   = props.mat;
dim        = props.geom.dim;
nstress    = props.geom.nstress;
numPatches = props.geom.numPatches;
ngaussX    = props.geom.ngaussX;
ngaussY    = props.geom.ngaussY;
NN         = dim*globdat.mesh.sizeBasis;

% Initialise cell for mesh refinement marker
markRef    = cell(numPatches,1);

% Initialise vector for global residual
fint       = zeros(NN,1);

% Initialise cells for local stiffness matrices
nelem      = globdat.mesh.numElements;
kUU        = cell(nelem,1);

% Set index counter to zero
elementCounter = 0;
indexCounter   = 0;

% Loop over patches
for indexPatch = 1:numPatches
    
    % Initialise vectors for elements in a patch
    markRef(indexPatch) = {zeros(1,length(globdat.PHTelem{indexPatch}))};
    
    % Loop over elements
    for i=1:length(globdat.PHTelem{indexPatch})
        if isempty(globdat.PHTelem{indexPatch}(i).children)
            
            % Extract element displacement and phase-field
            nument       = size(globdat.PHTelem{indexPatch}(i).C,1);
            elemNodes    = globdat.PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            elemUDofs    = reshape([2*elemNodes-1;2*elemNodes],1,dim*nument);
            elemPhiDofs  = elemNodes+NN;
            elemDisp     = globdat.state(elemUDofs);
            elemPhi      = globdat.state(elemPhiDofs);
            
            % Increment counter
            elementCounter = elementCounter + 1;
            indexCounter   = indexCounter + dim^2*nument^2;
            
            % Set local stiffness matrix, residual, displacement to zero
            localkUU      = zeros(dim*nument,dim*nument);
            localfintU    = zeros(dim*nument,1);
            
            % Loop over Gauss points
            kgauss = 0;
            for ii=1:ngaussX
                for jj=1:ngaussY
                    kgauss = kgauss+1;
                    
                    % Get phase-field and compute degradation function
                    phigp   = globdat.Basis.shape{indexPatch}{i}(kgauss,:) * elemPhi;
                    
                    fac1   = (1-phigp)^material.p;
                    fac2   = fac1 + material.a1*phigp + material.a1*material.a2*phigp*phigp + ...
                               material.a1*material.a2*material.a3*phigp*phigp*phigp;
                    gphi   = fac1/fac2;                    
                    
                    % Set refinement based on phase-field
                    if (phigp > props.geom.threshPhi) && (globdat.PHTelem{indexPatch}(i).level < props.geom.maxRefLevel)
                        markRef{indexPatch}(i)=1;
                    end
                    
                    % Get shape functions derivatives
                    [Bu,~] = getBMatrices(globdat.Basis.dgdx{indexPatch}{i},...
                                             nument,nstress,dim,kgauss);
                    
                    % Compute strain 
                    strain  = Bu*elemDisp;

                    volStrain     = 0.5*strain(1)+0.5*strain(2);
                    volStrainP    = MyHeaviside(volStrain,-1e-15);
                    volStrainN    = 1-volStrainP;

                    % Compute D tensor and stress
                    D = gphi * (material.K * volStrainP * globdat.Pvol + ...
                        2 * material.mu * globdat.Pdev ) + ...
                        material.K * volStrainN * globdat.Pvol;
                    
                    stress = D * strain;
                                            
                    % Compute local stiffness matrix
                    localkUU     = localkUU + (Bu'*D*Bu).* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    
                    % Compute local fintU
                    localfintU    = localfintU + Bu' * stress .* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    
                end % j Gauss
            end % i Gauss
            
            % Assemble global residual
            fint(elemUDofs)    = fint(elemUDofs) + localfintU;
            
            % Store element matrices in pre-allocated cells
            kUU{elementCounter}     = localkUU;
        end
    end
end
            
clear dsctrx

% COMMENT: Can we include the indexing in the previous loop?
II = zeros(1,indexCounter);
JJ = zeros(1,indexCounter);
S  = zeros(1,indexCounter);

% Set element and index counter to zero
elementCounter = 0;
indexCounter   = 0;

% TO-DO: Can we move this inside the previous loop? We will have to 
% pre-allocate II,JJ,S.
for indexPatch = 1:numPatches
    for i=1:length(globdat.PHTelem{indexPatch})
        if isempty(globdat.PHTelem{indexPatch}(i).children)
            elementCounter = elementCounter+1;
            nument    = size(globdat.PHTelem{indexPatch}(i).C,1);
            sctrx     = globdat.PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            dsctrx = reshape([2*sctrx-1;2*sctrx],1,dim*nument); 
            elmtStiff = [kUU{elementCounter}];
            II(indexCounter+1:indexCounter+dim^2*nument^2) = repmat(dsctrx,1,dim*nument);
            JJ(indexCounter+1:indexCounter+dim^2*nument^2) = reshape(repmat(dsctrx,dim*nument,1),1,dim^2*nument^2);
            S(indexCounter+1:indexCounter+dim^2*nument^2)  = reshape(elmtStiff,1,dim^2*nument^2);
            indexCounter = indexCounter +dim^2*nument^2;
        end
    end
end

% Assemble (sparse) stiffness matrix
K = sparse(II,JJ,S,NN,NN);

end


%=========================================================================%
% GET_MATRIX0_2
%=========================================================================%

function [K, fint] = getMatrix0_2(props,globdat)


% Some useful variables
material   = props.mat;
dim        = props.geom.dim;
nstress    = props.geom.nstress;
numPatches = props.geom.numPatches;
ngaussX    = props.geom.ngaussX;
ngaussY    = props.geom.ngaussY;
NN         = globdat.mesh.sizeBasis;
UEnd       = dim*globdat.mesh.sizeBasis;

% Initialise vector for global residual
fint        = zeros(NN,1);

% Initialise cells for local stiffness matrices
nelem      = globdat.mesh.numElements;
kPhiPhi    = cell(nelem,1);

% Set index counter to zero
elementCounter = 0;
indexCounter   = 0;

% Loop over patches
for indexPatch = 1:numPatches
    
    % Loop over elements
    for i=1:length(globdat.PHTelem{indexPatch})
        if isempty(globdat.PHTelem{indexPatch}(i).children)
            
            % Extract element displacement and phase-field
            nument       = size(globdat.PHTelem{indexPatch}(i).C,1);
            elemNodes    = globdat.PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            elemUDofs    = reshape([2*elemNodes-1;2*elemNodes],1,dim*nument);
            elemPhiDofs  = elemNodes+UEnd;
            elemDisp     = globdat.state(elemUDofs);
            elemPhi      = globdat.state(elemPhiDofs);
            elemPhi0     = globdat.state0(elemPhiDofs);
            
            % Increment counter
            elementCounter = elementCounter + 1;
            indexCounter   = indexCounter + nument^2;
            
            % Set local stiffness matrix, residual, displacement to zero
            localkPhiPhi  = zeros(nument,nument);
            localfintPhi  = zeros(nument,1);
            
            % Loop over Gauss points
            kgauss = 0;
            for ii=1:ngaussX
                for jj=1:ngaussY
                    kgauss = kgauss+1;
                    
                    % Get phase-field and compute degradation function
                    phigp   = globdat.Basis.shape{indexPatch}{i}(kgauss,:) * elemPhi;
                    phigp0  = globdat.Basis.shape{indexPatch}{i}(kgauss,:) * elemPhi0;
                    
                    fac1   = (1-phigp)^material.p;
                    dfac1  = -material.p*(1-phigp)^(material.p-1); 
                    ddfac1 = material.p*(material.p-1)*(1-phigp)^(material.p-2); 
                    
                    fac2   = fac1 + material.a1*phigp + material.a1*material.a2*phigp*phigp + ...
                               material.a1*material.a2*material.a3*phigp*phigp*phigp;
                    dfac2  = dfac1 + material.a1 + 2*material.a1*material.a2*phigp + ...
                               3*material.a1*material.a2*material.a3*phigp*phigp;
                    ddfac2 = ddfac1 + 2*material.a1*material.a2 + ...
                               6*material.a1*material.a2*material.a3*phigp;       
                           
                    % gphi     = fac1/fac2;
                    dgphi    = (dfac1*fac2-fac1*dfac2)/(fac2*fac2);
                    ddgphi   = ((ddfac1*fac2-fac1*ddfac2)*fac2-2*...
                               (dfac1*fac2-fac1*dfac2)*dfac2)/(fac2*fac2*fac2);
                           
                    % alpha  = material.eta*phigp + (1-material.eta)*phigp*phigp                           
                    dalpha   = material.eta + 2*(1-material.eta)*phigp;
                    ddalpha  = 2*(1-material.eta);
                    
                    % Penalisation stuff
                    if phigp - phigp0 > 1e-10
                        delphigp  = 0;
                    else
                        delphigp  = phigp0-phigp;
                    end                    
                    
                    % Get shape functions derivatives
                    [Bu,BPhi] = getBMatrices(globdat.Basis.dgdx{indexPatch}{i},...
                                             nument,nstress,dim,kgauss);
                    
                    % Compute gradients (strain and phase-field gradient)
                    strain  = Bu*elemDisp;
                    gradPhi = BPhi*elemPhi;

                    trueStrain    = Voigt2Tensor(strain,dim,true);
                    trStrainPos   = max(0,trace(trueStrain));
                    volStrain     = trStrainPos/dim;
                    devStrain     = trueStrain-eye(dim)*volStrain;
                    devStrain     = Tensor2Voigt(devStrain,dim);
                    Psi           = 0.5 * material.K * trStrainPos^2 + ...
                                    material.mu * dot(devStrain,devStrain);

                    % Get initial Psi (required for AT1 and PFCZM)

                    Psi = max(Psi,material.Psi0);
                        
                    % Compute local stiffness matrix
                    
                    localkPhiPhi = localkPhiPhi + 2*material.gc*material.l0/material.calpha* (BPhi'* BPhi).* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    localkPhiPhi = localkPhiPhi + (material.gc/(material.l0*material.calpha)*ddalpha + ddgphi * Psi + material.visc).*(globdat.Basis.shape{indexPatch}{i}(kgauss,:)'* globdat.Basis.shape{indexPatch}{i}(kgauss,:)) .* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    
                    % Add penalty term
                    
                    if phigp0 > phigp
                        localkPhiPhi = localkPhiPhi + material.pen*material.gc/material.l0.*(globdat.Basis.shape{indexPatch}{i}(kgauss,:)'* globdat.Basis.shape{indexPatch}{i}(kgauss,:)) .* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    end
                    
                    % Compute local residuals
                    
                    localfintPhi  = localfintPhi + 2*material.gc*material.l0/material.calpha* (BPhi'* gradPhi).* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    localfintPhi  = localfintPhi + (material.gc/(material.l0*material.calpha)*dalpha + dgphi * Psi - material.pen*material.gc/material.l0*delphigp + material.visc * (phigp-phigp0) ) .* globdat.Basis.shape{indexPatch}{i}(kgauss,:)' .* globdat.Basis.volume{indexPatch}{i}(kgauss);

                end % j Gauss point
            end     % i Gauss point
            
            % Assemble global residual
            
            fint(elemPhiDofs-UEnd)  =fint(elemPhiDofs-UEnd) + localfintPhi;
            
            % Store element matrices in pre-allocated cells
            
            kPhiPhi{elementCounter} = localkPhiPhi;

        end
    end
end
            
clear dsctrx

% COMMENT: Can we include the indexing in the previous loop?

II = zeros(1,indexCounter);
JJ = zeros(1,indexCounter);
S  = zeros(1,indexCounter);

% Set element and index counter to zero
elementCounter = 0;
indexCounter   = 0;

% TO-DO: Can we move this inside the previous loop? We will have to 
% pre-allocate II,JJ,S.
for indexPatch = 1:numPatches
    for i=1:length(globdat.PHTelem{indexPatch})
        if isempty(globdat.PHTelem{indexPatch}(i).children)
            elementCounter = elementCounter+1;
            nument    = size(globdat.PHTelem{indexPatch}(i).C,1);
            sctrx     = globdat.PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            elmtStiff = [kPhiPhi{elementCounter}];
            II(indexCounter+1:indexCounter+nument^2) = repmat(sctrx,1,nument);
            JJ(indexCounter+1:indexCounter+nument^2) = reshape(repmat(sctrx,nument,1),1,nument^2);
            S(indexCounter+1:indexCounter+nument^2)  = reshape(elmtStiff,1,nument^2);
            indexCounter = indexCounter+nument^2;
        end
    end
end

% Assemble (sparse) stiffness matrix
K = sparse(II,JJ,S,NN,NN);

end


%=========================================================================%
% GET_DISSIPATION
%=========================================================================%

function [h, g] = getDissipation(props,globdat)

% Some useful variables
material   = props.mat;
dim        = props.geom.dim;
numPatches = props.geom.numPatches;
nstress    = props.geom.nstress;
ngaussX    = props.geom.ngaussX;
ngaussY    = props.geom.ngaussY;
NN         = globdat.ndofs;
UEnd       = dim*globdat.mesh.sizeBasis;

% Initialise h
h = zeros(NN,1);
g = 0;

% Set index counter to zero
elementCounter = 0;
indexCounter   = 0;

% Loop over patches
for indexPatch = 1:numPatches
    
    % Loop over elements
    for i=1:length(globdat.PHTelem{indexPatch})
        if isempty(globdat.PHTelem{indexPatch}(i).children)
            
            % Extract element displacement and phase-field
            nument         = size(globdat.PHTelem{indexPatch}(i).C,1);
            elemNodes      = globdat.PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            elemPhiDofs    = elemNodes+UEnd;
            elemPhi        = globdat.state(elemPhiDofs);
            DelemPhi       = globdat.Dstate(elemPhiDofs);
            
            % Increment counter
            elementCounter = elementCounter + 1;
            indexCounter   = indexCounter + dim^2*nument^2;
            
            localh  = zeros(nument,1);
            
            % Loop over Gauss points
            kgauss = 0;
            for ii=1:ngaussX
                for jj=1:ngaussY
                    kgauss = kgauss+1;
                    
                    % Get phase-field and compute degradation function
                    phigp  = globdat.Basis.shape{indexPatch}{i}(kgauss,:) * elemPhi;
                    Dphigp = globdat.Basis.shape{indexPatch}{i}(kgauss,:) * DelemPhi;

                    % Compute local dissipation function
                    dalpha   = material.eta + 2*(1-material.eta)*phigp;
                    ddalpha  = 2*(1-material.eta);
                    
                    % Get shape functions derivatives
                    [~,BPhi] = getBMatrices(globdat.Basis.dgdx{indexPatch}{i},...
                                             nument,nstress,dim,kgauss);
                    
                    % Compute phase-field gradient
                    gradPhi    = BPhi*elemPhi;
                    DgradPhi   = BPhi*DelemPhi;
                                       
                    % Compute g
                    g = g + material.gc/(material.calpha*material.l0).* dalpha*Dphigp.* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    g = g + 2.0*(material.gc*material.l0/material.calpha).* dot(gradPhi,DgradPhi).* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    
                    % Compute h
                    localh = localh + material.gc/(material.calpha*material.l0)*(ddalpha*Dphigp+dalpha).*globdat.Basis.shape{indexPatch}{i}(kgauss,:)'.* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    localh = localh + 2.0*(material.gc*material.l0/material.calpha)*( BPhi'*(gradPhi+DgradPhi)).* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    
                end % j Gauss
            end % i Gauss
            
            % Assemble global residual
            h(elemPhiDofs)  = h(elemPhiDofs) + localh;

        end
    end
end

end




%=========================================================================%
% GET_AT2_DISSIPATION
%=========================================================================%

function [h, g] = getAT2Dissipation(props,globdat)

% Some useful variables
material   = props.mat;
dim        = props.geom.dim;
numPatches = props.geom.numPatches;
nstress    = props.geom.nstress;
ngaussX    = props.geom.ngaussX;
ngaussY    = props.geom.ngaussY;
NN         = globdat.ndofs;
UEnd       = dim*globdat.mesh.sizeBasis;

% Initialise h
h = zeros(NN,1);
g = 0;

% Set index counter to zero
elementCounter = 0;
indexCounter   = 0;

% Loop over patches
for indexPatch = 1:numPatches
    
    % Loop over elements
    for i=1:length(globdat.PHTelem{indexPatch})
        if isempty(globdat.PHTelem{indexPatch}(i).children)
            
            % Extract element displacement and phase-field
            nument         = size(globdat.PHTelem{indexPatch}(i).C,1);
            elemNodes      = globdat.PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            elemPhiDofs    = elemNodes+UEnd;
            elemPhi        = globdat.state(elemPhiDofs);
            DelemPhi       = globdat.Dstate(elemPhiDofs);
            
            % Increment counter
            elementCounter = elementCounter + 1;
            indexCounter   = indexCounter + dim^2*nument^2;
            
            localh  = zeros(nument,1);
            
            % Loop over Gauss points
            kgauss = 0;
            for ii=1:ngaussX
                for jj=1:ngaussY
                    kgauss = kgauss+1;
                    
                    % Get phase-field and compute degradation function
                    phigp  = globdat.Basis.shape{indexPatch}{i}(kgauss,:) * elemPhi;
                    Dphigp = globdat.Basis.shape{indexPatch}{i}(kgauss,:) * DelemPhi;
                    
                    % Get shape functions derivatives
                    [~,BPhi] = getBMatrices(globdat.Basis.dgdx{indexPatch}{i},...
                                             nument,nstress,dim,kgauss);
                    
                    % Compute phase-field gradient
                    gradPhi    = BPhi*elemPhi;
                    DgradPhi   = BPhi*DelemPhi;
                                       
                    % Compute g
                    g = g + material.gc/material.l0.* phigp*Dphigp.* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    g = g + material.gc*material.l0.* dot(gradPhi,DgradPhi).* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    
                    % Compute h
                    localh = localh + material.gc/material.l0*(phigp+Dphigp).*globdat.Basis.shape{indexPatch}{i}(kgauss,:)'.* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    localh = localh + material.gc*material.l0*( BPhi'*(gradPhi+DgradPhi)).* globdat.Basis.volume{indexPatch}{i}(kgauss);
                    
                end % j Gauss
            end % i Gauss
            
            % Assemble global residual
            h(elemPhiDofs)  = h(elemPhiDofs) + localh;

        end
    end
end

end




%=========================================================================%
% GET_MATRIX20
%=========================================================================%

function M = getMatrix20(props,globdat)

M       = [];

end


%=========================================================================%
% GET_MATRIX21
%=========================================================================%

function M = getMatrix21(props,globdat)

M       = [];

end


%-------------------------------------------------------------------------%
% B-matrices                                                              %
%-------------------------------------------------------------------------%

function [Bv, Bs] = getBMatrices(dgdx,nument,nstress,dim,kgauss)

Bv = zeros(nstress,dim*nument);
Bs = zeros(dim,nument);

for inode=1:nument

    % 2D B Matrix
    
    Bv(1,2*inode-1)  = dgdx(kgauss,1,inode);
    Bv(2,2*inode)    = dgdx(kgauss,2,inode);
    Bv(3,2*inode-1)  = dgdx(kgauss,2,inode);
    Bv(3,2*inode)    = dgdx(kgauss,1,inode);

    Bs(1,inode)=dgdx(kgauss,1,inode);
    Bs(2,inode)=dgdx(kgauss,2,inode); 

end

end


