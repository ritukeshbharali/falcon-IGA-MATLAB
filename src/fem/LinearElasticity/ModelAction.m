%-------------------------------------------------------------------------%
% ModelAction is a function that carries out an action requested by the 
% solver. Here is a list of actions for the Elasticity FE Model:
%    - INIT
%    - GET_MATRIX0
%    - GET_INT_FORCE
%    - GET_MATRIX2
%    - COMMIT
%-------------------------------------------------------------------------%

function globdat = ModelAction(props,globdat,action)

switch action

    case 'INIT'

        globdat = init(props,globdat);
        
    case 'GET_MATRIX0'
        
        [globdat.K, globdat.fint] = getMatrix0(props,globdat);

    case 'GET_INT_FORCE'
        
        globdat.fint = getIntForce(props,globdat);    
        
    case 'GET_MATRIX20'
        
        globdat.M = getMatrix20(props,globdat);
        
    case 'GET_MATRIX21'
        
        globdat.M = getMatrix21(props,globdat);     
        
    case 'COMMIT'    
        
        % Nothing to do

end

end


%-------------------------------------------------------------------------%
% INITIALIZATION                                                          %
%-------------------------------------------------------------------------%

function globdat = init(props,globdat)

globdat.dim   = props.geom.dim;
globdat.ndofs = globdat.dim * globdat.mesh.sizeBasis;

end


%-------------------------------------------------------------------------%
% GET_MATRIX0                                                             %
%-------------------------------------------------------------------------%

function [K, fint] = getMatrix0(props,globdat)

% Get data from globdat

Sol        = globdat.state;
PHTelem    = globdat.PHTelem;
MeshInfo   = globdat.MeshInfo;
Basis      = globdat.Basis;

% Get geometry and material data

geom       = props.geom;
mat        = props.mat;

% Some useful variables

dim        = geom.dim;
numPatches = geom.numPatches;
ngaussX    = geom.ngaussX;
ngaussY    = geom.ngaussY;
nstress    = geom.nstress;
ndofs      = globdat.ndofs;
nelems     = MeshInfo.numElements;

% Initialize data 

fint       = zeros(ndofs,1);
K          = cell(nelems,1);

% Set index counter to zero

elementCounter = 0;
indexCounter   = 0;

% Loop over patches

for idxPatch = 1:numPatches

    % Loop over elements

    for i=1:length(PHTelem{idxPatch})

        % Compute element matrices and internal forces for
        % elements which do not have children elements.

        if isempty(PHTelem{idxPatch}(i).children)
            
            % Extract element displacements

            nument  = size(PHTelem{idxPatch}(i).C,1);
            elNodes = PHTelem{idxPatch}(i).nodesGlobal(1:nument);

            switch dim

                case 1
                    elDofs = elNodes;

                case 2
                    elDofs  = reshape([2*elNodes-1;2*elNodes],...
                                       1,dim*nument);

                case 3
                    elDofs  = reshape([3*elNodes-2;3*elNodes-1;...
                                       3*elNodes],1,dim*nument);
            end

            elDisp  = Sol(elDofs);
            
            % Increment counter

            elementCounter = elementCounter + 1;
            indexCounter   = indexCounter   + dim^2*nument^2;
            
            % Set local stiffness matrix, residual, displacement to zero

            localK        = zeros(dim*nument,dim*nument);
            localfint     = zeros(dim*nument,1);

            % Loop over Gauss points

            ipoint = 0;
            
            for ii=1:ngaussX
                for jj=1:ngaussY

                    ipoint = ipoint+1;
                    
                    % Get shape functions derivatives

                    B = getBMatrix(globdat.Basis.dgdx{indexPatch}{i},...
                                   nument,nstress,dim,ipoint);

                    % Compute strain and stress

                    strain = B     * elDisp;
                    stress = mat.C * strain;

                    % Compute stiffness matrix and internal force

                    wip        = Basis.volume{idxPatch}{i}(ipoint);
                    localK     = localK + wip * ( B' * mat.C * B  );
                    localfint  = localfint + wip * ( B' * stress );
                    
                end % j Gauss
            end % i Gauss
            
            % Assemble global residual
            
            fint(elDofs)    = fint(elDofs) + localfint;
            
            % Store element matrices in pre-allocated cells
            
            K{elementCounter} = localK;

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

for idxPatch = 1:numPatches
    for i=1:length(PHTelem{idxPatch})
        if isempty(PHTelem{idxPatch}(i).children)
            
            elementCounter = elementCounter+1;
            nument         = size(PHTelem{idxPatch}(i).C,1);
            elNodes        = PHTelem{idxPatch}(i).nodesGlobal(1:nument);

            switch dim

                case 1
                    elDofs = elNodes;

                case 2
                    elDofs  = reshape([2*elNodes-1;2*elNodes],...
                                       1,dim*nument);

                case 3
                    elDofs  = reshape([3*elNodes-2;3*elNodes-1;...
                                       3*elNodes],1,dim*nument);
            end

            elemK = [K{elementCounter}];
            II(indexCounter+1:indexCounter+dim^2*nument^2) = repmat(elDofs,1,dim*nument);
            JJ(indexCounter+1:indexCounter+dim^2*nument^2) = reshape(repmat(elDofs,dim*nument,1),1,dim^2*nument^2);
            S(indexCounter+1:indexCounter+dim^2*nument^2)  = reshape(elemK,1,dim^2*nument^2);
            indexCounter = indexCounter +dim^2*nument^2;
        end
    end
end

% Assemble (sparse) stiffness matrix
K = sparse(II,JJ,S,ndofs,ndofs);

end


%-------------------------------------------------------------------------%
% GET_INT_FORCE                                                           %
%-------------------------------------------------------------------------%

function fint = getIntForce(props,globdat)

% Get data from globdat

Sol        = globdat.state;
PHTelem    = globdat.PHTelem;
Basis      = globdat.Basis;

% Get geometry and material data

geom       = props.geom;
mat        = props.mat;

% Some useful variables

dim        = geom.dim;
numPatches = geom.numPatches;
ngaussX    = geom.ngaussX;
ngaussY    = geom.ngaussY;
nstress    = geom.nstress;
ndofs      = globdat.ndofs;

% Initialize data 

fint       = zeros(ndofs,1);

% Set index counter to zero

elementCounter = 0;
indexCounter   = 0;

% Loop over patches

for idxPatch = 1:numPatches

    % Loop over elements

    for i=1:length(PHTelem{idxPatch})

        % Compute element matrices and internal forces for
        % elements which do not have children elements.

        if isempty(PHTelem{idxPatch}(i).children)
            
            % Extract element displacements

            nument  = size(PHTelem{idxPatch}(i).C,1);
            elNodes = PHTelem{idxPatch}(i).nodesGlobal(1:nument);

            switch dim

                case 1
                    elDofs = elNodes;

                case 2
                    elDofs  = reshape([2*elNodes-1;2*elNodes],...
                                       1,dim*nument);

                case 3
                    elDofs  = reshape([3*elNodes-2;3*elNodes-1;...
                                       3*elNodes],1,dim*nument);
            end

            elDisp  = Sol(elDofs);
            
            % Increment counter

            elementCounter = elementCounter + 1;
            indexCounter   = indexCounter   + dim^2*nument^2;
            
            % Set local internal force to zero

            localfint     = zeros(dim*nument,1);

            % Loop over Gauss points

            ipoint = 0;
            
            for ii=1:ngaussX
                for jj=1:ngaussY

                    ipoint = ipoint+1;
                    
                    % Get shape functions derivatives

                    B = getBMatrix(globdat.Basis.dgdx{indexPatch}{i},...
                                   nument,nstress,dim,ipoint);

                    % Compute strain and stress

                    strain = B     * elDisp;
                    stress = mat.C * strain;

                    % Compute stiffness matrix and internal force

                    wip        = Basis.volume{idxPatch}{i}(ipoint);
                    localfint  = localfint + wip * ( B' * stress );
                    
                end % j Gauss
            end % i Gauss
            
            % Assemble global residual
            
            fint(elDofs)    = fint(elDofs) + localfint;

        end
    end
end

end


%-------------------------------------------------------------------------%
% Assemble Lumped Mass                                                    %
%-------------------------------------------------------------------------%

function M = getMatrix20(props,globdat)

M       = [];

end




%-------------------------------------------------------------------------%
% Assemble Consistent Mass                                                %
%-------------------------------------------------------------------------%

function M = getMatrix21(props,globdat)

M       = [];

end


%-------------------------------------------------------------------------%
% B-matrix                                                                %
%-------------------------------------------------------------------------%

function B = getBMatrix(dgdx,nument,nstress,dim,kgauss)

B = zeros(nstress,dim*nument);

for inode=1:nument

    % 2D B Matrix
    
    B(1,2*inode-1)  = dgdx(kgauss,1,inode);
    B(2,2*inode)    = dgdx(kgauss,2,inode);
    B(3,2*inode-1)  = dgdx(kgauss,2,inode);
    B(3,2*inode)    = dgdx(kgauss,1,inode);
end

end
