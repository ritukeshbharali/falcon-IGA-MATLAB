function [dirichlet] = DirichletConstraints(PHTelem,geometry,NN,controlPts)
% Initialize the boundary conditions

p = geometry.p;
q = geometry.q;

down_nodes = 1:p+1;
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);

bottomEdge = [];
topEdge = [];
leftEdge = [];
boundaryEdge = [];

% Set the top boundary degrees of freedom
for patchIndex = 1
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            if isempty(PHTelem{patchIndex}(i).neighbor_up)
                topEdge = [topEdge, PHTelem{patchIndex}(i).nodesGlobal(up_nodes)];
                %                 topEdge_local = [topEdge_local, PHTelem{patchIndex}(i).nodes(up_nodes)];
            end
        end
    end
end

% Set the bottom boundary degrees of freedom
for patchIndex = 3
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            if isempty(PHTelem{patchIndex}(i).neighbor_up)
                bottomEdge = [bottomEdge, PHTelem{patchIndex}(i).nodesGlobal(up_nodes)];
                %                 bottomEdge_local = [bottomEdge_local, PHTelem{patchIndex}(i).nodes(up_nodes)];
                
            end
        end
    end
end

% Set the bottom boundary degrees of freedom
for patchIndex = 2
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            if isempty(PHTelem{patchIndex}(i).neighbor_up) && isempty(PHTelem{patchIndex}(i).neighbor_left)
                leftEdge = [leftEdge, PHTelem{patchIndex}(i).nodesGlobal(up_nodes)];
%                 leftEdge_local = [leftEdge_local, PHTelem{patchIndex}(i).nodes(up_nodes)];
                
                boundaryEdge = [boundaryEdge, PHTelem{patchIndex}(i).nodesGlobal(left_nodes)];
%                 boundaryEdge_local = [boundaryEdge_local, PHTelem{patchIndex}(i).nodes(left_nodes)];
            end
        end
    end
end


% Eliminate duplicate entries
Bottom = unique(bottomEdge,'stable');
Top = unique(topEdge,'stable');
cornerElem = intersect(leftEdge,boundaryEdge);

dofs = [2*Bottom-1, 2*Bottom ...
    2*Top, NN+Bottom, NN+cornerElem];

vals     = zeros(1,length(dofs));
l_bottom = length(Bottom);
l_top    = length(Top);
vals(2*l_bottom+1:2*l_bottom+l_top) = 1.0;

% Fill in the dirichlet dofs and their values
dirichlet = [dofs; vals]';

end