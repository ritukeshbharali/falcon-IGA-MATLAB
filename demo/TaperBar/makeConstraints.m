function [dirichlet] = DirichletConstraints(PHTelem,geometry,NN)

% Some useful variables
p = geometry.p;
q = geometry.q;

% Define side node indices
down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);

% TO-DO: Pre-allocate and then remove empty entries for speed up as the
% array re-sizing happens only once and not for every element on the
% boundary.

bottomEdge = [];
topEdge = [];
leftEdge = [];
rightEdge = [];
cornerElemDofs = [];

% Set the left boundary degrees of freedom
for patchIndex = 1:2
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            if isempty(PHTelem{patchIndex}(i).neighbor_left)
                leftEdge = [leftEdge, PHTelem{patchIndex}(i).nodesGlobal(left_nodes)];
            end
        end
    end
end

% Set the right boundary degrees of freedom
for patchIndex = 1:2
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            if isempty(PHTelem{patchIndex}(i).neighbor_right)
                rightEdge = [rightEdge, PHTelem{patchIndex}(i).nodesGlobal(right_nodes)];
            end
        end
    end
end

dofs      = [2*rightEdge-1, 2*leftEdge-1, 2*leftEdge];
vals      = zeros(1,length(dofs));
vals(1:length(rightEdge)) = 1.0;

% Fill in the dirichlet dofs and their values
dirichlet = [dofs; vals]';

end