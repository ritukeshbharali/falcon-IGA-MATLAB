function [dirichlet] = makeConstraints(PHTelem,geometry,NN)
% Initialize the boundary conditions

p = geometry.p;
q = geometry.q;

down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
cornerElem = [];
cornerElemDofs = [];
left_nodes = 1:(p+1):(1+(p+1)*q);

% Set the fixed boundary degrees of freedom
for patchIndex = 1
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            if isempty(PHTelem{patchIndex}(i).neighbor_down) && isempty(PHTelem{patchIndex}(i).neighbor_left)
                bottomEdge = PHTelem{patchIndex}(i).nodesGlobal(down_nodes);
                fixedSupport = bottomEdge(1);
                cornerElem = [cornerElem; patchIndex, i];
                cornerElemDofs = [cornerElemDofs, bottomEdge];
                break;
            end
        end
    end
end

% Set the fixed boundary degrees of freedom
for patchIndex = 4
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            if isempty(PHTelem{patchIndex}(i).neighbor_down) && isempty(PHTelem{patchIndex}(i).neighbor_right)
                bottomEdge = PHTelem{patchIndex}(i).nodesGlobal(down_nodes);
                rollerSupport = bottomEdge(end);
                cornerElem = [cornerElem; patchIndex, i];
                cornerElemDofs = [cornerElemDofs, bottomEdge];
                break;
            end
        end
    end
end

% Set the loading boundary degrees of freedom
for patchIndex = 2
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            if isempty(PHTelem{patchIndex}(i).neighbor_up) && isempty(PHTelem{patchIndex}(i).neighbor_right)
                topEdge = PHTelem{patchIndex}(i).nodesGlobal(up_nodes);
                loadingPt = topEdge(end);
                cornerElemDofs = [cornerElemDofs, topEdge];
                break;
            end
        end
    end
end


% Set the loading boundary degrees of freedom
for patchIndex = 3
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            if isempty(PHTelem{patchIndex}(i).neighbor_up) && isempty(PHTelem{patchIndex}(i).neighbor_left)
                topEdge = PHTelem{patchIndex}(i).nodesGlobal(up_nodes);
                loadingPt = [loadingPt, topEdge];
                cornerElemDofs = [cornerElemDofs, topEdge];
                break;
            end
        end
    end
end

dofs = [2*loadingPt, 2*fixedSupport-1, 2*fixedSupport ...
        2*rollerSupport, ...
        NN + cornerElemDofs];
		
vals                      = zeros(1,length(dofs));	
vals(1:length(loadingPt)) = 1.0;

% Fill in the dirichlet dofs and their values
dirichlet = [dofs; vals]';

end