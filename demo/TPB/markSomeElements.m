function [cornerElems] = markSomeElements(PHTelem)

% Get corner elements

cornerElems = [];

% Set the fixed boundary degrees of freedom
for patchIndex = 1
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            if isempty(PHTelem{patchIndex}(i).neighbor_down) && isempty(PHTelem{patchIndex}(i).neighbor_left)
                cornerElems = [cornerElems; patchIndex, i];
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
                cornerElems = [cornerElems; patchIndex, i];
                break;
            end
        end
    end
end

end