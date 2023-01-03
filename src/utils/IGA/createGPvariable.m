function GPvar = createGPvariable(PHTelem,geometry)

ngaussX    = geometry.ngaussX;
ngaussY    = geometry.ngaussY;
numPatches = length(PHTelem);
GPvar      = cell(1,numPatches);

for indexPatch = 1:length(PHTelem)
    GPvar(indexPatch) = {cell(1,length(PHTelem{indexPatch}))};
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            GPvar{indexPatch}{i} = zeros(1,ngaussX*ngaussY);
        end
    end
end
end












