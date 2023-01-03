%-----------------------------------------------------------------------------
%  doubleDotQuadixMat (C_ij = A_ijkl B_kl)
%-----------------------------------------------------------------------------

function C = doubleDotQuadixMat(A,B)

[in,jn,kn,ln] = size(A);
C             = zeros(in,jn);

for i = 1:in
    for j = 1:jn
        for k = 1:kn
            for l = ln
                C(i,j) = A(i,j,k,l) * B(k,l);
            end
        end
    end
end

end