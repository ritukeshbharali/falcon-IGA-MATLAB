%-----------------------------------------------------------------------------
%  otimesMatsU (C_ijkl = A_ik B_jl)
%-----------------------------------------------------------------------------

function C = otimesMatsU(A,B)

[arows,acols] = size(A);
[brows,bcols] = size(B);

C = zeros(arows,acols,brows,bcols);

for i = 1:arows
    for j = 1:acols
        for k = 1:brows
            for l = 1:bcols
                C(i,j,k,l) = A(i,k) * B(j,l);
            end
        end
    end
end

end