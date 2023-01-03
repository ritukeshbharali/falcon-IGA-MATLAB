%-----------------------------------------------------------------------------
%  otimesMats (C_ijkl = A_ij B_kl)
%-----------------------------------------------------------------------------

function C = otimesMats(A,B)

[arows,acols] = size(A);
[brows,bcols] = size(B);

C = zeros(arows,acols,brows,bcols);

for i = 1:arows
    for j = 1:acols
        for k = 1:brows
            for l = 1:bcols
                C(i,j,k,l) = A(i,j) * B(k,l);
            end
        end
    end
end

end