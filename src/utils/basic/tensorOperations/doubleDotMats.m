%-----------------------------------------------------------------------------
%  doubleDotMats (c = A_ij B_ij)
%-----------------------------------------------------------------------------

function c = doubleDotMats(A,B)

c = 0.0;

[arows,acols] = size(A);

for i = 1:arows
    for j = 1:acols
        c = c + A(i,j) * B(i,j);
    end
end

end