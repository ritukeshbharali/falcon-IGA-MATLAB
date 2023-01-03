%-----------------------------------------------------------------------------
%  getSymI4 (symmetric fourth-order identity tensor)
%-----------------------------------------------------------------------------

function C = getSymI4(rnk)

C = zeros(rnk,rnk,rnk,rnk);

for i = 1:rnk
    for j = 1:rnk
        for k = 1:rnk
            for l = 1:rnk
                C(i,j,k,l) = 0.5*((i==k)&&(j==l)) + 0.5*((i==l)&&(j==k));
            end
        end
    end
end

end