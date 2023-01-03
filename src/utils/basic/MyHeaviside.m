function [y] = MyHeaviside(x,tol)

if x>= tol
    y = 1;
else
    y = 0;
end

end

