function B = Quadix2SymMat(A,dim)

if dim == 2

    B = zeros(3,3);

    B(1,1) = A(1,1,1,1);
    B(2,1) = A(2,2,1,1);
    B(3,1) = 0.5 * ( A(1,2,1,1) + A(2,1,1,1) );

    B(1,2) = A(1,1,2,2);
    B(2,2) = A(2,2,2,2);
    B(3,2) = 0.5 * ( A(1,2,2,2) + A(2,1,2,2) );

    B(1,3) = 0.5 * ( A(1,1,1,2) + A(1,1,2,1) );
    B(2,3) = 0.5 * ( A(2,2,1,2) + A(2,2,2,1) );
    B(3,3) = 0.25 * ( A(1,2,1,2) + A(2,1,1,2) ...
                    + A(1,2,2,1) + A(2,1,2,1) );

else
    error('Not yet implemented!')
end