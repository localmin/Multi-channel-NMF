function [ H ] = Riccati( A, B )

    [M,~] = size( A );
    O = zeros( M, M );
    A = (A + A')/2;
    B = (B + B')/2;
    
    %Heamiltonian matrix
    Z = [ O -A ; -B O ];
    [P,D] = eig(Z);
    [~,indx] = sort( real(diag( D )) );
    P = P(:,indx);
    
    F = P(1:M,1:M);
    G = P(M+1:2*M,1:M);
    
    H = G * pinv(F);
    
end
