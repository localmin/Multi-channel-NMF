function X = make_fctrz( H, T, V )

[ I, J ] = size( T * V );
K = size( T, 2 );
M = size( H, 1 ); %H:(i,k)

X = zeros( M, M, R, C );
for i=1:R
  for j=1:C
    for k=1:K
      X(:,:,i,j) = X(:,:,i,j) + H(:,:,i,k) * T(i,k) * V(k,j);
    end
  end
end

end

