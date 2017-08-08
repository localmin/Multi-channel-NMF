function X = input_data( M, R, C, K )

H = zeros( M, M, R, K );
T = abs(randn( R, K ));
V = abs(randn( K, C ));

for i=1:R
  for k=1:K
    temp = randn( M, M ); temp = temp * temp';
    temp = temp / trace(temp);
    H(:,:,i,k) = temp;
  end
end

X = prodHTV( H, T, V );

end


