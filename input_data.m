function X = input_data( M, I, J, K )

H = zeros( M, M, I, K );
T = abs(randn( I, K ));
V = abs(randn( K, J ));

for i=1:I
  for k=1:K
    temp = randn( M, M ); temp = temp * temp';
    temp = temp / trace(temp);
    H(:,:,i,k) = temp;
  end
end

X = make_fctrz( H, T, V );

end


