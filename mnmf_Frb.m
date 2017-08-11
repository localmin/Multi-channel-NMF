function [ wrt, H, T, V ] = mnmf_Frb( X, K, Itr )

% time counting start
tic;
exe_time = 0;

% Getting size
[M,~,I,J] = size( X );

% initialization
T = abs(randn( I, K ));
V = abs(randn( K, J ));
H = zeros( M, M, I, K );
Id = eye( M, M ) / M;

for ii=1:I
  for kk=1:K
    H(:,:,ii,kk) = Id;
  end
end

% Scale Fitting
Xf = make_fctrz( H, T, V );
up = X(:)' * Xf(:);
low = Xf(:)' * Xf(:);
coef = sqrt( up / low );
T = T * coef;
V = V * coef;
Xf = Xf * coef * coef;

% Iteration by MU
wrt = zeros( Itr, 2 );

for lp=1:Itr
 % create new tmp variables
  tmpT = T;
  tmpV = V;
  tmpH = H; 

  for k=1:K
    
    xhu = zeros( I, J );
    xhl = zeros( I, J );
    for i=1:I
      for j=1:J
        vh = vec( H(:,:,i,k) );
        vx = vec( X(:,:,i,j)' );
        vxf = vec( Xf(:,:,i,j)' );
        xhu(i,j) = ( vx' * vb );
        xhl(i,j) = ( vxf' * vb );
      end
    end
    
    % update T
    up = xhu * V(k,:)' ;
    low = xhl * V(k,:)' ;
    nT(:,k) = T(:,k) .* (up ./ low);

    % update V
    up = T(:,k)' * xhu;
    low = T(:,k)' * xhl;
    tmpxnV(k,:) = V(k,:) .* (up ./ low);

    % update H
    for i=1:I
      
      A = zeros(M,M);
      B = zeros(M,M); 
      for j=1:J
        tmp1 = V(k,j) * Xf(:,:,i,j);
        A = A + tmp1;

        tmp2 = V(k,j) * X(:,:,i,j);
        B = B + tmp2;
      end
      tmp3 = inv(A);
     
      tmp = H(:,:,i,k) * tmp3 * B;
      tmp = ( tmp + tmp' ) / 2;
      [P,L] = eig( tmp );
      tmp = real( P * diag( max( diag(L), 0 ) ) * P' );
      tmp = tmp / trace( tmp );
      tmpH(:,:,i,k) = tmp;
      
    end
   
  end
  % Update variables
  H = tmpH;
  V = tmpV; 
  T = tmpT;
  Xf = make_fctrz( H, T, V );
  
  exe_time = exe_time + toc;

  D = reshape( X - Xf, M*M*I*J,1);
  % Extension Euclidean
  wrt(lp,:) = [ exe_time D'*D ];

  tic;
  
end

end