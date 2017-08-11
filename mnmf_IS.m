function [ wrt, H, T, V ] = mnmf_IS( X, K, Itr )

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
upper = X(:)' * Xf(:);
lower = Xf(:)' * Xf(:);
coef = sqrt( up / low );
T = T * coef;
V = V * coef;
Xf = Xf * coef * coef;

% Itration by IS
wrt = zeros( Itr, 2 );

Xfp = zeros( size( Xf ) );
for i=1:I
  for j=1:J
    Xfp(:,:,i,j) = pinv(Xf(:,:,i,j) );
  end
end

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
        vp = vec( Xfp(:,:,i,j)' );
        vxph = vec( X(:,:,i,j) * Xfp(:,:,i,j) * H(:,:,i,k) );
        vh = vec( H(:,:,i,k) );
        tru(i,j) = ( vp' * vxph );
        trl(i,j) = ( vp' * vh );
      end
    end
    
    % update T
    up = tru * V(k,:)' ;
    low = trl * V(k,:)' ;
    tmpT(:,k) = T(:,k) .* sqrt(up ./ low);

    % update V
    up = T(:,k)' * tru;
    low = T(:,k)' * trl;
    tmpV(k,:) = V(k,:) .* sqrt(up ./ low);

    % update H
    for i=1:I
      
      % update H
      A = zeros(M,M);
      BB = zeros(M,M);
      for j=1:J
        tmp = V(k,j) * Xfp(:,:,i,j);
        A = A + tmp;
        BB = BB + tmp * X(:,:,i,j) * Xfp(:,:,i,j);
      end
      B = H(:,:,i,k) * BB * H(:,:,i,k);
      
      %function of riccati
      tmp = Riccati(A,B);
      tmp = ( tmp + tmp' ) / 2;
      tmp = tmp / trace( tmp );
      tmpH(:,:,i,k) = tmp;
      
    end
   
  end
  % Update variables
  H = tmpH; 
  V = tmpV;
  T = tmpT;
  Xf = make_fctrz( H, T, V );

  for i=1:I
    for j=1:J
      Xfp(:,:,i,j) = pinv(Xf(:,:,i,j) );
    end
  end
  exe_time = exe_time + toc;

  
  % make IS
  IS = zeros( I, J );
  for i=1:I
    for j=1:J
        t = X(:,:,i,j) * Xfp(:,:,i,j);
        IS(i,j) = trace(t) - log(det(t)) - M;
    end
  end
 
  D = reshape( IS, R*C,1 );
  % Extension IS
  wrt(lp,:) = [ exe_time norm(D,1) ];

  tic;
  
end

end