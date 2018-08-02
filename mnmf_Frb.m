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
    
      vX = zeros(M,M,I,K);
		  vXf = zeros(M,M,I,K);
		  upperT = zeros(I,K);
	  	lowerT = zeros(I,K);

		  % Initialize T & update H
		  for i = 1:I
			  for k = 1:K

				  for j = 1:J
					  vX(:,:,i,k) =  vX(:,:,i,k) + V(k,j) * X(:,:,i,j);
					  vXf(:,:,i,k) =  vXf(:,:,i,k) + V(k,j) * Xf(:,:,i,j);
				  end
				  upperT(i,k) = vec(H(:,:,i,k))'* vec(vX(:,:,i,k));
				  lowerT(i,k) = vec(H(:,:,i,k))'* vec(vXf(:,:,i,k));

			  end
		  end	

	    % update T 
	    tmpT = T.*(upperT./lowerT);

   	    % InitializeV 
		  upperV = zeros(K,J);
		  lowerV = zeros(K,J);

		  for j = 1:J
			  for k = 1:K
				  for i = 1:I
					  upperV(k,j) = upperV(k,j) + T(i,k) * vec(H(:,:,i,k))'* vec(X(:,:,i,j));
					  lowerV(k,j) = lowerV(k,j) + T(i,k) * vec(H(:,:,i,k))'* vec(Xf(:,:,i,j));
				  end
			  end
		  end	 

	  	% update V
      tmpV = V.*(upperV./lowerV);

		  % H initialozed
		  % update H(Not  refactored)
		  for i = 1:I
			  for k = 1:K
            vXfi = inv(vXf(:,:,i,k));

	      		tmp = H(:,:,i,k) *  vXfi * vX(:,:,i,k);
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
