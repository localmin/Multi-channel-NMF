% plot error wrt time graphs
% Eu
figure(1)
plot( wrt_f(:,1), wrt_i(:,2), '-r');
title('evaluate speed of algorithms on Frb');
xlabel('Time[s]');
ylabel('Frb-norm');
legend( 'MU' );

% IS
figure(2)
plot( wrt_i(:,1), wrt_i(:,2), '-r');
title('evaluate speed of algorithms on IS');
xlabel('Time[s]');
ylabel('IS-divergence');
legend( 'MU' );

% view last values of cost functions
printf( "last value[Frb] = %f\n", wrt_f(itr,2) );
printf( "last value[IS] = %f\n", wrt_i(itr,2) );

% view last values of relative error
vec_X = reshape( X, M*M*I*J, 1 );

xf_f = reshape( Xf_frb, M*M*I*J, 1 );
xf_i = reshape( Xf_is, M*M*I*J, 1 );

RE_frb = (xf_f - vec_X)' * (xf_f - vec_X ) / (vec_X' *vec_X ); 
RE_IS = (xf_i - vec_X)' * (xf_i - vec_X ) / (vec_X' *vec_X ); 

printf( "Relative_error[Frb] = %f\n", RE_frb );
printf( "Relative_error[IS] = %f\n", RE_IS );