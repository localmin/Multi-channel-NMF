clear

addpath('./Frb');

parameter_setting;

X = input_data( M, I, J, K );

sd = rand *100;
% Frb-ciriterion
randn( 'seed', sd ); 
[wrt_f Hf, Tf, Vf] = mnmf_Frb( X, K, itr );

% IS-divergence
randn( 'seed', sd ); 
[wrt_i Hi, Ti, Vi] = mnmf_IS( X, K, itr );
 
% make factorized models
Xh_frb = mk_fctrz( Hf, Tf, Vf  );
Xh_is = mk_fctrz( Hi, Ti, Vi );

view_result;