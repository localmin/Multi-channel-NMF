# Multi-channel Non-negative Matrix Factorization

Multi-channel Non-negative Matrix Factorization by multiplicative update written in Octave.

## Features

You can execute Multi-channel NMF easily on Octave(of course MATLAB).

## Background of Multi-channel NMF

You should read the following article,
[http://www.kecl.ntt.co.jp/people/kameoka.hirokazu/publications/Sawada2013IEEETrans05-published.pdf](http://www.kecl.ntt.co.jp/people/kameoka.hirokazu/publications/Sawada2013IEEETrans05-published.pdf)
(I'm not an author of Multi-channel NMF.)

## Usage
`./exc_mnmf.m`

 ### Default Setting
  
 * This default program execute a simulation by real number. 
 
 If you want to deal with real data, you can change 
 codes, for example, "parameter_setting.m", "input_data.m" etc., for your use.
 
  * Size of the observation matrix: 12 × 20
  
  * Channel number: 2
  
  * Numbers of basis vectors: 4
 
 * Numbers of iterations: 100
 
 * Initial value setting method of update rules: random non-negative values for T & V, and identity matrix for H.
 
 You can change these default settings on "parameter_setting.m" like the following,

 

 ```
% parameters
% Number of channel
M = 2;
% size of observation matrix
I = 12;
J = 20;
% basis vectors
K = 2; 
% numbers of iterations
itr = 100;
 ```
 
 Data that you want to use can be also changed on "input_data.m"
 
 Initial value setting can be changed on each functions, "mnmf_Frb.m& mnmf_IS.m"
 
 ```
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
 ```

* I do NOT implement the NMF basis clustering and top-down clustering using the latent variables Z and bottom-up clustering in the reference paper.

 ## Caution

 Please read [Issue.](https://github.com/localmin/Multi-channel-NMF/issues)


 This code is provided without liability and warranty.
