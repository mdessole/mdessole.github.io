% 11/06/2020 
% M. Dessole, F. Marcuzzi, M. Vianello

clear all;

% initializing regression degree
n = 5;

% choose cardinality of the initial measure support
m = 100000;

%number of variables 
dim = 4;

% initializing G-efficiency and max iterations
tol=0.95;
maxit=10000;

% enable prints
verbose = 1;

% generate initial measure support
pts = haltonpts(m,dim,verbose);

% Lawson Hanson parameters
LHDM_options = struct( 'lsqnonneg', false, ... % NNLS is solved by Matlab's lsqnonneg when true, by LHDM otherwise
                       'init', false, ... % if true, initialization of Passive set via ULS is performed 
                       'k', ceil(nchoosek(2*n+dim,dim)/(n*(dim-1))), ... % parameter k in LHDM
                       'thres', 0.2222, ... % parameter thres in LHDM
                       'thres_w', 0.8 ); % parameter thtres_w in LHDM


% run test
[cpts,cw,geff,momerr]=dNORD(n,pts,tol,maxit,LHDM_options,verbose);


function [pts] = haltonpts(m,dim,verbose)
    if verbose
        fprintf('**********************************\n');
        fprintf('%d %d-dim Halton points test \n', m, dim);
        fprintf('**********************************\n');
    end
    pts = haltonseq(m,dim);
end
