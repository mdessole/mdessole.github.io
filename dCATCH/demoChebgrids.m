% 11/06/2020 
% M. Dessole, F. Marcuzzi, M. Vianello

clear all;

% initializing regression degree
n = 2;

%number of variables 
dim = 5;

% initializing G-efficiency and max iterations
tol=0.99;
maxit=10000;

% enable prints
verbose = 1;

% generate initial measure support and its cartinality according parameters
[pts,m] = dCHEBNODES(n,4,dim,verbose);

% Lawson Hanson parameters
LHDM_options = struct( 'lsqnonneg', false, ... % NNLS is solved by Matlab's lsqnonneg when true, by LHDM otherwise
                       'init', false, ... % if true, initialization of Passive set via ULS is performed 
                       'k', ceil(nchoosek(2*n+dim,dim)/(n*(dim-1))), ... % parameter k in LHDM
                       'thres', 0.2222, ... % parameter thres in LHDM
                       'thres_w', 0.8 ); % parameter thtres_w in LHDM


% run test
[cpts,cw,geff,momerr]=dNORD(n,pts,tol,maxit,LHDM_options,verbose);


function [pts, m] = dCHEBNODES(n,k,dim, verbose)
    if verbose
        fprintf('**********************************\n');
        fprintf('%d-dim Chebyshev test, k=%d, n=%d \n', dim, k, n);
        fprintf('**********************************\n');
    end
    m = 2*k*n;

    pts1 = cos((2*(1:m)-1)*pi/(2*m));
    ptsin = cell(dim,1);
    for i=1:dim
        ptsin{i} = pts1;
    end
    ptsout = cell(dim,1);
    [ptsout{:}] = ndgrid(ptsin{:});
    M=m^dim;
    pts=zeros(M,dim);
    for i=1:dim
        pts(:,i) = reshape(ptsout{i},M,1);
    end
end
