function [cpts,cw,geff,momerr] = dNORD(deg,X,gefftol,maxit,LHDM_options,verbose)

% computes a Near Optimal Regression Design with a given G-efficiency on   
% a discrete d-dim set X, by a basic multiplicative algorithm and 
% Caratheodory-Tchakaloff discrete measure compression  

% INPUT
% deg: polynomial regression degree
% X: d-column array of set points coordinates
% gefftol: G-efficiency threshold 
% maxit: maximum number of iterations of the multiplicative algorithm 
% LHDM_options: structure containing the values of optimization parameters

% OUTPUT
% cpts,cw: compressed near optimal design points and weights 
% geff: G-efficiency of the output design 
% momerr: moment reconstruction error by the compressed measure

% 11/06/2020 
% M. Dessole, F. Marcuzzi, M. Vianello

% FUNCTION BODY 

% initializing the probability weights
M=length(X(:,1));
w=ones(M,1)/M;

nit=0; go=1;

tic;
% multiplicative iteration up to the given G-efficiency 
while go==1
    nit=nit+1;
    if nit==1 
        [U,jvec]=dORTHVAND(deg,X,w);
        % dimension of the polynomial space on X 
        rk=length(jvec);
    else
        U=dORTHVAND(deg,X,w,jvec);
    end
    % updating the design  
    N=rk;
    K=sum((U.*U)');
    % computing G-efficiency of the design  
    geff=N/max(K);
    % updating the cicle exit condition
    go=(geff<gefftol & nit<maxit);
    if go==1
        % updating the design
        w=w.*K'/N;
    end
end
elapsed = toc;
fprintf('Titterington iterations = %d, elapsed time = %.6f s \n', nit, elapsed);

% Caratheodory-Tchakaloff design compression 
[cpts,cw,momerr]=dCATCH(2*deg,X,w,LHDM_options,verbose);

end

