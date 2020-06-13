function [pts,w, momerr] = dCATCH(deg,X,u,LHDM_options)

% Caratheodory-Tchakaloff d-variate discrete measure compression  
% for example probability measures (designs) or quadrature formulas;
% the moments are invariant (close to machine precision) up to degree deg;  
% adapts to the (numerical) dim of the polynomial space on the point set X; 
% works satisfactorily for low/moderate degrees

% INPUT:   
% deg: polynomial exactness degree
% X: d-column array of point coordinates
% u: 1-column array of positive weights
% LHDM_options: structure containing the values of optimization parameters

% OUTPUT:
% pts: d-column array of extracted mass points coordinates
% w: 1-column array of corresponding new positive weights 
% momerr: moment reconstruction error

% 11/06/2020 
% M. Dessole, F. Marcuzzi, M. Vianello

% FUNCTION BODY 


% Vandermonde-like matrix of a u-orthogonal polynomial basis on X  
U =dORTHVAND(deg,X,u);

fprintf("Vandermonde matrix's size = %d x %d \n", size(U,1), size(U,2));

%if length(U(1,:))<=length(X(1,:))
if size(U,1)<=size(U,2)
 
    fprintf("Vandermonde matrix is not underdetermined, nothing to compress \n");
    
    % no compression expected  
    pts=X;
    w=u;
    momerr = 0;
else

    % further orthogonalization to reduce the conditioning  
    [Q,~]=qr(U,0);

    % new moments 
    orthmom=Q'*u;
    [pts, w, momerr, ~]= NNLS(X, u, U, Q, orthmom, LHDM_options);

end

end

function [pts, w, momerr, e]= NNLS(X, u, U, Q, orthmom, options)
    % Caratheodory-Tchakaloff points and weights via accelerated NNLS 
    tic;
    if  isfield(options,'lsqnonneg')
        if options.lsqnonneg
            fprintf("Matlab's lsqnonneg \n");
            [weights,~,~,~,output] = lsqnonneg(Q',orthmom); 
            iter = output.iterations;
            cardP = [];
        else
            [weights,~,~,iter]=LHDM(Q',orthmom,options);
        end
    else
        [weights,~,~,iter]=LHDM(Q',orthmom,options);
    end   
    e = toc;
    
    % indexes of nonvanishing weights and compression  
    ind=find(abs(weights)>0);
    pts=X(ind,:);
    w=weights(ind);
    
    % moment reconstruction error  
    momerr=norm(U(ind,:)'*w-U'*u);

    % displaying results
    fprintf('number of outer iterations = %d \n', iter);
    fprintf('elapsed time = %.6f \n', e);
    fprintf('initial design cardinality = %4.0f \n',size(X,1));
    fprintf('concentrated support cardinality = %4.0f \n',length(w));
    fprintf('compression ratio = %4.0f \n',size(X,1)/length(w));
    fprintf('moment reconstruction error = %4.2e \n \n',momerr);
    
end
