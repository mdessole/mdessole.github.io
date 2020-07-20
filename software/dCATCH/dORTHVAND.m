function [U,jvec,Q,R] = dORTHVAND(deg,X,u,jvec,C)

% computes a Vandermonde-like matrix for degree n on a d-dim point set X 
% in an total-degree discrete orthogonal polynomial basis w.r.t. 
% the nonnegative weight array u 

% input
% deg: polynomial degree
% X: d-column array of point coordinates
% u: 1-column array of nonnegative weights, or nonnegative scalar in
%    case of equal weithgs
% jvec (optional): vector of column indexes, selects a polynomial basis  
% C (optional): Chebyshev-Vandermonde matrix on jvec basis

% output
% U: Vandermode-like matrix in a u-orthogonal polynomial basis on X  
% jvec: vector of column indexes, selects a polynomial basis (computed  
%       if the input jvec is empty)
% Q: orthogonal factor in the QR decomposition
%    diag(sqrt(u))*C(:,jvec)=Q*R where C=dCHEBVAND(n,X)
% R: triangular factor in the QR decomposition
%    diag(sqrt(u))*C(:,jvec)=Q*R where C=dCHEBVAND(n,X)

% 11/06/2020 
% M. Dessole, F. Marcuzzi, M. Vianello

% FUNCTION BODY


% total-degree Chebyshev-Vandermonde matrix on X and
% dimension of the polynomial space on X  
if (nargin == 4)
    C=[];
elseif (nargin == 3)
    C=[];
    jvec = [];
end

if isempty(C) 
    C=dCHEBVAND(deg,X);
end
if isempty(jvec)
    N=rank(C);
else
    N=length(jvec); 
end

% scaling the matrix rows by the sqrt of the weights 
B = zeros(size(C));
if isscalar(u)
    B=C*sqrt(u);
else
    for k=1:length(C(1,:))
        B(:,k)=C(:,k).*sqrt(u);
    end
end


% polynomial basis orthogonalization
if N<length(C(1,:)) 
    if isempty(jvec)
        try
           [Q,R0,pm]=qr(B,'vector');
        catch ME
            fprinf('Q is too large. Computing economy-size QR ...\n');
            [Q,R0,pm]=qr(B,0);%'vector');
        end 
        jvec=pm(1:N);
        R=R0(1:N,1:N);
    else
        [Q,R]=qr(B(:,jvec),0);
    end
else 
    [Q,R]=qr(B,0);
    jvec=(1:N);
end

U=C(:,jvec)/R;


end

