% 11/06/2020
function [U,jvec,Q,R] = dORTHVAND(deg,X,u,jvec,C)

% computes a Vandermonde-like matrix for degree n on a d-dim point set X 
% in an total-degree discrete orthogonal polynomial basis w.r.t. 
% the nonnegative weight array u 

% input
% deg: polynomial degree
% X: d-column array of point coordinates 
% u: 1-column array of nonnegative weights 
% jvec (optional): vector of column indexes, selects a polynomial basis  
% C (optional): Chebyshev-Vandermonde matrix on jvec basis

% output
% U: Vandermode-like matrix in a u-orthogonal polynomial basis on X  
% jvec: vector of column indexes, selects a polynomial basis (computed  
% if the input jvec is empty)
% Q:
% R:

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
    C=dchebvand(deg,X);
end
if isempty(jvec)
    N=rank(C);
else
    N=length(jvec); 
end

% scaling the matrix rows by the sqrt of the weights 
B = zeros(size(C));
for k=1:length(C(1,:))
    B(:,k)=C(:,k).*sqrt(u);
end

% polynomial basis orthogonalization
if N<length(C(1,:)) 
    if isempty(jvec)
        try
           [Q,R0,pm]=qr(B,'vector');
        catch ME
            fprinf("Q is too large. Computing economy-size QR ...\n");
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

% function [U,jvec,Q,R] = dORTHVAND(deg,X,u,varargin)
% 
% % computes a Vandermonde-like matrix for degree n on a d-dim point set X 
% % in an total-degree discrete orthogonal polynomial basis w.r.t. 
% % the nonnegative weight array u 
% 
% % input
% % deg: polynomial degree
% % X: d-column array of point coordinates 
% % u: 1-column array of nonnegative weights 
% % jvec: vector of column indexes, selects a polynomial basis  
% 
% % output
% % U: Vandermode-like matrix in a u-orthogonal polynomial basis on X  
% % jvec: vector of column indexes, selects a polynomial basis (computed  
% % if the input jvec is empty)
% 
% % FUNCTION BODY
% % total-degree Chebyshev-Vandermonde matrix on X and
% % dimension of the polynomial space on X  
% 
% global opt_orth
% 
% if (nargin == 5)
%     C = varargin{5}; 
%     jvec = varargin{4};
% elseif (nargin == 4)
%     jvec = varargin{4};
%     C=[];
% else
%     C=[];
%     jvec = [];
% end
% if isempty(C) 
%     C=dchebvand(deg,X);
% end
% if isempty(jec)
%     N=rank(C);
% else
%     N=length(jvec); 
% end
% 
% % scaling the matrix rows by the sqrt of the weights 
% for k=1:length(C(1,:))
%     B(:,k)=C(:,k).*sqrt(u);
% end
% 
% % polynomial basis orthogonalization
% if N<length(C(1,:)) 
%     if isempty(jvec)
%         [Q,R0,pm]=qr(B,'vector');
%         jvec=pm(1:N);
%         R=R0(1:N,1:N);
%     else
%         [Q,R]=qr(B(:,jvec),0);
%     end
% else 
%     [Q,R]=qr(B,0);
%     jvec=(1:N);
% end
% 
% if (opt_orth == 1)
%     U = Q;
% else
%     U=C(:,jvec)/R;
% end
% 
% end 
% 
% 
