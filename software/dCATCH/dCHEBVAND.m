function V = dCHEBVAND(n,X)

% computes the Chebyshev-Vandermonde matrix for degree n on a d-dim mesh X

% input
% n: polynomial degree
% X: d-column array of mesh points

% output
% V: Chebyshev-Vandermonde matrix for degree n on the mesh X  

% 11/06/2020 
% M. Dessole, F. Marcuzzi, M. Vianello

% FUNCTION BODY

% box containing the mesh
d=size(X,2); 
a=min(X);
b=max(X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d-uples of indices with sum less or equal to n
% graded lexicographical order

N = nchoosek(n+d,d);
duples = zeros(N,d);
for i=2:N
    duples(i,:) = mono_next_grlex(d,duples(i-1,:));
end

% mapping the mesh in the hypercube [-1,1]^d
map = zeros(size(X));
for i=1:d
    map(:,i)=(2*X(:,i)-b(i)-a(i))/(b(i)-a(i));
end

% Chebyshev-Vandermonde matrix on the mesh
T=chebpolys(n,map(:,1));
V=T(:,duples(:,1)+1);
for i=2:d
    T=chebpolys(n,map(:,i));
    V=V.*T(:,duples(:,i)+1);
end
end


function T=chebpolys(n,x)

% computes the univariate Chebyshev-Vandermonde matrix by recurrence 
% for degree n on an array of 1d points

% input
% n: polynomial degree
% x: 1-column array of 1d points

% output
% T: Chebyshev-Vandermonde matrix for degree n on x

% FUNCTION BODY

T=zeros(length(x),n+1);
t0=ones(length(x),1);
T(:,1)=t0;
t1=x;
T(:,2)=t1;

for j=2:n
    t2=2*x.*t1-t0;
    T(:,j+1)=t2;
    t0=t1;
    t1=t2;
end
end


