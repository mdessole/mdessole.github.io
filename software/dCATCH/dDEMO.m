% 11/06/2020
clear all;

% initializing regression degree
n = 8;

% choose nb of discretization points per axis
m = 100;

%number of variables 
dim = 10;

% initializing G-efficiency and max iterations
tol=0.95;
maxit=10000;

% enable prints
verbose = 1;

% generate initial measure support according to the chosen example 
example = 1;
if (example == 1)
    % dim is set to 2
    [pts, dim] = france_setup(m);
elseif (example == 1)
    pts = haltonpts(m,dim);
elseif (example == 1)
    % tol is set to 0.99
    % m is computed according to example's parameters
    [pts, tol, m] = dCHEBNODES(n,4,dim);
elseif (example == 1)
    %certer points of 5-bubble
    X=[3 3 3;1 1 1; -1 4 -1; 0 2 -3; -3 3 2.5];
    %radii of 5-bubble
    r = [2;3;2.5;1;1.5];
    %generate Halton points inside the 5-bubble
    [pts] = multiball(X,r,m,dim);
else
    fprintf('Example choice %d not implemented. Running Halton points test...', example);
end

% Lawson Hanson parameters
LHDM_options = struct( 'lsqnonneg', false, ... % NNLS is solved by Matlab's lsqnonneg when true, by LHDM otherwise
                       'init', false, ... % if true, initialization of Passive set via ULS is performed 
                       'k', ceil(nchoosek(2*n+dim,dim)/(n*(dim-1))), ... % parameter k in LHDM
                       'thres', 0.2222, ... % parameter thres in LHDM
                       'thres_w', 0.8 ); % parameter thtres_w in LHDM


% run test
[cpts,cw,geff,momerr]=dNORD(n,pts,tol,maxit, LHDM_options, verbose);

function [pts] = multiball(X,r,m,dim)

    fprintf("**********************************\n");
    fprintf("%d %d-dim multiball test \n", m, dim);
    fprintf("**********************************\n");

    pts = haltonseq(m,dim);
    %computing the minimal box
    a = min(X - repmat(r,1,dim),[],1); % d-dim array
    b = max(X + repmat(r,1,dim),[],1); % d-dim array
    %map halton points into the d-dim box [a,b]
    for i=1:dim
        pts(:,i) = a(i) + pts(:,i)*(b(i)-a(i));
    end
    %map halton points into the d-dim box [a,b]
    cond = false(m,1);
    for i=1:size(X,1)
        sum = zeros(m,1);
        for j=1:dim
            sum = sum + (pts(:,j) - X(i,j)).^2;
        end
        cond = cond | (sum < r(i)^2);
    end
    pts(find(~cond),:) = [];
end

function [pts, tol, m] = dCHEBNODES(n,k,dim)

    fprintf("**********************************\n");
    fprintf("%d-dim Chebyshev test, k=%d, n=%d \n", dim, k, n);
    fprintf("**********************************\n");

    m = 2*k*n;
    tol = 0.99;
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


function [pts] = haltonpts(m,dim)

    fprintf("**********************************\n");
    fprintf("%d %d-dim Halton points test \n", m, dim);
    fprintf("**********************************\n");
    test_name = strcat('halton_d', num2str(dim));
    pts = haltonseq(m,dim);
end

function [pts, dim] = france_setup(m)

    fprintf("***********************************\n");
    fprintf("France test on a  %d x %d grid\n", m, m);
    fprintf("***********************************\n");

    dim = 2;

    % EXAMPLE: France-shaped polygon
    Q=[0 0;
        8.2 -2;
        8 0;10 0.8;13 -1;15.5 1;14 4;15 4.8;14.5 7.2;13 6.5; ...
        15 10;16 10.3;17 14.5;12.5 15.7;7.2 19.6;6 19.2;6 17.5;3 16.1; ...
        3.2 15.7; 1.1 15.5; 1 16.1; 0 16.1;0.3 13.7;-4 13.8;-4.5 12; ...
        -1 10.5; 1 8;0 0];

    
    % minimal rectangle containing Q  
    %m=100;
    a=min(Q(:,1)); b=max(Q(:,1));
    c=min(Q(:,2)); d=max(Q(:,2));
    q1=(a+b)/2;  
    q2=(c+d)/2;
    r=max(b-a,c-d)/2;
    
    % intersecting an mxm grid with Q 
    u=linspace(q1-r,q1+r,m);
    v=linspace(q2-r,q2+r,m);
    [xx,yy]=meshgrid(u,v);
    p=[xx(:) yy(:)];
    ind=find(inpolygon(p(:,1),p(:,2),Q(:,1),Q(:,2))); %True if point is in polygon
    pts=p(ind,:); %points of the grid lying in polygon
    
end

