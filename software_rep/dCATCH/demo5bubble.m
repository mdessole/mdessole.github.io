% 11/06/2020
% M. Dessole, F. Marcuzzi, M. Vianello

clear all;

% initializing regression degree
n = 10;

%number of variables 
dim = 3;

% initializing G-efficiency and max iterations
tol=0.95;
maxit=10000;

%certer points of 5-bubble
X=[3 3 3; 1 1 1; -1 4 -1; 0 2 -3; -3 3 2.5];
%radii of 5-bubble
r = [2;3;2.5;1;1.5];
% generate initial measure support inside the 5-bubble
% starting from m Halton points
m = 40^3;
[pts] = multibubble(X,r,m,dim);


% Lawson Hanson parameters
LHDM_options = struct( 'lsqnonneg', false, ... % NNLS is solved by Matlab's lsqnonneg when true, by LHDM otherwise
                       'init', false, ... % if true, initialization of Passive set via ULS is performed 
                       'k', ceil(nchoosek(2*n+dim,dim)/(n*(dim-1))), ... % parameter k in LHDM
                       'thres', 0.2222, ... % parameter thres in LHDM
                       'thres_w', 0.8 ); % parameter thtres_w in LHDM


% run test
[cpts,cw,geff,momerr]=dNORD(n,pts,tol,maxit, LHDM_options);

function [pts] = multibubble(X,r,m,dim)

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

