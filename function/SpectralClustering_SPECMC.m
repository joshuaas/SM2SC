% SpectralClustering.m
% This is a modification of the ordinary spectral clustering which allows the initialization from previous clustering
% results.
%
%--------------------------------------------------------------------------
% This function takes an adjacency matrix of a graph and computes the 
% clustering of the nodes using the spectral clustering algorithm of 
% Ng, Jordan and Weiss.
% CMat: NxN adjacency matrix
% n: number of groups for clustering
% groups: N-dimensional vector containing the memberships of the N points 
% to the n groups obtained by spectral clustering
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
%--------------------------------------------------------------------------
% Modified by Chunguang Li, Nov. 7, 2014.
function [groups, kerNS] = SpectralClustering_SPECMC(CKSym,n,idx0 ,opt)
if (nargin < 4)
    opt =0; % unnormalization for L
end
if nargin < 3
    warmstart =0;
else
    warmstart =1;
end
warning off;
%N = size(CKSym,1);
N = size(CKSym{1},1);


MAXiter = 1000; % Maximum number of iterations for KMeans 
%REPlic = 20; % Number of replications for KMeans
REPlic = 25; % Number of replications for KMeans

% Normalized spectral clustering according to Ng & Jordan & Weiss
% using Normalized Symmetric Laplacian L = I - D^{-1/2} W D^{-1/2}
if (~opt)
     V= numel(CKSym);
     DN = cell(V,1);
     LaplaceN = cell(V,1);
     
    for v=1:V
        DN{v} = diag( 1./sqrt(sum(CKSym{v})+eps) );
        LaplaceN{v} = speye(N) - DN{v} * CKSym{v} * DN{v};
    end
    
    LapN = 0;
    
    for v=1:V
        LapN = LapN +  LaplaceN{v};
    end
    
   %LapN = LapN/V;
    
%     DN = diag( 1./sqrt(sum(CKSym)+eps) );
%     LapN = speye(N) - DN * CKSym * DN;
    [~, ~,vN] = svd(LapN);
    kerN = vN(:, N-n+1:N );
    kerNS = zeros( size(kerN) );
    
    for i = 1:N
        kerNS(i,:) = kerN(i,:) ./ norm(kerN(i,:)+eps);
    end
    
else
    %%  Unnormalized Laplacian
    V= numel(CKSym);
    Lap = 0;
    for v=1:V
        D{v} = diag(sum(CKSym{v}));
        Lapalace{v} = D{v} - CKSym{v};
        Lap = Lap +  Lapalace{v};
    end

%     D =diag(sum(CKSym));
%     Lap = D - CKSym;
    
    %% SVD
    
    [~,~,V] = svd(Lap);
    kerN = V(:,N-n+1:N);
    
%     %% Re-normalization
%     for i = 1:N
%         kerNS(i,:) = kerN(i,:) ./ norm(kerN(i,:)+eps);
%     end    
    kerNS =kerN;
end

 if (~warmstart)
    [groups]  = kmeans(kerNS,n,'start', 'sample', 'maxiter',MAXiter,'replicates',REPlic,'emptyaction','singleton');
%     s = silhouette(kerNS, groups);
%     s = sum(s) ;
%     if s > gbest.s
%         gbest.groups = groups;
%         gbest.s      = s ;
%     end
else
    [groups,~,loss]= kmeans(kerNS,n,'start', 'warmstart', 'initiallabel',idx0, 'maxiter',MAXiter,'replicates',800,'emptyaction','singleton');
%     s = silhouette(kerNS, groups);
%     s = sum(loss) ;
%     if s < gbest.s
%         gbest.groups = groups;
%         gbest.s      = s ;
%     end
end