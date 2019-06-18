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
function [groups, kerN] = SpectralClustering_SPECMC(CKSym,n, weight, idx0,opt)
if (nargin < 5)
    opt =0; % unnormalization for L
end
if nargin < 4
    warmstart =0;
else
    warmstart =1;
end
warning off;
%N = size(CKSym,1);
N = size(CKSym{1},1);


MAXiter = 1000; % Maximum number of iterations for KMeans 
%REPlic = 20; % Number of replications for KMeans
REPlic = 10; % Number of replications for KMeans

% Normalized spectral clustering according to Ng & Jordan & Weiss
% using Normalized Symmetric Laplacian L = I - D^{-1/2} W D^{-1/2}
if (~opt)
     V= size(CKSym,2);
    
    for v=1:V
        DN{v} = diag( 1./sqrt(sum(CKSym{v})+eps) );
        LaplaceN{v} = speye(N) - DN{v} * CKSym{v} * DN{v};
    end
    LapN = 0;
    for v=1:V
        LapN = LapN + weight(v)*LaplaceN{v};
    end
    
   %LapN = LapN/V;

    
%     DN = diag( 1./sqrt(sum(CKSym)+eps) );
%     LapN = speye(N) - DN * CKSym * DN;
    [~,~,vN] = svd(LapN);
    kerN = vN(:,N-n+1:N);
    for i = 1:N
        kerNS(i,:) = kerN(i,:) ./ norm(kerN(i,:)+eps);
    end
    
else
    %%  Unnormalized Laplacian
    V= size(CKSym,2);
    Lap = 0;
    for v=1:V
        D{v} = diag(sum(CKSym{v}));
        Lapalace{v} = D{v} - CKSym{v};
        Lap = Lap + Lapalace{v};
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
else
    [groups]= kmeans(kerNS,n,'start', 'warmstart', 'initiallabel',idx0, 'maxiter',MAXiter,'replicates',REPlic,'emptyaction','singleton');
end