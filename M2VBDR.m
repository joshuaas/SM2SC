function [CD, CDmean,Z, Zmean, C,loss,dlambda ] = M2VBDR(X,gt,param,extparam) 

lambda1 = param(1) ;
lambda2 = param(2) ;
lambda3 = param(3) ;
thresh1 = param(4) ;
thresh2 = param(5) ;
rseed = -1 ;
ex_flag= 0;
if nargin > 3
    
    if length(extparam) == 2
        ex_flag  = 1 ;
        D_min  = extparam(1) ;
        C_max  = extparam(2) ;
    else
        rseed = extparam ; 
    end
end
% min_{Z,B,W} 0.5*||X-XZ||_F^2+0.5*lambda*||Z-B||_F^2+gamma*<Diag(B1)-B,W>
% s.t. diag(B)=0, B>=0, B=B^T, 
%      0<=W<=I, Tr(W)=k.
k = max(unique(gt)) ;
if nargin < 5
    display = 0;
end
N = cellfun(@(x) size(x,2), X,'UniformOutput', false);
tol = 1e-3;
maxIter = 10; 
V = numel(X) ;
XtX = cellfun(@(x) x'*x, X,'UniformOutput', false);
invXtXI = cellfun(@(x,n) inv( x + lambda1 * eye(n) ) , XtX , N,'UniformOutput', false);
XtXprod  = cellfun(@(x,y) x * y, invXtXI,XtX, 'UniformOutput', false ) ;
loss = zeros(maxIter, 1) ; 
if rseed >  0
 Z = cellfun(@(n)rseed * randn(n),N, 'UniformOutput', false) ;
else
 Z = cellfun(@(n)zeros(n),N, 'UniformOutput', false) ;
end
% Z = CellFunc( @(x)constructW_PKN(x),X) ;
Zk = Z ;
W = eye(N{1});
% W = Z{1} ;
C =  ones(N{1});
C  = C - diag(diag(C)) ;
L = CalLaplacian(C);    
dlambda = [] ;
D = cellfun(@(n) zeros(n),N, 'UniformOutput', false) ;
iter = 0;
diffZ=zeros(maxIter,1) ;
diffC = diffZ ;
diffD = zeros(maxIter,1) ;
diffW = diffD;
loss = zeros(maxIter,1);
while iter < maxIter
    iter = iter + 1;        
     
% W = VV * VV';
    


    % update Z
    Zk = Z ;
   
%     CDv = CellFunc(@(x)C .* x, D ) ;
%     if(iter > 1)
%        for ii = 1:V
%             [Ud,sd,Vd] = svdecon(CDv{ii}) ;
%             Deld{ii} = (invXtXI{ii} * Ud) * sd * Vd' ; 
%        end
%     end

    

    Z = cellfun(@(x,y,z) x +  lambda1 *  y* (C.*z) , XtXprod, invXtXI, D,'UniformOutput', false) ;
    
    % update D
    Dk = D ;
    D = cellfun(@(z) (lambda1 *C.* z) ./ (lambda1 *C.^2 + lambda2) ,Z,'UniformOutput', false) ;
    
    % update C
    if ex_flag
    	D = cellfun(@(x)max(x, D_min),D, 'UniformOutput', false) ;
    else
    	D = cellfun(@(x)max(x, 1e-4),D, 'UniformOutput', false) ;
    end
% D = CellFunc(@(x) ProjD(x),D) ;
    Ck = C;

%     tW =  diag(W) * ones(1,N{1})-W ;
%      for i =1:N{1}
%          for j = (i+1):N{1}
%              d1 = 0;d2 =0;
%                 for v = 1:V
%                          d1 = d1 + D{v}(i,j)^2 + D{v}(j,i)^2 ;
%                          d2 = d2 + D{v}(i,j)*Z{v}(i,j) + D{v}(j,i)*Z{v}(j,i) ;
%                 end 
%                 d3 = tW(i,j) + tW(j,i) ;
%                 C(i,j) = (lambda1 * d2 - lambda3 * d3) ./(lambda1 * d1 +0.01 ) ;
%                 C(j,i)  = C(i,j) ;
%          end 
%      end
  
S1 = zeros(size(C)) ;
S2 = zeros(size(C)) ;
for v = 1:V
    S1 = S1 + D{v}.* Z{v} + (D{v}') .* (Z{v}');
    S2 = S2 + D{v}.^2 + (D{v}').^2;
end
tW =  diag(W) * ones(1,N{1})-W ;
tW = tW + tW';
C  = (lambda1 * S1 - lambda3 *tW) ./(lambda1 * S2+thresh1) ;
C = max(thresh2, (C+C')/2);
if ex_flag
	C = min(C, C_max) ;
else
	C = min(C,10) ;
end
C = C - diag(diag(C));
L = CalLaplacian(C);  
 Wold = W ;
 if iter > 1
    DDold = DD ;
 end
%      nonz = [ones(N{1} - k,1) ; zeros(k,1)] ;
%      [VV, DD] = eig(L + 0.01 *diag(nonz));
     [VV, DD] = eig(L );

%    [VV,~] = laneig(L,k,'AS') ;
    DD = diag(DD);
     [~, ind] = sort(DD,'ascend'); 
%     [~, k ] = max(diff(DD)) ; 
W = VV(:,ind(1:k)) * VV(:,ind(1:k))';       
    % update W
    CaldiffZ ;
    CaldiffC ;
    CaldiffD ;
    CaldiffW ;

    stopC = max([diffZ(iter),diffC(iter), diffD(iter)]);
  CalLoss ;    
    if stopC < tol 
        break;
    end
end
diff.Z = diffZ ;
diff.C = diffC ;
diff.D = diffD ;
diff.W = diffW;
CD   = CellFunc(@(x) C.* x, D) ;
CD   = CellFunc(@(x) (x + x')/2, CD ) ;
Z     = CellFunc(@(z) (z + z')/2, Z) ;
CDmean = MergeAdjacentMatrix(CD) ;
Zmean = MergeAdjacentMatrix(Z) ;
 
%     function  Dp  = ProjD(DD)
%           T1   = DD(DD>=0) ;
%           T2   =  DD(DD<0)  ;
%           DD(DD>=0) =  min(T1, 1e-4) ;
%           DD(DD<0)    = max(T2,-1e-4) ;
%           Dp = DD;
%     end


function L =  CalLaplacian(B)
   L = (diag( B * ones(size(B,2), 1) )) - B;
end


function CaldiffZ
%     diffZ(iter)= max(cellfun(@(X,Y) max( max(abs( X - Y ) ) ), Z, Zk )) ;
diffZ(iter) =  sum(cellfun(@(X,Y) sumsqr( X - Y ) , Z, Zk )) ;
end

function CaldiffC
%     diffC(iter) = max(max(abs(C - Ck)));    
      diffC(iter) = sumsqr(C - Ck) ;    
end
function CaldiffW
%     diffW(iter) = max( max( abs(W- Wold)));    
diffW(iter) = sumsqr(W - Wold) ;
end

function CaldiffD
%      diffD(iter)= max(cellfun(@(X,Y)  max(max(abs(X - Y ))), D, Dk )) ;
diffD(iter) =  sum(cellfun(@(X,Y) sumsqr( X - Y ) , D, Dk )) ;
end

function CalLoss
   l1 =0.5 * sum(cellfun(@(x,y)sumsqr(x-x*y), X, Z)) ;
   l2 = 0.5 * lambda1 * sum(cellfun(@(x,y) sumsqr(x-C.*y), Z, D));
   l3  = 0.5 * lambda2 * sumsqr(D) ;
   l4 = 0.5 * lambda3 * sum(sum(W.*CalLaplacian(C)));
   loss(iter) = l1 +l2 +l3+l4 ;
end


end