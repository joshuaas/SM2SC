function [Acc, NMI, Fscore, grps,CKSym,Ccell, embed,loss,diff] = evalBDRMV(X, gt,bestparam, rho,opt,exparam,rseed)
	nCluster  = max(unique(gt)) ;      
    
    if nargin >5
            [~,CDmean,Z, Zmean,C,loss,diff] = M2VBDR(X,gt,bestparam,exparam) ;
    else
            [~,CDmean,Z, Zmean,C,loss,diff] = M2VBDR(X,gt,bestparam) ;
    end
% 	[B,Z,loss] = BDR_solver(X,nCluster,lambda,gamma) ;
         if strcmp(opt, 'Z')
	 	CKSym = BuildAdjacency(thrC(Zmean,rho));
	 elseif strcmp(opt, 'B')
	 	 	CKSym = BuildAdjacency(thrC(CDmean,rho));
	 elseif strcmp(opt, 'C')
	 		CKSym = BuildAdjacency(thrC(C,rho));
         end
     Ccell =  cellfun(@BuildAdjacency, Z, 'UniformOutput',false) ;
     
         
%         CKSym = abs(Z'*Z./(nX'*nX));
        [grps, embed] = SpectralClustering(CKSym,nCluster);
%      BinG = FormBin(grps);
	[~ ,NMI, ~] = compute_nmi(gt,grps);
         Acc = Accuracy(grps,double(gt));
         Fscore = compute_f(gt,grps);
end