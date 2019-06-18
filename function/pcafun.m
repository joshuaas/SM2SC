function [U,Xpred] = pcafun( X,ratio )
%PCAFUN Summary of this function goes here
%   Detailed explanation goes here
if nargin ==1
    ratio = 0.95;
end
if ratio == int16(ratio)
    mod = 'dim' ;
else
    mod = 'rat'  ;
end

[n, d]  = size(X) ;
H =  speye(n) - (1/n) * ones(n);
X = H * X;
variance =  X'*X/(n - 1);
[U,sigma] = eigs(variance, d) ;
 if  strcmp(mod,  'rat')  
    U = U(:, fliplr( 1 : d )  );
    sigma = diag(sigma) ;
    sigma =sigma( end : -1 : 1) ;
    rat   = cumsum(sigma) /sum(sigma);
    k = sum(rat < ratio);
 else
     k = ratio ;
 end

    U = U(:, 1:k);

if nargout ==  2
    Xpred = (X *U);
end

end

