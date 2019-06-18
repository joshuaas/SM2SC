clear all
acc = zeros(30,1) ;
nmi = zeros(30,1) ;
f = zeros(30,1) ;
addpath(genpath('./ClusteringMeasure'))
addpath(genpath('./function'))
% addpath(genpath('./kmeans'))
addpath(genpath('./CodefromSSC'))
addpath('./data_test')
%% bcc/msra
rng(1111131789)
load('best_msra.mat')
disp('testing msra .....')
for i = 1:30
    [acc(i), nmi(i), f(i)] = evalBDRMV(X, gt,bestparam, 1, 'B') ;
 %   disp(i)
end
disp('algorithm finished .....')
fprintf('MSRA acc: %.4f, nmi: %.4f, fsore: %.4f\n', mean(acc), mean(nmi), mean(f)) ;

load('best_bcc.mat')
rng(1111131789)
disp('testing bbc......')
for i = 1:30
    [acc(i), nmi(i), f(i)] = evalBDRMV(X, gt,bestparam, 1, 'B') ;
 %   disp(i)
end
disp('algorithm finished .....')
fprintf('BBC Sport acc: %.4f, nmi: %.4f, fsore: %.4f\n', mean(acc), mean(nmi), mean(f)) ;



%%orl

rng(46828391)
disp('testing orl......')
load('best_orl.mat')

for i = 1:30
    [acc(i), nmi(i), f(i)] = evalBDRMV(X, gt,bestparam, 1, 'B') ;
  %  disp(i)
end
disp('algorithm finished .....')
fprintf('ORL acc: %.4f, nmi: %.4f, fsore: %.4f\n', mean(acc), mean(nmi), mean(f)) ;
