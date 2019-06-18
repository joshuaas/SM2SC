function Zmean =  MergeAdjacentMatrix(Z)
 Zmean =zeros(size(Z{1}));
 for i = 1:numel(Z)
 	Zmean  =  Zmean + (abs(Z{i}) + abs(Z{i}'))/2;
 end
 Zmean = Zmean/ numel(Z);
end