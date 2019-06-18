function res = CellFunc(handle,varargin)
	res = cellfun(handle,varargin{:}, 'UniformOutput', false) ;
end