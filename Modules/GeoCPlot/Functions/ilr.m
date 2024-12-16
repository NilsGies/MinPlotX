	function ilr=ilr(X)
	%ILR makes ILR-TRANSFORMATION 
	%of a D-compositional data set X 
	% returns a (D-1)-matrix of coeficients

	[fil,col]=size(X);
	
	l=log(X);
	s=l(:,1);
	ilr=[];
	for k=2:col;
	 s=[s,s(:,k-1)+l(:,k)];
	end;
	for k=1:(col-1)
	ilr=[ilr,(1/sqrt(k*(k+1)))*(s(:,k)-k*l(:,k+1))];
	end;

