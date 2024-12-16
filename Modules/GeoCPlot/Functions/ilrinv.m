	function ilrinv=ilrinv(Y)
	%ilrinv makes the inverse of
	% ILR-TRANSFORMATION of (D-1)-dimensional
	% real data set Y 
	% returns a D-compositional data set

	[r c]=size(Y);	
	D=c+1;
	q=zeros(r,D);
	
	for k=1:r 
        q(k,1)=1;
        q(k,2)=exp(-sqrt(2)*Y(k,1));
        for i=3:D 
            q(k,i)=(exp(sqrt((i-2)*(i-1))*Y(k,i-2))/exp(sqrt((i-1)*i)*Y(k,i-1)))^(1/(i-1)); 
        end; 
    end;
	sumq=zeros(r,1);
	for k=1:r 
        for i=1:D sumq(k)=sumq(k)+prod(q(k,1:i)); 
        end; 
    end;
	ilrinv=zeros(r,D);
	for k=1:r
        ilrinv(k,1)=1/sumq(k);
        for i=2:D
        ilrinv(k,i)=q(k,i)*ilrinv(k,i-1);
        end;
    end;
	%ilrinv=raw2comp(ilrinv);
	

