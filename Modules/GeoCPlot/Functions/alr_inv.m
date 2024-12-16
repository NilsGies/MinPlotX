    function inv=alrinv(X,pos)
   % inverse of alr-transfo
    
    [m,n]=size(X);
    if nargin<2 pos=n+1; end;
    X_d=ones(m,1);
    X=exp(X);
    X_d=X_d+sum(X')';
    X_d=1./X_d;
    for i=1:n
        X(:,i)=X_d.*X(:,i);
    end;
    X=[X(:,1:(pos-1)),X_d,X(:,(pos:end))];
    inv=X;

	
