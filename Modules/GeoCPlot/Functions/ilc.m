function ilc=ilc(Z);

%Inverse finction of clr-transfo

[fil,col]=size(Z);
X=exp(Z);
S=sum(X,2);
for f=1:fil;
    ilc(f,:)=X(f,:)/S(f);
end;
