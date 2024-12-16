function clr=clr(X)
%CLR calcula la TRANSFORMACIO CLR del conjunt X de dades

[fil,col]=size(X);
l=log(X);
clr=zeros(size(X));
for f=1:fil
    m=mean(l(f,:));
    clr(f,:)=l(f,:)-m;
end

