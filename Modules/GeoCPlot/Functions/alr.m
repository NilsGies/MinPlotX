function alr=alr(X,d)
% alr-transfo.
m=size(X,1);
if nargin<2
    d=size(X,2);
end
l=log(X);
for i=1:m
    alr(i,:)=l(i,:)-l(i,d);
end
alr(:,d)=[];

