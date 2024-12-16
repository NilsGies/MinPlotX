function plot_profiles_NG(X,Y,value,Id,mask,wn,in,Image,ImageLimits,xlimit)
% X=X(mask)-min(X);
% Y=Y(mask)-min(Y);
% X=linspace(min(X),max(X),numel(X))';
% 
% X=X';
% Y=Y';
%Plot of the spot positions
figure
hold on
plot(X,Y,'o');
text(X,Y,num2str(Id(mask)))
axis image
xlabel('X position')
ylabel('Y position')
title('Spot positions')

%mask=[]

[p,~]=polyfit(X(mask),flipud(Y(mask)),1);
a=p(1);
b=p(2);


%%
figure;hold on;
Y2=linspace(min(Y),max(Y),numel(Y))';
plot(X(mask),a*Y(mask)+b,'-r');
plot(X,Y,'kx')
%%
hold on
plot(X(mask),a*Y(mask)+b,'-r');

for i=1:numel(mask)
    c=Y(i)-a*(X(i));
    xH2(i)=(-b+Y(i)+1/a*X(i))/(a+1/a);
    yH2(i)=a*xH2(i)+b;
    plot(xH2(i),yH2(i),'.r')
    Projection=[X(i),xH2(i)];
    plot(Projection,-1/a*Projection+Y(i)+1/a*X(i),'--k')
end

%saveas(gcf,[String,'_Projections.pdf'])

xH2=xH2(xH2>0);

%Profile
Xmin=min(xH2);
Xprofile=xH2-Xmin+1; %scale from 1 µm

[~,Idx]=sort(Xprofile);

xx=Xprofile;
yy=value(mask);

figure
hold on

plot(xx(Idx),yy(Idx),'Color','k','Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k');

xx1=xx(Idx);
yy1=yy(Idx);
name1=mask(Idx)-1;
for j=1:length(xx1)
    text(xx1(j),yy1(j),num2str(name1(j)))
end
xlabel('Profile µm')
ylabel('Water content (ppm H2O)')
%title(String)
%saveas(gcf,[String,'_Profile.pdf'])

p=polyfit(xx1,yy1',2);
hold on
plot([xx1(1):0.1:xx1(end)],polyval(p,[xx1(1):0.1:xx1(end)]),'-')

%% Plot3D spectra

ColorCode1=[0 0.4470 0.7410];

figure
nexttile
hold on

for i=1:numel(mask)
  %  hold on
    plot3(wn(:),Xprofile(i).*ones(length(wn),1),in(:,i),'-','Color',ColorCode1);
end

%title(String);
xlabel('Wavenumbers [cm^{-1}]');
ylabel('Profile [µm]');
zlabel('Absorbance [cm^{-1}]');
set(gca,'XDir','reverse');

% zlim([-0.1 1])
%xlim(sort(xlimit))
view(45,45)
grid on

%saveas(gcf,[String,'_Stack.pdf'])

% nexttile
% hold on
% 
% for i=1:numel(mask)
% 
%     in(:,i)=linCorrect(in(:,i));
%   %  hold on
%     plot3(wn(:),Xprofile(i).*ones(length(wn),1),in(:,i),'-','Color',ColorCode1);
% end
% 
% %title(String);
% xlabel('Wavenumbers [cm^{-1}]');
% ylabel('Profile [µm]');
% zlabel('Absorbance [cm^{-1}]');
% set(gca,'XDir','reverse');
% 
% % zlim([-0.1 1])
% view(45,45)
% grid on

nexttile
hold on

%
range=wn<max(xlimit)&wn>min(xlimit);
in=in(range,:);
wn=wn(range);

for i=1:numel(mask)

    in(:,i)=linCorrect(in(:,i));
  %  hold on
    plot3(wn(:),Xprofile(i).*ones(length(wn),1),in(:,i),'-','Color',ColorCode1);
end

%title(String);
xlabel('Wavenumbers [cm^{-1}]');
ylabel('Profile [µm]');
zlabel('Absorbance [cm^{-1}]');
set(gca,'XDir','reverse');

% zlim([-0.1 1])
xlim(sort(xlimit))
view(45,45)
grid on

%saveas(gcf,[String,'_Stack.pdf'])


%% Plot position on image
figure
hold on
axis on
if not (isempty(Image)) && not (isempty(ImageLimits)) 
image(flipud(Image),'XData', ImageLimits(1,:), 'YData',ImageLimits(2,:));

end

set(gca,'YDir','normal')
axis image


%Plot the elements
hold on
plot(X(mask),Y(mask),'Marker','o','MarkerEdgeColor','w','MarkerFaceColor',[0.4660 0.6740 0.1880],'MarkerSize',7,'LineStyle','none');

plot(X(mask),a*X(mask)+b,'-r','LineWidth',1);
for i=1:numel(mask)
    c=Y(i)-a*(X(i));

    xH(i)=(-b+Y(i)+1/a*X(i))/(a+1/a);
    yH(i)=a*xH(i)+b;
    plot(xH(i),yH(i),'.r','MarkerSize',10)
    Projection=[X(i),xH(i)];
    plot(Projection,-1/a*Projection+Y(i)+1/a*X(i),'--k','LineWidth',0.75)
end

%title(String)
%saveas(gcf,[String,'_Img.pdf'])

end


