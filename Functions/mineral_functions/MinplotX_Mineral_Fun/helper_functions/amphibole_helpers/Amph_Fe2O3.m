function apfu=Amph_Fe2O3(D,W,O2_Ti)

% All ferrous iron, 23 O2 normalization

[m,n]=size(D); %finds the x and y size of the input data matrix

%moles of cations
D_Fe2O3=D;
D_Fe2O3(:,5)=D(:,5)+D(:,6)*(159.6882./(2*71.8444)); %convert FeO to Fe2O3
D_Fe2O3(:,6)=D(:,6)-D(:,6); %remove FeO

%moles of cations
MC=D_Fe2O3./(W); 

%Moles of O2
O2(:,1)=MC(:,1)*2; %SiO2
O2(:,2)=MC(:,2)*2; %TiO2
O2(:,3)=MC(:,3)*3; %Al2O3
O2(:,4)=MC(:,4)*3; %Cr2O3
O2(:,5)=MC(:,5)*3; %Fe2O3
O2(:,6)=MC(:,6); %FeO
O2(:,7)=MC(:,7); %MnO
O2(:,8)=MC(:,8); %MgO
O2(:,9)=MC(:,9); %CaO
O2(:,10)=MC(:,10); %Na2O
O2(:,11)=MC(:,11); %K2O

O2total=sum(O2,2); %O2 totals

O2_N=O2_Ti./O2total; %O2 normalization factor

%moles of cations
N_Ox=O2.*O2_N;

%atoms per formula unit
apfu(:,1)=N_Ox(:,1)./2; %SiO2, Criteria 1-1 (see below)
apfu(:,2)=N_Ox(:,2)./2; %TiO2
apfu(:,3)=N_Ox(:,3).*(2/3); %Al2O3
apfu(:,4)=N_Ox(:,4).*(2/3); %Cr2O3
apfu(:,5)=N_Ox(:,5).*(2/3); %Fe2O3
apfu(:,6)=N_Ox(:,6); %FeO
apfu(:,7)=N_Ox(:,7); %MnO
apfu(:,8)=N_Ox(:,8); %MgO
apfu(:,9)=N_Ox(:,9); %CaO
apfu(:,10)=N_Ox(:,10).*2; %Na2O
apfu(:,11)=N_Ox(:,11).*2; %K2O 

%Fe3+ estimation criteria 
apfu(:,12)=8./apfu(:,1); %Criteria 1-1: Si cannot be more than 8 cations

apfu(:,13)=16./sum(apfu(:,1:1:11),2); %Criteria 1-2: Should be <= 16 apfu

apfu(:,14)=15./sum(apfu(:,1:1:9),2); %Criteria 1-3: Should be <=15 apfu

apfu(:,15)=8./(apfu(:,1)+apfu(:,3)); % Criteria 2-1: Si + Al = 8

apfu(:,16)=15./sum(apfu(:,1:1:10),2); %Criteria 2-2: Should be equal to 15 apfu

apfu(:,17)=13./sum(apfu(:,1:1:8),2); %Criteria 2-3: Should be equal to 13 apfu

apfu(:,18)=36./(46-sum(apfu(:,1:4),2)); %Criteria 2-4: tetrahedral sites completely filled by 3+ and 4+ cation

apfu(:,19)=O2_N(:,1); 

end

