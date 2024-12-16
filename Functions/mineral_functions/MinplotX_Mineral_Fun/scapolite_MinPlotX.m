%%scapolite
%last modified 09.07.2024

function [StrctFrm, apfu,options_definition]=scapolite_MinPlotX(data,headers,options)
%% define empty output variables
StrctFrm=[];
apfu=[];
%% options definition
options_definition=struct();

%% Check if the data or header is empty 
if not(exist('headers','var')) || not(exist('data','var')) || isempty(data) || isempty(headers)
    return % exit the function and return options_definition
end

%% If no options 
if not(exist('options','var')) || isempty(options)
options=options_definition; % use default values generated in the options definition block
end
%% start calculation
[m,~]=size(data); %finds the x and y size of the input data matrix


%% Molecular weights

SiO2_mw=60.083;
TiO2_mw=79.865;
Al2O3_mw=101.961;
Cr2O3_mw=151.989;
Fe2O3_mw=159.6874;
Y2O3_mw=225.809;
NiO_mw=74.692;
ZnO_mw=81.381;
FeO_mw=71.8442;
MnO_mw=70.937;
Mn2O3_mw=157.873;
MgO_mw=40.304;
CaO_mw=56.0774;
SrO_mw=103.619;
SO3_mw=80.0594;
Na2O_mw=61.979;
K2O_mw=94.195;
BaO_mw=153.329;
F_mw=18.998;
Cl_mw=35.45;
H2O_mw=18.015;
CO2_mw=44.009;

%% Moles of cations

MC=zeros(m,14);%create columns of zeros if optional data are not included

MC(:,1)=data(:,strcmp(headers,'SiO2'))./SiO2_mw; % for SiO2
MC(:,2)=(data(:,strcmp(headers,'Al2O3'))./Al2O3_mw).*2; % for Al2O3

%calculates for FeO and/or Fe2O3
if any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) %if FeO and Fe2O3 are both included as inputs 

     %If FeO and Fe2O3 are both given, Fe2O3 is converted to FeO and
     %combined with FeO
     MC(:,3)=(data(:,strcmp(headers,'Fe2O3')).*((2*FeO_mw)./Fe2O3_mw)+data(:,strcmp(headers,'FeO')))./FeO_mw;
else
    if any(strcmp(headers,'FeO')) %if FeO is FeO total

        MC(:,3)=data(:,strcmp(headers,'FeO'))./FeO_mw;

    elseif any(strcmp(headers,'Fe2O3')) %if Fe2O3 total is given, converted to FeO total

        %calculate moles of cations
        MC(:,3)=(data(:,strcmp(headers,'Fe2O3')).*((2*FeO_mw)./Fe2O3_mw))./FeO_mw;

    else
        %Fe2O3 and FeO are not included, Fe2 = 0
    end
end

%calculates for MnO and/or Mn2O3
if any(strcmp(headers,'MnO')) && any(strcmp(headers,'Mn2O3')) %if MnO and Mn2O3 are both included as inputs 

     %If MnO and Mn2O3 are both given, Mn2O3 is converted to MnO and
     %combined with MnO
     MC(:,4)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw)+data(:,strcmp(headers,'MnO')))./MnO_mw;
else
    if any(strcmp(headers,'MnO')) %if MnO is MnO total

        MC(:,4)=data(:,strcmp(headers,'MnO'))./MnO_mw;

    elseif any(strcmp(headers,'Mn2O3'))%if Mn2O3 total is given, converted to MnO total

        %calculate moles of cations
        MC(:,4)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw))./MnO_mw;
    else
        %Mn2O3 and MnO are not included, Mn2 = 0
    end
end

%calculates for MgO if it is included in the analysis
if any(strcmp(headers,'MgO'))
    MC(:,5)=data(:,strcmp(headers,'MgO'))./MgO_mw; %for MgO
end

MC(:,6)=data(:,strcmp(headers,'CaO'))./CaO_mw; %for CaO
MC(:,7)=(data(:,strcmp(headers,'Na2O'))./Na2O_mw).*2; %for Na2O
MC(:,8)=(data(:,strcmp(headers,'K2O'))./K2O_mw).*2; %for K2O

%calculates for BaO if it is included in the analysis
if any(strcmp(headers,'BaO'))
    MC(:,9)=data(:,strcmp(headers,'BaO'))./BaO_mw; %for BaO
end

%calculates for SrO if it is included in the analysis
if any(strcmp(headers,'SrO'))
    MC(:,10)=data(:,strcmp(headers,'SrO'))./SrO_mw; %for SrO
end

%calculates for F if it is included in the analysis
if any(strcmp(headers,'F'))
    MC(:,11)=data(:,strcmp(headers,'F'))./F_mw; %for F
end

%calculates for Cl if it is included in the analysis
if any(strcmp(headers,'Cl'))
    MC(:,12)=data(:,strcmp(headers,'Cl'))./Cl_mw; %for Cl
end

%calculates for SO3 if it is included in the analysis
if any(strcmp(headers,'SO3'))
    MC(:,13)=data(:,strcmp(headers,'SO3'))./SO3_mw; %for SO3
end

%calculates for CO2 if it is included in the analysis
if any(strcmp(headers,'CO2'))
    MC(:,14)=(data(:,strcmp(headers,'CO2'))./CO2_mw); %for CO2
end

%% Moles of O2 units
O2(:,1)=MC(:,1).*(2); %SiO2
O2(:,2)=MC(:,2).*(3/2); %Al2O3
O2(:,3)=MC(:,3); %FeO
O2(:,4)=MC(:,4); %MnO
O2(:,5)=MC(:,5); %MgO
O2(:,6)=MC(:,6); %CaO
O2(:,7)=MC(:,7).*(1/2); %Na2O
O2(:,8)=MC(:,8).*(1/2); %K2O
O2(:,9)=MC(:,9); %BaO
O2(:,10)=MC(:,10); %SrO
O2(:,11)=MC(:,11); %F
O2(:,12)=MC(:,12); %Cl
O2(:,13)=MC(:,13).*4; %SO4
O2(:,14)=MC(:,14).*3; %CO3

if any(strcmp(headers,'CO2')) %if CO2 is known

    O2_Sum=sum(O2(:,1:14),2)-0.5.*(O2(:,11)+O2(:,12)); %sum of O2, excluding F, S1-, S2-, and Cl

    O2_N=(27)./O2_Sum; %normalization factor

    %normalized moles of anions
    N_Ox=O2.*O2_N;

else %method of Ketcham (2015) if H2O is not known

    O2_Sum=sum(O2(:,1:14),2)-O2(:,12); %sum of O2, excluding F, S1-, S2-, and Cl

    O2_N=(24)./O2_Sum; %normalization factor

    %normalized moles of anions
    N_Ox=O2.*O2_N;
end

%% Moles of cations (Oxygen Normalization)

apfu_O(:,1)=N_Ox(:,1)./2; %Si
apfu_O(:,2)=N_Ox(:,2).*(2/3); %Al
apfu_O(:,3)=N_Ox(:,3); %Fe2
apfu_O(:,4)=N_Ox(:,4); %Mn2
apfu_O(:,5)=N_Ox(:,5); %Mg
apfu_O(:,6)=N_Ox(:,6); %Ca
apfu_O(:,7)=N_Ox(:,7).*2; %Na
apfu_O(:,8)=N_Ox(:,8).*2; %K
apfu_O(:,9)=N_Ox(:,9); %Ba
apfu_O(:,10)=N_Ox(:,10); %Sr
apfu_O(:,11)=N_Ox(:,11); %F
apfu_O(:,12)=N_Ox(:,12); %Cl
apfu_O(:,13)=N_Ox(:,13).*(1/4); %SO3
apfu_O(:,14)=N_Ox(:,14).*(1/3); %CO2

%% normalization to 12 Si + Al

%Normmalization to Si + Al 
Nfact=12./(apfu_O(:,1)+apfu_O(:,2)); %normalization

apfu_N=apfu_O.*Nfact; %normalization of the cations


%if CO2 is known, the known value is taken
if not(any(strcmp(headers,'CO2')))
    apfu_N(:,14)=1-apfu_N(:,11)-apfu_N(:,12)-apfu_N(:,13); %CO2
end

apfu_N(:,15)=sum(apfu_N(:,1:10),2); %calculations the cation total
apfu=apfu_N;


%% endmember calculations
XEqAn=(apfu(:,2)-3)./3; %calculation of equivalent anorthite
XMe=(apfu(:,6)+apfu(:,5)+apfu(:,10)+apfu(:,9)+apfu(:,4)+apfu(:,3))/4; %sum of divalent cations divided by 4

%limit on significant digits (eliminates rounding noise)
apfu(apfu<1e-6) = 0;
XEqAn(XEqAn<1e-3) = 0; %limit on endmember noise (cannot be less than a fraction of a percent)
XMe(XMe<1e-3) = 0; %limit on endmember noise (cannot be less than a fraction of a percent)


all=[apfu XEqAn XMe];

StrctFrm=array2table(all,'VariableNames',{'Si_T','Al_T','Fe2_M','Mn2_M','Mg_M','Ca_M','Na_M','K_M','Ba_M','Sr_M','F','Cl','S6','C','Cation_Sum','Xeqan','Xme'});
apfu=array2table(apfu,'VariableNames',{'Si','Al','Fe2','Mn2','Mg','Ca','Na','K','Ba','Sr','F','Cl','S','C','Cation_Sum'});
end



