%%feldspar structural formula
% last modified 03.07.2024

function [StrctFrm, apfu,options_definition]=feldspar_MinPlotX(data,headers,options)
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

Opfu=8.0; %oxygens per formula unit

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
Na2O_mw=61.979;
K2O_mw=94.195;
BaO_mw=153.329;
F_mw=18.998;
Cl_mw=35.45;

%% moles of cations
MC=zeros(m,10);

MC(:,1)=data(:,strcmp(headers,'SiO2'))./SiO2_mw; %for SiO2
MC(:,2)=(data(:,strcmp(headers,'Al2O3'))./Al2O3_mw).*2; %for Al2O3

%calculates for FeO and/or Fe2O3

if any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) %if FeO and Fe2O3 are both included as inputs

    %If FeO and Fe2O3 are both given, Fe2O3 is converted to FeO and
    %combined with FeO
    MC(:,3)=(data(:,strcmp(headers,'Fe2O3')).*((2.*FeO_mw)./Fe2O3_mw)+data(:,strcmp(headers,'FeO')))./FeO_mw;
else
    if any(strcmp(headers,'FeO')) %if FeO is FeO total

        MC(:,3)=data(:,strcmp(headers,'FeO'))./FeO_mw;

    elseif any(strcmp(headers,'Fe2O3')) %if Fe2O3 total is given, converted to FeO total

        %calculate moles of cations
        MC(:,3)=(data(:,strcmp(headers,'Fe2O3')).*((2.*FeO_mw)./Fe2O3_mw))./FeO_mw;

    else
        % if FeO total or Fe2O3 total is not included, Fe = 0
    end
end


%calculates for MnO and/or Mn2O3

if any(strcmp(headers,'MnO')) && any(strcmp(headers,'Mn2O3')) %if MnO and Mn2O3 are both included as inputs

    %If MnO and Mn2O3 are both given, Mn2O3 is converted to MnO and
    %combined with MnO
    MC(:,4)=(data(:,strcmp(headers,'Mn2O3')).*((2.*MnO_mw)./Mn2O3_mw)+data(:,strcmp(headers,'MnO')))./MnO_mw;
else
    if any(strcmp(headers,'MnO')) %if MnO is MnO total

        MC(:,4)=data(:,strcmp(headers,'MnO'))./MnO_mw;

    elseif any(strcmp(headers,'Mn2O3')) %if Mn2O3 total is given, converted to MnO total

        %calculate moles of cations
        MC(:,4)=(data(:,strcmp(headers,'Mn2O3')).*((2.*MnO_mw)./Mn2O3_mw))./MnO_mw;

    else
        % if MnO total or Mn2O3 total is not included, Mn = 0
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

%% Oxygen Units
O2(:,1)=MC(:,1).*2; %for SiO2
O2(:,2)=MC(:,2).*(3/2); %for Al2O3
O2(:,3)=MC(:,3); %for FeO
O2(:,4)=MC(:,4); %for MnO
O2(:,5)=MC(:,5); %for MgO
O2(:,6)=MC(:,6); %for CaO
O2(:,7)=MC(:,7)*(1/2); %for Na2O
O2(:,8)=MC(:,8)*(1/2); %for K2O
O2(:,9)=MC(:,9); %for BaO
O2(:,10)=MC(:,10); %for SrO

O2total=sum(O2,2); %O2 totals

MCnormfact=Opfu./sum(O2,2); %normalization factor

%% Structural Formula
StrctFrm=MCnormfact.*MC; %creates a matrix of normalized cations
StrctFrm(:,11)=sum(StrctFrm,2); %calculations the total, which should be close to 5

%% endmember calculation

Endmembers(:,1)=StrctFrm(:,6)./(StrctFrm(:,6)+StrctFrm(:,7)+StrctFrm(:,8)); %anorthite 
Endmembers(:,2)=StrctFrm(:,7)./(StrctFrm(:,6)+StrctFrm(:,7)+StrctFrm(:,8)); %albite
Endmembers(:,3)=StrctFrm(:,8)./(StrctFrm(:,6)+StrctFrm(:,7)+StrctFrm(:,8)); %orthoclase

Endmembers(Endmembers<1e-3) = 0; %limit on endmember noise (cannot be less than a fraction of a percent)

%% Outputs

all=[StrctFrm Endmembers];

%limit on significant digits (eliminates rounding noise)
all(all<1e-6) = 0;
all(isnan(all))=0;

StrctFrm=array2table(all,'VariableNames',{'Si_T','Al_T','Fe2_M1','Mn2_M1','Mg_M1','Ca_M1','Na_M1','K_M1','Ba_M1','Sr_M1','Cation_Sum','Xan','Xab','Xor'}); 
apfu=array2table(all,'VariableNames',{'Si','Al','Fe2','Mn2','Mg','Ca','Na','K','Ba','Sr','Cation_Sum','Xan','Xab','Xor'}); 
end
