%% Ilmenite structural formula
% last modified 01.08.2024
function [StrctFrm, apfu,options_definition]=ilmenite_Fe3unknown(data,headers,options)
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

cat=2.0; %cations per formula unit
Opfu=3.0; %oxygens per formula unit

%% Molecular weights

SiO2_mw=60.083;
TiO2_mw=79.865;
Al2O3_mw=101.961;
Cr2O3_mw=151.989;
Fe2O3_mw=159.6874;
Mn2O3_mw=157.873;
V2O3_mw=149.881;
Y2O3_mw=225.809;
NiO_mw=74.692;
ZnO_mw=81.381;
FeO_mw=71.8442;
MnO_mw=70.937;
MgO_mw=40.304;
CaO_mw=56.0774;
Na2O_mw=61.979;
K2O_mw=94.195;
BaO_mw=153.329;
F_mw=18.998;
Cl_mw=35.45;

%% Calculate cations units

MC=zeros(m,9); %create columns of zeros if optional data are not included

MC(:,1)=data(:,strcmp(headers,'TiO2'))./TiO2_mw; %for TiO2

%adds a column of zeros if SiO2 is not included in the calculation
if any(strcmp(headers,'SiO2'))
    MC(:,2)=data(:,strcmp(headers,'SiO2'))./SiO2_mw; %for SiO2
end

MC(:,3)=(data(:,strcmp(headers,'Al2O3'))./Al2O3_mw).*2; %for Al2O3
MC(:,4)=(data(:,strcmp(headers,'Cr2O3'))./Cr2O3_mw).*2; %for Cr2O3

%adds a column of zeros if V2O3 is not included in the calculation
if any(strcmp(headers,'V2O3'))
    MC(:,5)=(data(:,strcmp(headers,'V2O3'))./V2O3_mw).*2; %for V2O3
end

%calculates for FeO and/or Fe2O3

if any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) %if FeO and Fe2O3 are both included as inputs 

     %If FeO and Fe2O3 are both given, Fe2O3 is converted to FeO and
     %combined with FeO
     MC(:,6)=(data(:,strcmp(headers,'Fe2O3')).*((2*FeO_mw)./Fe2O3_mw)+data(:,strcmp(headers,'FeO')))./FeO_mw;
else
    if any(strcmp(headers,'FeO')) %if FeO is FeO total

        MC(:,6)=data(:,strcmp(headers,'FeO'))./FeO_mw;

    else %if Fe2O3 total is given, converted to FeO total

        %calculate moles of cations
        MC(:,6)=(data(:,strcmp(headers,'Fe2O3')).*((2*FeO_mw)./Fe2O3_mw))./FeO_mw;

    end
end

%calculates for MnO and/or Mn2O3

if any(strcmp(headers,'MnO')) && any(strcmp(headers,'Mn2O3')) %if MnO and Mn2O3 are both included as inputs 

     %If MnO and Mn2O3 are both given, Mn2O3 is converted to MnO and
     %combined with MnO
     MC(:,7)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw)+data(:,strcmp(headers,'MnO')))./MnO_mw;
else
    if any(strcmp(headers,'MnO')) %if MnO is MnO total

        MC(:,7)=data(:,strcmp(headers,'MnO'))./MnO_mw;

    else %if Mn2O3 total is given, converted to MnO total

        %calculate moles of cations
        MC(:,7)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw))./MnO_mw;

    end
end

%adds a column of zeros if MgO is not included in the calculation
if any(strcmp(headers,'MgO'))
    MC(:,8)=data(:,strcmp(headers,'MgO'))./MgO_mw; %for MgO
end

%adds a column of zeros if CaO is not included in the calculation
if any(strcmp(headers,'CaO'))
    MC(:,9)=data(:,strcmp(headers,'CaO'))./CaO_mw; %for CaO
end

MCnormfact=cat./sum(MC,2); %normalization factor

%% Calculate normalized cations units

MCnorm=MCnormfact.*MC; %creates a matrix of normalized cation

%% Calculate Oxygen Units

O2(:,1)=MCnorm(:,1).*2; %for TiO2
O2(:,2)=MCnorm(:,2).*2; %for TiO2
O2(:,3)=MCnorm(:,3).*(3/2); %for Al2O3
O2(:,4)=MCnorm(:,4).*(3/2); %for Cr2O3
O2(:,5)=MCnorm(:,5).*(3/2); %for V2O3
O2(:,6)=MCnorm(:,6); %for FeO
O2(:,7)=MCnorm(:,7); %for MnO
O2(:,8)=MCnorm(:,8); %for MgO
O2(:,9)=MCnorm(:,9); %for CaO

O2total=sum(O2,2); %O2 totals

%% Atoms PFU
apfu=zeros(m,12); 

apfu(:,1)=MCnorm(:,1); %for Ti
apfu(:,2)=MCnorm(:,2); %for Si
apfu(:,3)=MCnorm(:,3); %for Al
apfu(:,4)=MCnorm(:,4); %for Cr
apfu(:,5)=MCnorm(:,5); %for V

apfu(:,9)=MCnorm(:,7); %for Mn2+
apfu(:,10)=MCnorm(:,8); %for Mg
apfu(:,11)=MCnorm(:,9); %for Ca

%calculation of Fe3+ from stoichiometry and charge balance
%the following if statement firsts checks if totalO2 = 3
%if so, then there is no Fe3+
%if totalO2 < 3, then we assume that the deficiency is caused by the
%assumption Fetotal = Fe2+
%in the nested if statement, if FeTotal > 2*(3-totalO2) then the amount
%of Fe3+ = 2*(3-totalO2), if false then, all Fe is Fe3+

for c=1:m
    if (Opfu-O2total(c,1)) >= 0
        if MCnorm(c,6) > 2*(Opfu-O2total(c,1))
            apfu(c,6)=2*(Opfu-O2total(c,1)); 
        else
            apfu(c,6)=MCnorm(c,6);
        end
    else
        apfu(c,6)=0;
    end
end

apfu(:,8)=MCnorm(:,6)-apfu(:,6); %the apfu of Fe2+ equals totalFe-Fe3+

apfu(:,12)=sum(apfu,2); %calculations the total, which should be 3

% Oxygen deficiency 
apfu(:,13)=Opfu-O2total; %must be greater than zero

%% Endmembers 

XHem=apfu(:,6)./(apfu(:,6)+apfu(:,1));
XIlGkPy=1-XHem;
XIlm=XIlGkPy.*(apfu(:,8)./(apfu(:,8)+apfu(:,9)+apfu(:,10)));
XPph=XIlGkPy.*(apfu(:,9)./(apfu(:,8)+apfu(:,9)+apfu(:,10)));
XGk=XIlGkPy.*(apfu(:,10)./(apfu(:,8)+apfu(:,9)+apfu(:,10)));

all=[apfu XIlm XGk XPph XHem];

all(all<1e-6) = 0;
all(isnan(all))=0;
all(:,13)=Opfu-O2total; %must be greater than zero

StrctFrm=array2table(all,'VariableNames',{'Ti','Si','Al','Cr','V','Fe3','Mn3','Fe2','Mn2','Mg','Ca','Sum','O2_deficiency','Xilm','Xgk','Xpph','Xhem'});
apfu=array2table(all,'VariableNames',{'Ti','Si','Al','Cr','V','Fe3','Mn3','Fe2','Mn2','Mg','Ca','Sum','O2_deficiency','Xilm','Xgk','Xpph','Xhem'});

end 

