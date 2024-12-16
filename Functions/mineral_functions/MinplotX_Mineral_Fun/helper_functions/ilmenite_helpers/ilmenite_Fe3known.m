%% Ilmenite structural formula
% last modified 01.08.2024

function [StrctFrm, apfu,options_definition]=ilmenite_Fe3known(data,headers,options)
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

%% Checks for Fe3+/Fetotal ratio
%assigns Fe3/Fetotal if it is included in the analysis

if isfield(options,'Fe3_ratio') && isfield(options,'UseKnownFe3') && options.UseKnownFe3.Value==true
    Fe_rat=options.Fe3_ratio.Value;
    options.RecalcTotalFe.Value=true;
else
    Fe_rat(:,1)=zeros(m,1); %Fe3+/FeTotal ratio is 0
    options.RecalcTotalFe.Value=false;
end
%% Checks for Mn3+/Mntotal ratio

%assigns Mn3/Mntotal if it is included in the analysis
if isfield(options,'Mn3_ratio') && isfield(options,'UseKnownMn3') && options.UseKnownMn3.Value==true
    Mn_rat=options.Mn3_ratio.Value;
    options.RecalcTotalMn.Value=true;
else
    Mn_rat(:,1)=ones(m,1); %Mn3+/MnTotal ratio is 1
    options.RecalcTotalMn.Value=false;
end


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

MC=zeros(m,11); %create columns of zeros if optional data are not included

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

%calculates for FeO and Fe2O3

if any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) && isfield(options,'RecalcTotalFe') && options.RecalcTotalFe.Value==true

    %recalculate total FeO based on FeO+Fe2O3:
    FeOT=(data(:,strcmp(headers,'Fe2O3')).*((2*FeO_mw)./Fe2O3_mw)+data(:,strcmp(headers,'FeO')));

    %first calculate the weight percent of Fe2O3 and FeO from FeO total
    Fe3_wtper(:,1)=(FeOT.*Fe_rat)./((2*FeO_mw)./Fe2O3_mw); %calculates Fe2O3 from Fe3+ ratio and total FeO
    Fe2_wtper(:,1)=FeOT-(FeOT.*Fe_rat); %FeO = Total FeO - Fe3+ in Total FeO

    %calculate moles of cations
    MC(:,6)=(Fe3_wtper(:,1)./Fe2O3_mw).*2;
    MC(:,8)=Fe2_wtper(:,1)./FeO_mw;

elseif any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) %if FeO and Fe2O3 are both included as inputs

     MC(:,6)=(data(:,strcmp(headers,'Fe2O3'))./Fe2O3_mw).*2;
     MC(:,8)=data(:,strcmp(headers,'FeO'))./FeO_mw;
else
    if any(strcmp(headers,'FeO')) %if FeO is FeO total

        %first calculate the weight percent of Fe2O3 and FeO from FeO total
        Fe3_wtper(:,1)=(data(:,strcmp(headers,'FeO')).*Fe_rat)./((2*FeO_mw)./Fe2O3_mw); %calculates Fe2O3 from Fe3+ ratio and total FeO
        Fe2_wtper(:,1)=data(:,strcmp(headers,'FeO'))-(data(:,strcmp(headers,'FeO')).*Fe_rat); %FeO = Total FeO - Fe3+ in Total FeO

        %calculate moles of cations
        MC(:,6)=(Fe3_wtper(:,1)./Fe2O3_mw).*2;
        MC(:,8)=Fe2_wtper(:,1)./FeO_mw;

    else %if Fe2O3 is Fe2O3 total

        %first calculate the weight percent of Fe2O3 and FeO from Fe2O3 total
        Fe3_wtper(:,1)=data(:,strcmp(headers,'Fe2O3'))-(data(:,strcmp(headers,'Fe2O3')).*(1-Fe_rat)); %Fe2O3 = Total Fe2O3 - Fe2+ in Total Fe2O3 
        Fe2_wtper(:,1)=(data(:,strcmp(headers,'Fe2O3')).*(1-Fe_rat)).*((2*FeO_mw)./Fe2O3_mw); %calculates FeO from Fe3+ ratio and total Fe2O3

        %calculate moles of cations
        MC(:,6)=(Fe3_wtper(:,1)./Fe2O3_mw).*2;
        MC(:,8)=Fe2_wtper(:,1)./FeO_mw;

    end
end

%calculates for MnO and Mn2O3

if any(strcmp(headers,'MnO')) && any(strcmp(headers,'Mn2O3')) && isfield(options,'RecalcTotalMn') && options.RecalcTotalMn.Value==true

    %recalculate total MnO based on MnO+Mn2O3:
    MnOT=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw)+data(:,strcmp(headers,'MnO')));

    %first calculate the weight percent of Mn2O3 and MnO from MnO total
    Mn3_wtper(:,1)=(MnOT.*Mn_rat)./((2*MnO_mw)./Mn2O3_mw); %calculates Mn2O3 from Mn3+ ratio and total MnO
    Mn2_wtper(:,1)=MnOT-(MnOT.*Mn_rat); %MnO = Total MnO - Mn3+ in Total MnO

    %calculate moles of cations
    MC(:,7)=(Mn3_wtper(:,1)./Mn2O3_mw).*2;
    MC(:,9)=Mn2_wtper(:,1)./MnO_mw;

elseif any(strcmp(headers,'MnO')) && any(strcmp(headers,'Mn2O3')) %if MnO and Mn2O3 are both included as inputs

     MC(:,7)=(data(:,strcmp(headers,'Mn2O3'))./Mn2O3_mw).*2;
     MC(:,9)=data(:,strcmp(headers,'MnO'))./MnO_mw;
else
    if any(strcmp(headers,'MnO')) %if MnO is MnO total

        %first calculate the weight percent of Mn2O3 and MnO from MnO total
        Mn3_wtper(:,1)=(data(:,strcmp(headers,'MnO')).*Mn_rat)./((2*MnO_mw)./Mn2O3_mw); %calculates Mn2O3 from Mn3+ ratio and total MnO
        Mn2_wtper(:,1)=data(:,strcmp(headers,'MnO'))-(data(:,strcmp(headers,'MnO')).*Mn_rat); %MnO = Total MnO - Mn3+ in Total MnO

        %calculate moles of cations
        MC(:,7)=(Mn3_wtper(:,1)./Mn2O3_mw).*2;
        MC(:,9)=Mn2_wtper(:,1)./MnO_mw;

    elseif any(strcmp(headers,'Mn2O3')) %if Mn2O3 is Mn2O3 total

        %first calculate the weight percent of Fe2O3 and FeO from Fe2O3 total
        Mn3_wtper(:,1)=data(:,strcmp(headers,'Mn2O3'))-(data(:,strcmp(headers,'Mn2O3')).*(1-Mn_rat)); %Mn2O3 = Total Mn2O3 - Mn2+ in Total Mn2O3 
        Mn2_wtper(:,1)=(data(:,strcmp(headers,'Mn2O3')).*(1-Mn_rat)).*((2*MnO_mw)./Mn2O3_mw); %calculates MnO from Mn3+ ratio and total Mn2O3

        %calculate moles of cations
        MC(:,7)=(Mn3_wtper(:,1)./Mn2O3_mw).*2;
        MC(:,9)=Mn2_wtper(:,1)./MnO_mw;
    else
        %Mn2O3 and MnO are not included, Mn2 and Mn3 = 0
    end
end

%adds a column of zeros if MgO is not included in the calculation
if any(strcmp(headers,'MgO'))
    MC(:,10)=data(:,strcmp(headers,'MgO'))./MgO_mw; %for MgO
end

%adds a column of zeros if CaO is not included in the calculation
if any(strcmp(headers,'CaO'))
    MC(:,11)=data(:,strcmp(headers,'CaO'))./CaO_mw; %for CaO
end

%% Calculate Oxygen Units

O2(:,1)=MC(:,1).*2; %for TiO2
O2(:,2)=MC(:,2).*2; %for SiO2
O2(:,3)=MC(:,3).*(3/2); %for Al2O3
O2(:,4)=MC(:,4).*(3/2); %for Cr2O3
O2(:,5)=MC(:,5).*(3/2); %for V2O3
O2(:,6)=MC(:,6).*(3/2); %for Fe2O3
O2(:,7)=MC(:,7).*(3/2); %for Mn2O3
O2(:,8)=MC(:,8); %for FeO
O2(:,9)=MC(:,9); %for MnO
O2(:,10)=MC(:,10); %for MgO
O2(:,11)=MC(:,11); %for CaO

O2total=sum(O2,2); %O2 totals
MCnormfact=Opfu./O2total; %normalization factor

apfu=MCnormfact.*MC; %creates a matrix of normalized cations
apfu(:,12)=sum(apfu,2); %calculations the total, which should be close to 2
apfu(:,13)=zeros(m,1); %zero in the case of no Fe estimation

%% Endmembers 

XHem=apfu(:,6)./(apfu(:,6)+apfu(:,1));
XIlGkPy=1-XHem;
XIlm=XIlGkPy.*(apfu(:,8)./(apfu(:,8)+apfu(:,9)+apfu(:,10)));
XPph=XIlGkPy.*(apfu(:,9)./(apfu(:,8)+apfu(:,9)+apfu(:,10)));
XGk=XIlGkPy.*(apfu(:,10)./(apfu(:,8)+apfu(:,9)+apfu(:,10)));

all=[apfu XIlm XGk XPph XHem];

all(all<1e-6) = 0;
all(isnan(all))=0;

StrctFrm=array2table(all,'VariableNames',{'Ti','Si','Al','Cr','V','Fe3','Mn3','Fe2','Mn2','Mg','Ca','Sum','O2_deficiency','Xilm','Xgk','Xpph','Xhem'});
apfu=array2table(all,'VariableNames',{'Ti','Si','Al','Cr','V','Fe3','Mn3','Fe2','Mn2','Mg','Ca','Sum','O2_deficiency','Xilm','Xgk','Xpph','Xhem'});

end 

