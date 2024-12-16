%%chloritoid Structural Formula
% last modified 03.07.2024

function [StrctFrm, apfu]=chloritoid_Fe3known(data,headers,options)
%% define empty output variables
StrctFrm=[];
apfu=[];

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

cat=8.0; %cations per formula unit
Opfu=12.0; %oxygens per formula unit

%% Checks for Fe3+/Fetotal ratio

%assigns Fe3/Fetotal if it is included in the analysis
if isfield(options,'Fe3_ratio') && isfield(options,'UseKnownFe3') && options.UseKnownFe3.Value==true
    Fe_rat=options.Fe3_ratio.Value;
    options.RecalcTotalFe.Value=true;
else
    Fe_rat(:,1)=zeros(m,1); %Fe3+/FeTotal ratio is 0
    options.RecalcTotalFe.Value=false;
end


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
H2O_mw=18.015;

%% Calculate moles of cations

MC=zeros(m,10); %create columns of zeros if optional data are not included

MC(:,1)=data(:,strcmp(headers,'SiO2'))./SiO2_mw; %for SiO2

%calculates for TiO2 if it is included in the analysis
if any(strcmp(headers,'TiO2'))
    MC(:,2)=data(:,strcmp(headers,'TiO2'))./TiO2_mw; %for TiO2
end

MC(:,3)=(data(:,strcmp(headers,'Al2O3'))./Al2O3_mw).*2; %for Al2O3

%calculates for FeO and Fe2O3

if any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) && isfield(options,'RecalcTotalFe') && options.RecalcTotalFe.Value==true

    %recalculate total FeO based on FeO+Fe2O3:
    FeOT=(data(:,strcmp(headers,'Fe2O3')).*((2*FeO_mw)./Fe2O3_mw)+data(:,strcmp(headers,'FeO')));

    %first calculate the weight percent of Fe2O3 and FeO from FeO total
    Fe3_wtper(:,1)=(FeOT.*Fe_rat)./((2*FeO_mw)./Fe2O3_mw); %calculates Fe2O3 from Fe3+ ratio and total FeO
    Fe2_wtper(:,1)=FeOT-(FeOT.*Fe_rat); %FeO = Total FeO - Fe3+ in Total FeO

    %calculate moles of cations
    MC(:,4)=(Fe3_wtper(:,1)./Fe2O3_mw).*2;
    MC(:,5)=Fe2_wtper(:,1)./FeO_mw;

elseif any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) %if FeO and Fe2O3 are both included as inputs

     MC(:,4)=(data(:,strcmp(headers,'Fe2O3'))./Fe2O3_mw).*2;
     MC(:,5)=data(:,strcmp(headers,'FeO'))./FeO_mw;
else
    if any(strcmp(headers,'FeO')) %if FeO is FeO total

        %first calculate the weight percent of Fe2O3 and FeO from FeO total
        Fe3_wtper(:,1)=(data(:,strcmp(headers,'FeO')).*Fe_rat)./((2*FeO_mw)./Fe2O3_mw); %calculates Fe2O3 from Fe3+ ratio and total FeO
        Fe2_wtper(:,1)=data(:,strcmp(headers,'FeO'))-(data(:,strcmp(headers,'FeO')).*Fe_rat); %FeO = Total FeO - Fe3+ in Total FeO

        %calculate moles of cations
        MC(:,4)=(Fe3_wtper(:,1)./Fe2O3_mw).*2;
        MC(:,5)=Fe2_wtper(:,1)./FeO_mw;

    else %if Fe2O3 is Fe2O3 total

        %first calculate the weight percent of Fe2O3 and FeO from Fe2O3 total
        Fe3_wtper(:,1)=data(:,strcmp(headers,'Fe2O3'))-(data(:,strcmp(headers,'Fe2O3')).*(1-Fe_rat)); %Fe2O3 = Total Fe2O3 - Fe2+ in Total Fe2O3 
        Fe2_wtper(:,1)=(data(:,strcmp(headers,'Fe2O3')).*(1-Fe_rat)).*((2*FeO_mw)./Fe2O3_mw); %calculates FeO from Fe3+ ratio and total Fe2O3

        %calculate moles of cations
        MC(:,4)=(Fe3_wtper(:,1)./Fe2O3_mw).*2;
        MC(:,5)=Fe2_wtper(:,1)./FeO_mw;

    end
end

%calculates for MnO and/or Mn2O3
if any(strcmp(headers,'MnO')) && any(strcmp(headers,'Mn2O3')) %if MnO and Mn2O3 are both included as inputs 

     %If MnO and Mn2O3 are both given, Mn2O3 is converted to MnO and
     %combined with MnO
     MC(:,6)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw)+data(:,strcmp(headers,'MnO')))./MnO_mw;
else
    if any(strcmp(headers,'MnO')) %if MnO is MnO total

        MC(:,6)=data(:,strcmp(headers,'MnO'))./MnO_mw;

    elseif any(strcmp(headers,'Mn2O3'))%if Mn2O3 total is given, converted to MnO total

        %calculate moles of cations
        MC(:,6)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw))./MnO_mw;
    else
        %Mn2O3 and MnO are not included, Mn2 = 0
    end
end

MC(:,7)=data(:,strcmp(headers,'MgO'))./MgO_mw; %for MgO

%calculates for CaO if it is included in the analysis
if any(strcmp(headers,'CaO'))
    MC(:,8)=data(:,strcmp(headers,'CaO'))./CaO_mw; %for CaO
end

%calculates for Na2O if it is included in the analysis
if any(strcmp(headers,'Na2O'))
    MC(:,9)=(data(:,strcmp(headers,'Na2O'))./Na2O_mw).*2; %for Na2O
end

%calculates for H2O if it is included in the analysis
if any(strcmp(headers,'H2O'))
    MC(:,10)=(data(:,strcmp(headers,'H2O'))./H2O_mw).*2; %for H2O
end

%% Calculate Oxygen Units

O2(:,1)=MC(:,1).*2; %for SiO2
O2(:,2)=MC(:,2).*2; %for TiO2
O2(:,3)=MC(:,3).*(3/2); %for Al2O3
O2(:,4)=MC(:,4).*(3/2); %for Fe2O3
O2(:,5)=MC(:,5); %for FeO
O2(:,6)=MC(:,6); %for MnO
O2(:,7)=MC(:,7); %for MgO
O2(:,8)=MC(:,8); %for CaO
O2(:,9)=MC(:,9)./2; %for Na2O
O2(:,10)=MC(:,10)./2; %H2O

O2_N=zeros(m,1);

%selects normalization scheme
for c=1:m
    if O2(c,10) > 0 %checks if H2O content is > 0
        O2_N(c,1)=(14)./sum(O2(c,1:10),2); %hydrous normalization
    else
        O2_N(c,1)=(12)./sum(O2(c,1:9),2); %anhydrous normalization
    end
end

%normalized moles of anions
N_Ox=O2.*O2_N;

%% atoms pfu

apfu=zeros(m,12);

apfu(:,1)=N_Ox(:,1)./2; %for SiO2
apfu(:,2)=N_Ox(:,2)./2; %for TiO2
apfu(:,3)=N_Ox(:,3).*(2/3); %for Al2O3
apfu(:,4)=N_Ox(:,4).*(2/3); %for Fe2O3
apfu(:,5)=N_Ox(:,5); %for FeO
apfu(:,6)=N_Ox(:,6); %for MnO
apfu(:,7)=N_Ox(:,7); %for MgO
apfu(:,8)=N_Ox(:,8); %for CaO
apfu(:,9)=N_Ox(:,9).*2; %for Na2O

%selects OH calculation scheme
for c=1:m
    if O2(c,10) > 0 %checks if H2O content is > 0
        apfu(c,10) = N_Ox(c,10).*2; %H, %hydrous normalization
    else
        apfu(c,10) = 4;%fixes OH the value to 4
    end
end

apfu(:,11)=sum(apfu(:,1:9),2); %calculations the cation total



%% Structural Formula

StrctFrm=zeros(m,13);

%T SITE
%Si
for c=1:m
    if apfu(c,1)<2.000
        StrctFrm(c,1)=apfu(c,1); %If Si < 2, then Si(T) = the measured Si content
    else
        StrctFrm(c,1)=2; %If Si is in excess, then Si(T) = 2
    end
end

%Al2O3 Layer (L2)
for c=1:m
    if apfu(c,3)<=3.000
        StrctFrm(c,2)=apfu(c,3); %If Al =< 3, then Al(3) = the measured Al content
    else
        StrctFrm(c,2)=3; %If Al>3 is in excess, then Al(L1) = 3
    end
end

%(Al, Ti, Fe3+) + (Fe, Mg, Mn) layer (L1)
StrctFrm(:,3)=apfu(:,3)-StrctFrm(:,2); %Al
StrctFrm(:,4)=apfu(:,2); %Ti
StrctFrm(:,5)=apfu(:,4); %Fe3+
StrctFrm(:,6)=apfu(:,5); %Fe2+
StrctFrm(:,7)=apfu(:,6); %Mn
StrctFrm(:,8)=apfu(:,7); %Mg
StrctFrm(:,9)=apfu(:,8); %Ca
StrctFrm(:,10)=apfu(:,9); %Na
StrctFrm(:,11)=sum(StrctFrm,2); %sum
StrctFrm(:,12)=apfu(:,10); %OH

XMg=zeros(m,1);
for c=1:m
    if StrctFrm(c,8) + StrctFrm(c,6) > 0 %checks if Mg content is > 0
        XMg(c,1)=StrctFrm(c,8)./(StrctFrm(c,8)+StrctFrm(c,6)); %XMg = Mg/(Mg+Fe)
    else
        %if Fe + Mg = 0, then the XMg calculation results in division by
        %zero, prevents NaN output
        XMg(c,1) = 0; 
    end
end

%% output

all=[StrctFrm XMg];

%limit on significant digits (eliminates rounding noise)
all(all<1e-6) = 0;
apfu(apfu<1e-6) = 0;

StrctFrm=array2table(all,'VariableNames',{'Si_T','Al_L2','Al_L1','Ti_L1','Fe3_L1','Fe2_L1','Mn2_L1','Mg_L1','Ca_L1','Na_L1','Cation_Sum','OH_H','O2_deficiency','XMg'});
apfu=array2table(apfu,'VariableNames',{'Si','Ti','Al','Fe3','Fe2','Mn2','Mg','Ca','Na','OH','Cation_Sum','O2_deficiency'});

end


