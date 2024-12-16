%lawsonite structural formula
% last modified 18.07.2024

function [StrctFrm, apfu,options_definition]=lawsonite_MinPlotX(data,headers,options)
%% define empty output variables
StrctFrm=[];
apfu=[];

%% options definition
options_definition=struct();
%Parameter 1
options_definition.UseKnownFe3.question='Do you want to use a known Fe3 ratio';
options_definition.UseKnownFe3.Value=true;

%Parameter 2
options_definition.Fe3_ratio.question='Enter Fe3_ratio';
options_definition.Fe3_ratio.Value=1;
options_definition.Fe3_ratio.limits=[0 1];

%Parameter 3
options_definition.UseKnownMn3.question='Do you want to use a known Mn3 ratio';
options_definition.UseKnownMn3.Value=false;

%Parameter 4
options_definition.Mn3_ratio.question='Enter Mn3_ratio';
options_definition.Mn3_ratio.Value=1;
options_definition.Mn3_ratio.limits=[0 1];

%Parameter 5
options_definition.RecalcTotalFe.question='Recalculate FeO and Fe2O3?';
options_definition.RecalcTotalFe.Value=false;

%Parameter 6
options_definition.RecalcTotalMn.question='Recalculate MnO and Mn2O3?';
options_definition.RecalcTotalMn.Value=false;

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

Opfu=8; %oxygens per formula unit
%% Checks for Fe3+/Fetotal ratio

%assigns Fe3/Fetotal if it is included in the analysis

%assigns Fe3/Fetotal if it is included in the analysis
if isfield(options,'Fe3_ratio') && isfield(options,'UseKnownFe3') && options.UseKnownFe3.Value==true
    Fe_rat=options.Fe3_ratio.Value;
    options.RecalcTotalFe.Value=true;
else
    Fe_rat(:,1)=ones(m,1); %Fe3+/FeTotal ratio is 1´
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

MC=zeros(m,14);%create columns of zeros if optional data are not included

MC(:,1)=data(:,strcmp(headers,'SiO2'))./SiO2_mw; % for SiO2

%calculates for TiO2 if it is included in the analysis
if any(strcmp(headers,'TiO2'))
    MC(:,2)=data(:,strcmp(headers,'TiO2'))./TiO2_mw; %for TiO2
end

MC(:,3)=(data(:,strcmp(headers,'Al2O3'))./Al2O3_mw).*2; % for Al2O3


%calculates for Cr2O3 if it is included in the analysis
if any(strcmp(headers,'Cr2O3'))
    MC(:,4)=(data(:,strcmp(headers,'Cr2O3'))./Cr2O3_mw).*2; %for Cr2O3
end

%calculates for FeO and Fe2O3

if any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) && isfield(options,'RecalcTotalFe') && options.RecalcTotalFe.Value==true

    %recalculate total FeO based on FeO+Fe2O3:
    FeOT=(data(:,strcmp(headers,'Fe2O3')).*((2*FeO_mw)./Fe2O3_mw)+data(:,strcmp(headers,'FeO')));

    %first calculate the weight percent of Fe2O3 and FeO from FeO total
    Fe3_wtper(:,1)=(FeOT.*Fe_rat)./((2*FeO_mw)./Fe2O3_mw); %calculates Fe2O3 from Fe3+ ratio and total FeO
    Fe2_wtper(:,1)=FeOT-(FeOT.*Fe_rat); %FeO = Total FeO - Fe3+ in Total FeO

    %calculate moles of cations
    MC(:,5)=(Fe3_wtper(:,1)./Fe2O3_mw).*2;
    MC(:,7)=Fe2_wtper(:,1)./FeO_mw;

elseif any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) %if FeO and Fe2O3 are both included as inputs

     MC(:,5)=(data(:,strcmp(headers,'Fe2O3'))./Fe2O3_mw).*2;
     MC(:,7)=data(:,strcmp(headers,'FeO'))./FeO_mw;
else
    if any(strcmp(headers,'FeO')) %if FeO is FeO total

        %first calculate the weight percent of Fe2O3 and FeO from FeO total
        Fe3_wtper(:,1)=(data(:,strcmp(headers,'FeO')).*Fe_rat)./((2*FeO_mw)./Fe2O3_mw); %calculates Fe2O3 from Fe3+ ratio and total FeO
        Fe2_wtper(:,1)=data(:,strcmp(headers,'FeO'))-(data(:,strcmp(headers,'FeO')).*Fe_rat); %FeO = Total FeO - Fe3+ in Total FeO

        %calculate moles of cations
        MC(:,5)=(Fe3_wtper(:,1)./Fe2O3_mw).*2;
        MC(:,7)=Fe2_wtper(:,1)./FeO_mw;

    elseif any(strcmp(headers,'Fe2O3'))%if Fe2O3 is Fe2O3 total

        %first calculate the weight percent of Fe2O3 and FeO from Fe2O3 total
        Fe3_wtper(:,1)=data(:,strcmp(headers,'Fe2O3'))-(data(:,strcmp(headers,'Fe2O3')).*(1-Fe_rat)); %Fe2O3 = Total Fe2O3 - Fe2+ in Total Fe2O3 
        Fe2_wtper(:,1)=(data(:,strcmp(headers,'Fe2O3')).*(1-Fe_rat)).*((2*FeO_mw)./Fe2O3_mw); %calculates FeO from Fe3+ ratio and total Fe2O3

        %calculate moles of cations
        MC(:,5)=(Fe3_wtper(:,1)./Fe2O3_mw).*2;
        MC(:,7)=Fe2_wtper(:,1)./FeO_mw;

    else
        %Fe2O3 and FeO are not included, Fe2 and Fe3 = 0
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
    MC(:,6)=(Mn3_wtper(:,1)./Mn2O3_mw).*2;
    MC(:,8)=Mn2_wtper(:,1)./MnO_mw;

elseif any(strcmp(headers,'MnO')) && any(strcmp(headers,'Mn2O3')) %if MnO and Mn2O3 are both included as inputs

     MC(:,6)=(data(:,strcmp(headers,'Mn2O3'))./Mn2O3_mw).*2;
     MC(:,8)=data(:,strcmp(headers,'MnO'))./MnO_mw;
else
    if any(strcmp(headers,'MnO')) %if MnO is MnO total

        %first calculate the weight percent of Mn2O3 and MnO from MnO total
        Mn3_wtper(:,1)=(data(:,strcmp(headers,'MnO')).*Mn_rat)./((2*MnO_mw)./Mn2O3_mw); %calculates Mn2O3 from Mn3+ ratio and total MnO
        Mn2_wtper(:,1)=data(:,strcmp(headers,'MnO'))-(data(:,strcmp(headers,'MnO')).*Mn_rat); %MnO = Total MnO - Mn3+ in Total MnO

        %calculate moles of cations
        MC(:,6)=(Mn3_wtper(:,1)./Mn2O3_mw).*2;
        MC(:,8)=Mn2_wtper(:,1)./MnO_mw;

    elseif any(strcmp(headers,'Mn2O3')) %if Mn2O3 is Mn2O3 total

        %first calculate the weight percent of Fe2O3 and FeO from Fe2O3 total
        Mn3_wtper(:,1)=data(:,strcmp(headers,'Mn2O3'))-(data(:,strcmp(headers,'Mn2O3')).*(1-Mn_rat)); %Mn2O3 = Total Mn2O3 - Mn2+ in Total Mn2O3 
        Mn2_wtper(:,1)=(data(:,strcmp(headers,'Mn2O3')).*(1-Mn_rat)).*((2*MnO_mw)./Mn2O3_mw); %calculates MnO from Mn3+ ratio and total Mn2O3

        %calculate moles of cations
        MC(:,6)=(Mn3_wtper(:,1)./Mn2O3_mw).*2;
        MC(:,8)=Mn2_wtper(:,1)./MnO_mw;
    else
        %Mn2O3 and MnO are not included, Mn2 and Mn3 = 0
    end
end

if any(strcmp(headers,'MgO'))
    MC(:,9)=data(:,strcmp(headers,'MgO'))./MgO_mw; % for MgO
end

MC(:,10)=data(:,strcmp(headers,'CaO'))./CaO_mw; % for CaO


if any(strcmp(headers,'SrO'))
    MC(:,11)=data(:,strcmp(headers,'SrO'))./SrO_mw; % for SrO
end

%calculates for Na2O if it is included in the analysis
if any(strcmp(headers,'Na2O'))
    MC(:,12)=(data(:,strcmp(headers,'Na2O'))./Na2O_mw).*2; %for Na2O
end

%calculates for K2O if it is included in the analysis
if any(strcmp(headers,'K2O'))
    MC(:,13)=(data(:,strcmp(headers,'K2O'))./K2O_mw).*2; %for K2O
end

%calculates for H2O if it is included in the analysis
if any(strcmp(headers,'H2O'))
    MC(:,14)=(data(:,strcmp(headers,'H2O'))./H2O_mw).*2; %for H2O
end


%% Calculate Oxygen Units

O2=zeros(m,14);

O2(:,1)=MC(:,1).*2; %for SiO2
O2(:,2)=MC(:,2).*2; %for TiO2
O2(:,3)=MC(:,3).*(3/2); %for Al2O3
O2(:,4)=MC(:,4).*(3/2); %for Cr2O3
O2(:,5)=MC(:,5).*(3/2); %for Fe2O3
O2(:,6)=MC(:,6).*(3/2); %for Mn2O3
O2(:,7)=MC(:,7); %for FeO
O2(:,8)=MC(:,8); %for MnO
O2(:,9)=MC(:,9); %for MgO
O2(:,10)=MC(:,10); %for CaO
O2(:,11)=MC(:,11); %for SrO
O2(:,12)=MC(:,12).*(1/2); %for Na2O
O2(:,13)=MC(:,13).*(1/2); %for K2O
O2(:,14)=MC(:,14).*(1/2); %H2O

O2_N=zeros(m,1);

%selects normalization scheme
for c=1:m
    if O2(c,14) > 0 %checks if H2O content is > 0
        O2_N(c,1)=(10)./sum(O2(c,1:14),2); %hydrous normalization
    else
        O2_N(c,1)=(8)./sum(O2(c,1:13),2); %anhydrous normalization
    end
end

%normalized moles of anions
N_Ox=O2.*O2_N;

%% atoms pfu

apfu=zeros(m,15);

apfu(:,1)=N_Ox(:,1).*(1/2); %for SiO2
apfu(:,2)=N_Ox(:,2).*(1/2); %for TiO2
apfu(:,3)=N_Ox(:,3).*(2/3); %for Al2O3
apfu(:,4)=N_Ox(:,4).*(2/3); %for Cr2O3
apfu(:,5)=N_Ox(:,5).*(2/3); %for Fe2O3
apfu(:,6)=N_Ox(:,6).*(2/3); %for Mn2O3
apfu(:,7)=N_Ox(:,7); %for FeO
apfu(:,8)=N_Ox(:,8); %for MnO
apfu(:,9)=N_Ox(:,9); %for MgO
apfu(:,10)=N_Ox(:,10); %for CaO
apfu(:,11)=N_Ox(:,11); %for SrO
apfu(:,12)=N_Ox(:,12).*2; %for Na2O
apfu(:,13)=N_Ox(:,13).*2; %for K2O

%selects OH calculation scheme
for c=1:m
    if O2(c,14) > 0 %checks if H2O content is > 0
        apfu(c,14) = N_Ox(c,14).*2; %H, %hydrous normalization
    else
        apfu(c,14) = 4; %fixes H the value to 4
    end
end

apfu(:,15)=sum(apfu(:,1:13),2); %calculations the total, which should be close to 8


%% Structural Formula
%T site - Si, Al (sums to 2)
%M sites - Ti, Al, Cr3+, Fe3+, Mn3+, Mg2+, Fe2+, Mn2+ (Sums to 2)
%A site - Ca, Sr, Na, K, Mg, Fe2+, Mn2+ (Sums to 1)

StrctFrm=zeros(m,21);

%T SITE
%Si
for c=1:m
    if apfu(c,1)<2.000
        StrctFrm(c,1)=apfu(c,1); %If Si < 2, then Si(T) = the measured Si content
    else
        StrctFrm(c,1)=2; %If Si is in excess, then Si(T) = 2
    end
end

%Al(T)
for c=1:m
    if 2-StrctFrm(c,1)>0 %Is 2-Si > 0? If y, then some Al goes into T
        StrctFrm(c,2)=2-StrctFrm(c,1);
    else
        StrctFrm(c,2)=0; %if Si=2, then no Al goes into T
    end
end


StrctFrm(:,3)=StrctFrm(:,2)+StrctFrm(:,1); %T site sum


StrctFrm(:,4)=apfu(:,3)-StrctFrm(:,2); %Al (M)
StrctFrm(:,5)=apfu(:,2);%Ti (M)
StrctFrm(:,6)=apfu(:,4); %Cr (M)
StrctFrm(:,7)=apfu(:,5); %Fe3+ (M)
StrctFrm(:,8)=apfu(:,6); %Mn3+ (M)

%Mg (M)
for c=1:m
    if (StrctFrm(c,4)+StrctFrm(c,5)+StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8))<2.000 %Mg is only considered if the M site is not yet filled
        if (2-(StrctFrm(c,4)+StrctFrm(c,5)+StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8))) > apfu(c,9)
            StrctFrm(c,9)=apfu(c,9); %all Mg goes into the M site
        else
            StrctFrm(c,9)=2-(StrctFrm(c,4)+StrctFrm(c,5)+StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8)); % only some Mg goes into the M site
        end
    else
        StrctFrm(c,9)=0; % no Mg goes into the M site (it's already filled)
    end
end

%Fe2+ (M)
for c=1:m
    if (StrctFrm(c,4)+StrctFrm(c,5)+StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8)+StrctFrm(c,9))<2.000 %Fe2+ is only considered if the M site is not yet filled
        if (2-(StrctFrm(c,4)+StrctFrm(c,5)+StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8)+StrctFrm(c,9))) > apfu(c,7)
            StrctFrm(c,10)=apfu(c,7); %all Fe2+ goes into the M site
        else
            StrctFrm(c,10)=2-(StrctFrm(c,4)+StrctFrm(c,5)+StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8)+StrctFrm(c,9)); % only some Fe2+ goes into the M site
        end
    else
        StrctFrm(c,10)=0; % no Fe2+ goes into the octahedral site (it's already filled)
    end
end

%Mn2+ (M)
for c=1:m
    if (StrctFrm(c,4)+StrctFrm(c,5)+StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8)+StrctFrm(c,9)+StrctFrm(c,10))<2.000 %Mn is only considered if the M site is not yet filled
        if (2-(StrctFrm(c,4)+StrctFrm(c,5)+StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8)+StrctFrm(c,9)+StrctFrm(c,10))) > apfu(c,8)
            StrctFrm(c,11)=apfu(c,8); %all Mn2+ goes into M site
        else
            StrctFrm(c,11)=2-(StrctFrm(c,4)+StrctFrm(c,5)+StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8)+StrctFrm(c,9)+StrctFrm(c,10)); % only some Mn goes into the M site
        end
    else
        StrctFrm(c,11)=0; % no Mn goes into the M site (it's already filled)
    end
end

StrctFrm(:,12)=sum(StrctFrm(:,4:11),2); %M site sum

%A site 
StrctFrm(:,13)=apfu(:,9)-StrctFrm(:,9); %Mg (A) = Mg total - Mg (M)
StrctFrm(:,14)=apfu(:,7)-StrctFrm(:,10); %Fe2 (A) = Fe2 total - Fe2 (M)
StrctFrm(:,15)=apfu(:,8)-StrctFrm(:,11); %Mn2 (A) = Mn2 total - Mn2 (M)
StrctFrm(:,16)=apfu(:,10); %Ca (A)
StrctFrm(:,17)=apfu(:,11); %Sr (A)
StrctFrm(:,18)=apfu(:,12); %Na (A)
StrctFrm(:,19)=apfu(:,13); %K (A)
StrctFrm(:,20)=sum(StrctFrm(:,13:19),2); %A site sum

%OH site
for c=1:m
    if apfu(c,14) > 2
        %if H is > 2 then OH site is full
        StrctFrm(c,21)=2; %OH
    else
        %otherwise all H is on the OH site
        StrctFrm(c,21)=apfu(c,14); %OH
    end
end

%H2O site
for c=1:m
    if apfu(c,14) > 2
        %if H is > 2, the remainder of H goes to H2O, divided by 2
        %to give moles of H2O from moles of H
        StrctFrm(c,22)=(apfu(c,14)-StrctFrm(c,21))/2; %H2O
    else
        %if all H goes into OH site, 
        StrctFrm(c,22)=0; %H2O
    end
end


Endmembers(:,1)=(StrctFrm(:,4)./StrctFrm(:,12)).*(StrctFrm(:,16)./(StrctFrm(:,20))); %XLws = (Al/(M site sum)) * (Ca/(A site sum))

all=[StrctFrm Endmembers];

%limit on significant digits (eliminates rounding noise)
all(all<1e-6) = 0;
apfu(apfu<1e-6) = 0;

StrctFrm=array2table(all,'VariableNames',{'Si_T','Al_T','T_sum','Al_M','Ti_M','Cr_M','Fe3_M','Mn3_M','Mg_M','Fe2_M','Mn2_M','M_Sum','Mg_A','Fe2_A','Mn2_A','Ca_A','Sr_A','Na_A','K_A','A_Sum','OH','H2O','Xlws'});
apfu=array2table(apfu,'VariableNames',{'Si','Ti','Al','Cr','Fe3','Mn3','Fe2','Mn2','Mg','Ca','Sr','Na','K','OH','Cation_Sum'});

end


