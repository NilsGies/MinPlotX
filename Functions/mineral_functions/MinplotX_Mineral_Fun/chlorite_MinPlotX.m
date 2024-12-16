%%chlorite Structural Formula
% last modified 03.07.2024

function [StrctFrm, apfu,options_definition]=chlorite_MinPlotX(data,headers,options)
%% define empty output variables
StrctFrm=[];
apfu=[];

%% options definition
options_definition=struct();

%Parameter 1
options_definition.UseKnownFe3.question='Is Fe3+/FeTotal Ratio known?';
options_definition.UseKnownFe3.Value=false;

%Parameter 2
options_definition.Fe3_ratio.question='Fe3+/FeTotal Ratio?';
options_definition.Fe3_ratio.Value=0;
options_definition.Fe3_ratio.limits=[0 1];

%Parameter 3 
options_definition.UseKnowntetra_Fe3.question='Use known Fe3+ Tetra Ratio?';
options_definition.UseKnowntetra_Fe3.Value=false;

%Parameter 4
options_definition.tetra_Fe3.question='Enter tetra_Fe3'; 
options_definition.tetra_Fe3.Value=0;
options_definition.tetra_Fe3.limits=[0 1]; 

%Parameter 5
options_definition.RecalcTotalFe.question='Recalculate FeO and Fe2O3?';
options_definition.RecalcTotalFe.Value=false;

%% Check if the data or header is empty 
if not(exist('headers','var')) || not(exist('data','var')) || isempty(data) || isempty(headers)
    return % exit the function and return options_definition
end

%% If no options 
if not(exist('options','var')) || isempty(options)
options=options_definition; % use default values generated in the options definition block
end

%% Start Calculation
[m,~]=size(data); %finds the x and y size of the input data matrix

%% Checks for Fe3+/Fetotal ratio

%assigns Fe3/Fetotal if it is included in the analysis
if isfield(options,'Fe3_ratio') && isfield(options,'UseKnownFe3') && options.UseKnownFe3.Value==true
    Fe_rat=options.Fe3_ratio.Value;
    options.RecalcTotalFe.Value=true;
else
    Fe_rat(:,1)=zeros (m,1); %Fe3+/FeTotal ratio is 
    options.RecalcTotalFe.Value=false;
end

%% Checks for tetra Fe3
if isfield(options,'UseKnowntetra_Fe3') && options.UseKnowntetra_Fe3.Value&& isfield(options,'tetra_Fe3') 
    tetra_Fe3=options.tetra_Fe3.Value;
else
    tetra_Fe3=zeros(m,1); %default Fe3+ on T ratio is 0´
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
%create columns of zeros if optional data are not included
MC=zeros(m,10);
apfu=zeros(m,10);

MC(:,1)=data(:,strcmp(headers,'SiO2'))./SiO2_mw; %for SiO2

%calculates for TiO2 if it is included in the analysis
if any(strcmp(headers,'TiO2'))
    MC(:,2)=data(:,strcmp(headers,'TiO2'))./TiO2_mw; %for TiO2
end

MC(:,3)=(data(:,strcmp(headers,'Al2O3'))./Al2O3_mw).*2; %for Al2O3

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
    MC(:,6)=Fe2_wtper(:,1)./FeO_mw;

elseif any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) %if FeO and Fe2O3 are both included as inputs

     MC(:,5)=(data(:,strcmp(headers,'Fe2O3'))./Fe2O3_mw).*2;
     MC(:,6)=data(:,strcmp(headers,'FeO'))./FeO_mw;
else
    if any(strcmp(headers,'FeO')) %if FeO is FeO total

        %first calculate the weight percent of Fe2O3 and FeO from FeO total
        Fe3_wtper(:,1)=(data(:,strcmp(headers,'FeO')).*Fe_rat)./((2*FeO_mw)./Fe2O3_mw); %calculates Fe2O3 from Fe3+ ratio and total FeO
        Fe2_wtper(:,1)=data(:,strcmp(headers,'FeO'))-(data(:,strcmp(headers,'FeO')).*Fe_rat); %FeO = Total FeO - Fe3+ in Total FeO

        %calculate moles of cations
        MC(:,5)=(Fe3_wtper(:,1)./Fe2O3_mw).*2;
        MC(:,6)=Fe2_wtper(:,1)./FeO_mw;

    else %if Fe2O3 is Fe2O3 total

        %first calculate the weight percent of Fe2O3 and FeO from Fe2O3 total
        Fe3_wtper(:,1)=data(:,strcmp(headers,'Fe2O3'))-(data(:,strcmp(headers,'Fe2O3')).*(1-Fe_rat)); %Fe2O3 = Total Fe2O3 - Fe2+ in Total Fe2O3 
        Fe2_wtper(:,1)=(data(:,strcmp(headers,'Fe2O3')).*(1-Fe_rat)).*((2*FeO_mw)./Fe2O3_mw); %calculates FeO from Fe3+ ratio and total Fe2O3

        %calculate moles of cations
        MC(:,5)=(Fe3_wtper(:,1)./Fe2O3_mw).*2;
        MC(:,6)=Fe2_wtper(:,1)./FeO_mw;

    end
end


%calculates for NiO if it is included in the analysis
if any(strcmp(headers,'NiO'))
    MC(:,7)=data(:,strcmp(headers,'NiO'))./NiO_mw; %for NiO
end

%calculates for MnO and/or Mn2O3
if any(strcmp(headers,'MnO')) && any(strcmp(headers,'Mn2O3')) %if MnO and Mn2O3 are both included as inputs

    %If MnO and Mn2O3 are both given, Mn2O3 is converted to MnO and
    %combined with MnO
    MC(:,8)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw)+data(:,strcmp(headers,'MnO')))./MnO_mw;
else
    if any(strcmp(headers,'MnO')) %if MnO is MnO total

        MC(:,8)=data(:,strcmp(headers,'MnO'))./MnO_mw;

    elseif any(strcmp(headers,'Mn2O3')) %if Mn2O3 total is given, converted to MnO total

        %calculate moles of cations
        MC(:,8)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw))./MnO_mw;

    else
        %Mn2O3 and MnO are not included, Mn2 = 0
    end
end

MC(:,9)=data(:,strcmp(headers,'MgO'))./MgO_mw; %for MgO

%calculates for H2O if it is included in the analysis
if any(strcmp(headers,'H2O'))
    MC(:,10)=(data(:,strcmp(headers,'H2O'))./H2O_mw).*2; %for H2O
end

%% Oxygen Units

O2(:,1)=MC(:,1).*2; %for SiO2
O2(:,2)=MC(:,2).*2; %for TiO2
O2(:,3)=MC(:,3).*(3/2); %for Al2O3
O2(:,4)=MC(:,4).*(3/2); %for Cr2O3
O2(:,5)=MC(:,5).*(3/2); %for Fe2O3
O2(:,6)=MC(:,6); %for FeO
O2(:,7)=MC(:,7); %for NiO
O2(:,8)=MC(:,8); %for MnO
O2(:,9)=MC(:,9); %for MgO
O2(:,10)=MC(:,10).*(1/2); %H2O

O2_N=zeros(m,1);

%selects normalization scheme
for c=1:m
    if O2(c,10) > 0 %checks if H2O content is > 0
        O2_N(c,1)=(18)./sum(O2(c,1:10),2); %hydrous normalization
    else
        O2_N(c,1)=(14)./sum(O2(c,1:9),2); %anhydrous normalization
    end
end

%normalized moles of anions
N_Ox=O2.*O2_N;

%% atoms pfu

apfu(:,1)=N_Ox(:,1)./2; %Si
apfu(:,2)=N_Ox(:,2)./2; %Ti
apfu(:,3)=N_Ox(:,3).*(2/3); %Al
apfu(:,4)=N_Ox(:,4).*(2/3); %Cr
apfu(:,5)=N_Ox(:,5).*(2/3); %Fe3
apfu(:,6)=N_Ox(:,6); %Fe2
apfu(:,7)=N_Ox(:,7); %Ni
apfu(:,8)=N_Ox(:,8); %Mn
apfu(:,9)=N_Ox(:,9); %Mg

%selects OH calculation scheme

for c=1:m
    if O2(c,10) > 0 %checks if H2O content is > 0
        apfu(c,10) = N_Ox(c,10).*2; %H, %hydrous normalization
    else
        apfu(c,10) = 8; %fixes H the value to 8
    end
end

apfu(:,11)=sum(apfu(:,1:9),2); %calculations the cation total

%% Structural formula

StrctFrm=zeros(m,15);

%T site
StrctFrm(:,1)=apfu(:,1); %Si (T)

StrctFrm(:,3)=apfu(:,5).*tetra_Fe3; %Fe3+ (T)

%Al (T)
for c=1:m
    if 4-StrctFrm(c,1)-StrctFrm(c,3) > 0
        %checks if T site is filled by Si + Fe3
        if 4-StrctFrm(c,1)-StrctFrm(c,3) > apfu(c,3)
            %if the T leftover space in the T site is greater than the Al
            %content all Al goes on T (unlikely for chlorite!)
            StrctFrm(c,2)=apfu(c,3);
        else
            StrctFrm(c,2)=4-StrctFrm(c,1)-StrctFrm(c,3);
        end
    else
        %If the T site is filled by Si + Fe3 then Al(T) = 0
    end
end

%T site sum
StrctFrm(:,4)=sum(StrctFrm(:,1:1:3),2); %should sum to 4

%Al (M)
StrctFrm(:,5)=apfu(:,3)-StrctFrm(:,2); %Al total - Al on the T site
StrctFrm(:,6)=apfu(:,2); %Ti (M)
StrctFrm(:,7)=apfu(:,4); %Cr (M)
StrctFrm(:,8)=apfu(:,5)-StrctFrm(:,3); %Fe3 (M) = Fe3total - Fe3(T)
StrctFrm(:,9)=apfu(:,6); %Fe2 (M)
StrctFrm(:,10)=apfu(:,7); %Ni (M)
StrctFrm(:,11)=apfu(:,8); %Mn2 (M)
StrctFrm(:,12)=apfu(:,9); %Mg (M)

%Vacancies on M1
for c=1:m
    if (StrctFrm(c,5)-StrctFrm(c,2))/2 + (StrctFrm(c,8)/2)> 0 %Vacancies on M1 = (Alvi-Aliv)/2 + Fe3vi/2 is > 0
        StrctFrm(c,13)=(StrctFrm(c,5)-StrctFrm(c,2))/2 + (StrctFrm(c,8)/2); %Vacancies on M1 = (Alvi-Aliv)/2
    else
        StrctFrm(c,13) = 0; 
    end
end

StrctFrm(:,14)=sum(StrctFrm(:,1:1:12),2); %Total Cations
StrctFrm(:,15)=apfu(:,10); %OH

XMg=zeros(m,1);
for c=1:m
    if StrctFrm(c,12) + StrctFrm(c,9) > 0 %checks if Mg content is > 0
        XMg(c,1)=StrctFrm(c,12)./(StrctFrm(c,12)+StrctFrm(c,9)); %XMg = Mg/(Mg+Fe)
    else
        %if Fe + Mg = 0, then the XMg calculation results in division by
        %zero, prevents NaN output
        XMg(c,1) = 0; 
    end
end

all=[StrctFrm XMg];

%limit on significant digits (eliminates rounding noise)
all(all<1e-6) = 0;
apfu(apfu<1e-6) = 0;

StrctFrm=array2table(all,'VariableNames',{'Si_T','Al_T','Fe3_T','Sum_T','Al_M','Ti_M','Cr_M','Fe3_M','Fe2_M','Ni_M','Mn2_M','Mg_M','vac_M1','Cation_Sum','OH','XMg'});

% StrctFrm=array2table(all,'VariableNames',{'Si_T','Al_T','Fe3_T','Sum_T','Al_M','Ti_M','Cr_M','Fe3_M','Fe2_M','Ni_M','Mn2_M','Mg_M','vac_M1','Cation_Sum','OH'});
% StrctFrm.XMg(:)=0;
% StrctFrm.XMg(:)=StrctFrm.Mg_M(StrctFrm.Mg_M>0)/(StrctFrm.Mg_M(StrctFrm.Mg_M>0)+StrctFrm.Fe2_M(StrctFrm.Mg_M>0));


apfu=array2table(apfu,'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Ni','Mn2','Mg','OH','Cation_Sum'});

%apfu.Fe3_FeTot=apfu.Fe3./(apfu.Fe3+apfu.Fe2);

end



