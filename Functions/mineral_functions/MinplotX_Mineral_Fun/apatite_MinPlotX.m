%apatite structural formula
% last modified 19.07.2024

function [StrctFrm,APFU,options_definition]=apatite_MinPlotX(data,headers,options)
%% define empty output variables
StrctFrm=[];
APFU=[];

%% options definition
%Parameter 6
options_definition.carbon_species.question='Is carbon analyzed??';
options_definition.carbon_species.Value=false;

%Parameter 7
options_definition.Z_ratio.question='Enter ratio of C on the anion site';
options_definition.Z_ratio.Value=0;
options_definition.Z_ratio.limits=[0 1];

%Parameter 1
options_definition.Sulfur_species.question='Is sulfur analyzed??';
options_definition.Sulfur_species.Value=false;

%Parameter 2
options_definition.S6_ratio.question='Enter S6+/S';
options_definition.S6_ratio.Value=1;
options_definition.S6_ratio.limits=[0 1];

%Parameter 3
options_definition.S4_ratio.question='Enter S4+/S';
options_definition.S4_ratio.Value=0;
options_definition.S4_ratio.limits=[0 1];

%Parameter 4
options_definition.S1minus_ratio.question='Enter S1-/S';
options_definition.S1minus_ratio.Value=0;
options_definition.S1minus_ratio.limits=[0 1];

%Parameter 5
options_definition.S2minus_ratio.question='Enter S2-/S';
options_definition.S2minus_ratio.Value=0;
options_definition.S2minus_ratio.limits=[0 1];

%% Check if the data or header is empty 
if not(exist('headers','var')) || not(exist('data','var')) || isempty(data) || isempty(headers)
    return % exit the function and return options_definition
end

%% If no options 
if not(exist('options','var')) || isempty(options)
options=options_definition; % use default values generated in the options definition block
end
%% start calculation

%choose the speciation of sulfur:
if options.Sulfur_species.Value==true

    Sratio = [options.S6_ratio.Value;
        options.S4_ratio.Value;
        options.S1minus_ratio.Value;
        options.S2minus_ratio.Value];

else
    Sratio =[options_definition.S6_ratio.Value;
        options_definition.S4_ratio.Value;
        options_definition.S1minus_ratio.Value;
        options_definition.S2minus_ratio.Value];
end

%choose the allocation of C:
if options.carbon_species.Value==true
    Z_ratio = [options.Z_ratio.Value];
else
    Z_ratio = [options_definition.Z_ratio.Value];
end


[m,~]=size(data); %finds the x and y size of the input data matrix

%finds the column position oxide headers in any order to assign data to the correct array
%positions

%% Molecular weights

SiO2_mw=60.083;
TiO2_mw=79.865;
Al2O3_mw=101.961;
Cr2O3_mw=151.989;
Fe2O3_mw=159.6874;
Y2O3_mw=225.809;
Ce2O3_mw=328.229;
La2O3_mw=325.80794;
NiO_mw=74.692;
ZnO_mw=81.381;
FeO_mw=71.8442;
MnO_mw=70.937;
Mn2O3_mw=157.873;
MgO_mw=40.304;
CaO_mw=56.0774;
Na2O_mw=61.979;
K2O_mw=94.195;
BaO_mw=153.329;
SrO_mw=103.619;
SO3_mw=80.0594;
F_mw=18.998;
Cl_mw=35.45;
P2O5_mw=141.942524;
H2O_mw=18.015;
CO2_mw=44.009;

%% Moles of oxides

MC=zeros(m,19); %create columns of zeros if optional data are not included

%calculates for SiO2 if it is included in the analysis
if any(strcmp(headers,'SiO2'))
    MC(:,1)=data(:,strcmp(headers,'SiO2'))./SiO2_mw; %for SiO2
end

MC(:,2)=(data(:,strcmp(headers,'P2O5'))./P2O5_mw).*2; %for P2O5

%calculates for TiO2 if it is included in the analysis
if any(strcmp(headers,'TiO2'))
    MC(:,3)=data(:,strcmp(headers,'TiO2'))./TiO2_mw; %for TiO2
end

%calculates for Al2O3 if it is included in the analysis
if any(strcmp(headers,'Al2O3'))
    MC(:,4)=(data(:,strcmp(headers,'Al2O3'))./Al2O3_mw).*2; %for Al2O3
end

%calculates for FeO if it is included in the analysis
if any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) %if FeO and Fe2O3 are both included as inputs 

     %If FeO and Fe2O3 are both given, Fe2O3 is converted to FeO and
     %combined with FeO
     MC(:,5)=(data(:,strcmp(headers,'Fe2O3')).*((2*FeO_mw)./Fe2O3_mw)+data(:,strcmp(headers,'FeO')))./FeO_mw;
else
    if any(strcmp(headers,'FeO')) %if FeO is FeO total

        MC(:,5)=data(:,strcmp(headers,'FeO'))./FeO_mw;

    elseif any(strcmp(headers,'Fe2O3'))%if Fe2O3 total is given, converted to FeO total

        %calculate moles of cations
        MC(:,5)=(data(:,strcmp(headers,'Fe2O3')).*((2*FeO_mw)./Fe2O3_mw))./FeO_mw;
    else
        %Fe2O3 and FeO are not included, Fe2 = 0
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

%calculates for MgO if it is included in the analysis
if any(strcmp(headers,'MgO'))
    MC(:,7)=data(:,strcmp(headers,'MgO'))./MgO_mw; %for MgO
end

MC(:,8)=data(:,strcmp(headers,'CaO'))./CaO_mw; %for CaO

%calculates for BaO if it is included in the analysis
if any(strcmp(headers,'BaO'))
    MC(:,9)=data(:,strcmp(headers,'BaO'))./BaO_mw; %for BaO
end

%calculates for SrO if it is included in the analysis
if any(strcmp(headers,'SrO'))
    MC(:,10)=data(:,strcmp(headers,'SrO'))./SrO_mw; %for SrO
end

%calculates for K2O if it is included in the analysis
if any(strcmp(headers,'K2O'))
    MC(:,11)=(data(:,strcmp(headers,'K2O'))./K2O_mw).*2; %for K2O
end

%calculates for Na2O if it is included in the analysis
if any(strcmp(headers,'Na2O'))
    MC(:,12)=(data(:,strcmp(headers,'Na2O'))./Na2O_mw).*2; %for Na2O
end

%calculates for Ce2O3 if it is included in the analysis
if any(strcmp(headers,'Ce2O3'))
    MC(:,13)=(data(:,strcmp(headers,'Ce2O3'))./Ce2O3_mw).*2; %for Ce2O3
end

%calculates for La2O3 if it is included in the analysis
if any(strcmp(headers,'La2O3'))
    MC(:,14)=(data(:,strcmp(headers,'La2O3'))./La2O3_mw).*2; %for La2O3
end

if any(strcmp(headers,'SO3')) 
MC(:,15)=data(:,strcmp(headers,'SO3'))./SO3_mw; %for SO3 
end 

if any(strcmp(headers,'F')) 
MC(:,16)=data(:,strcmp(headers,'F'))./F_mw; %for F
end 

if any(strcmp(headers,'Cl')) 
MC(:,17)=data(:,strcmp(headers,'Cl'))./Cl_mw; %for Cl
end 

%calculates for H2O if it is included in the analysis
if any(strcmp(headers,'H2O'))
    MC(:,18)=(data(:,strcmp(headers,'H2O'))./H2O_mw).*2; %for H2O
end

%calculates for CO2 if it is included in the analysis
if any(strcmp(headers,'CO2'))
    MC(:,19)=(data(:,strcmp(headers,'CO2'))./CO2_mw); %for CO2
end


%% Moles of O2 units
O2(:,1)=MC(:,1).*2; %SiO2
O2(:,2)=MC(:,2).*(5/2); %P2O5
O2(:,3)=MC(:,3).*2; %TiO2
O2(:,4)=MC(:,4).*(3/2); %Al2O3
O2(:,5)=MC(:,5); %FeO
O2(:,6)=MC(:,6); %MnO
O2(:,7)=MC(:,7); %MgO
O2(:,8)=MC(:,8); %CaO
O2(:,9)=MC(:,9); %BaO
O2(:,10)=MC(:,10); %SrO
O2(:,11)=MC(:,11).*(1/2); %K2O
O2(:,12)=MC(:,12).*(1/2); %Na2O
O2(:,13)=MC(:,13).*(3/2); %Ce2O3
O2(:,14)=MC(:,14).*(3/2); %La2O3
O2(:,15)=(MC(:,15).*Sratio(1,1)).*3; %S6+
O2(:,16)=(MC(:,15).*Sratio(2,1)).*2; %S4+
O2(:,17)=MC(:,15).*Sratio(3,1); %S1-
O2(:,18)=MC(:,15).*Sratio(4,1); %S2-
O2(:,19)=MC(:,16); %F
O2(:,20)=MC(:,17); %Cl
O2(:,21)=MC(:,18).*(1/2); %H2O
O2(:,22)=MC(:,19).*3; %CO3

if any(strcmp(headers,'H2O')) %if H2O is known

    O2_Sum=sum(O2(:,1:22),2)-0.5.*(O2(:,17)+O2(:,18)+O2(:,19)+O2(:,20)); %sum of O2, excluding F, S1-, S2-, and Cl

    O2_N=(26)./O2_Sum; %normalization factor

    %normalized moles of anions
    N_Ox=O2.*O2_N;

else %method of Ketcham (2015) if H2O is not known

    O2_CO3_P=O2(:,22)-Z_ratio.*O2(:,22); %CO3 substituting for P2O5
    O2_Sum=sum(O2(:,1:16),2)+O2_CO3_P; %sum of O2, excluding F, S1-, S2-, Cl, and CO3 on the anion site

    O2_N=(25)./O2_Sum; %normalization factor

    %normalized moles of anions
    N_Ox=O2.*O2_N;
end

%% atoms pfu

APFU(:,1)=N_Ox(:,1).*(1/2); %SiO2
APFU(:,2)=N_Ox(:,2).*(2/5); %P2O5
APFU(:,3)=N_Ox(:,3).*(1/2); %TiO2
APFU(:,4)=N_Ox(:,4).*(2/3); %Al2O3
APFU(:,5)=N_Ox(:,5); %FeO
APFU(:,6)=N_Ox(:,6); %MnO
APFU(:,7)=N_Ox(:,7); %MgO
APFU(:,8)=N_Ox(:,8); %CaO
APFU(:,9)=N_Ox(:,9); %BaO
APFU(:,10)=N_Ox(:,10); %SrO
APFU(:,11)=N_Ox(:,11).*2; %K2O
APFU(:,12)=N_Ox(:,12).*2; %Na2O
APFU(:,13)=N_Ox(:,13).*(2/3); %Ce2O3
APFU(:,14)=N_Ox(:,14).*(2/3); %La2O3
APFU(:,15)=N_Ox(:,15).*(1/3); %S6+
APFU(:,16)=N_Ox(:,16).*(1/2); %S4+
APFU(:,17)=N_Ox(:,17); %S1-
APFU(:,18)=N_Ox(:,18); %S2-
APFU(:,19)=N_Ox(:,19); %F
APFU(:,20)=N_Ox(:,20); %Cl
APFU(:,22)=N_Ox(:,22).*(1/3); %CO3

if any(strcmp(headers,'H2O')) %if H2O is known
    APFU(:,21)=N_Ox(:,21).*2; %H2O
else %if H2O is estimated
    O2_CO3_Z=Z_ratio.*N_Ox(:,22); %CO3 substituting for on Z
    APFU(:,21)=2-(N_Ox(:,20)+N_Ox(:,19)+N_Ox(:,17)+(2.*N_Ox(:,18))+((1/3).*O2_CO3_Z)); %H = 2 - (F + Cl + S^1- + 2 * S^2- + CO3^2-)
end

%% Structural Formula

%T site
StrctFrm(:,1)=APFU(:,2); %P (T)
StrctFrm(:,2)=APFU(:,1); %Si (T)
StrctFrm(:,3)=APFU(:,15); %S6+ (T)
StrctFrm(:,4)=APFU(:,16); %S4+ (T)
StrctFrm(:,5)=APFU(:,22)-APFU(:,22).*Z_ratio; %C4+ (T)

%M sites
StrctFrm(:,6)=APFU(:,3); %Ti (M)
StrctFrm(:,7)=APFU(:,4); %Al (M)
StrctFrm(:,8)=APFU(:,5); %Fe (M)
StrctFrm(:,9)=APFU(:,6); %Mn (M)
StrctFrm(:,10)=APFU(:,7); %Mg (M)
StrctFrm(:,11)=APFU(:,8); %Ca (M)
StrctFrm(:,12)=APFU(:,9); %Ba (M)
StrctFrm(:,13)=APFU(:,10); %Sr (M)
StrctFrm(:,14)=APFU(:,11); %K (M)
StrctFrm(:,15)=APFU(:,12); %Na (M)
StrctFrm(:,16)=APFU(:,13); %Ce (M)
StrctFrm(:,17)=APFU(:,14); %La (M)
StrctFrm(:,18)=sum(StrctFrm(:,1:17),2); %cation sum

%Z site
StrctFrm(:,19)=APFU(:,21); %OH
StrctFrm(:,20)=APFU(:,19); %F
StrctFrm(:,21)=APFU(:,20); %Cl
StrctFrm(:,22)=APFU(:,17); %S1-
StrctFrm(:,23)=APFU(:,18); %S2-
StrctFrm(:,24)=APFU(:,22).*Z_ratio; %C4+
StrctFrm(:,25)=sum(StrctFrm(:,19:24),2); %anion sum

%limit on significant digits (eliminates rounding noise)
StrctFrm(StrctFrm<1e-6) = 0;
APFU(APFU<1e-6) = 0;

StrctFrm=array2table(StrctFrm,'VariableNames',{'P_T','Si_T','S6_T','S4_T','C4_T','Ti_M','Al_M','Fe2_M' ...
    ,'Mn2_M','Mg_M','Ca_M','Ba_M','Sr_M','K_M','Na_M','Ce_M','La_M','Cation_Sum','OH_Z', ...
    'F_Z','Cl_Z','S1-_Z','S2-_Z','C4_Z','Anion_Sum'});

APFU=array2table(APFU,'VariableNames',{'Si','P','Ti','Al','Fe','Mn','Mg','Ca','Ba' ...
    ,'Sr','K','Na','Ce','La','S6','S4','S1-','S2-','F','Cl','OH','C'});


end