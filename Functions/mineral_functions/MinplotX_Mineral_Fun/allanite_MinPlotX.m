%allanite structural formula
% last modified 03.07.2024

function [StrctFrm, apfu,options_definition]=allanite_MinPlotX(data,headers,options)
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

Opfu=12.5; %oxygens per formula unit

%% Checks for Fe3+/Fetotal ratio

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

H2O_mw=18.015;
Na2O_mw=61.979;
MgO_mw=40.304;
Al2O3_mw=101.961;
SiO2_mw=60.083;
P2O5_mw=141.943;
K2O_mw=94.195;
CaO_mw=56.0774;
TiO2_mw=79.865;
Cr2O3_mw=151.989;
MnO_mw=70.937;
Mn2O3_mw=157.873;
Fe2O3_mw=159.6874;
FeO_mw=71.8442;
NiO_mw=74.692;
ZnO_mw=81.381;
SrO_mw=103.619;
Y2O3_mw=225.809;
ZrO2_mw=123.222;
BaO_mw=153.329;
La2O3_mw=325.808;
Ce2O3_mw=328.229;
Pr2O3_mw=329.812;
Nd2O3_mw=336.481;
Sm2O3_mw=348.72;
Eu2O3_mw=351.925;
Gd2O3_mw=362.50;
Tb2O3_mw=365.848;
Dy2O3_mw=373.0;
Ho2O3_mw=377.858;
Er2O3_mw=382.515;
Tm2O3_mw=385.865;
Yb2O3_mw=394.087;
Lu2O3_mw=397.931;
PbO_mw=223.2;
ThO2_mw=264.036;
UO2_mw=270.027;
F_mw=18.998;
Cl_mw=35.45;

%% Calculate moles of cations

MC=zeros(m,36);%create columns of zeros if optional data are not included

MC(:,1)=data(:,strcmp(headers,'SiO2'))./SiO2_mw; % for SiO2

%calculates for TiO2 if it is included in the analysis
if any(strcmp(headers,'TiO2'))
    MC(:,2)=data(:,strcmp(headers,'TiO2'))./TiO2_mw; %for TiO2
end

%calculates for ZrO2 if it is included in the analysis
if any(strcmp(headers,'ZrO2'))
    MC(:,3)=data(:,strcmp(headers,'ZrO2'))./ZrO2_mw; %for ZrO2
end

MC(:,4)=(data(:,strcmp(headers,'Al2O3'))./Al2O3_mw).*2; % for Al2O3

%calculates for Cr2O3 if it is included in the analysis
if any(strcmp(headers,'Cr2O3'))
    MC(:,5)=(data(:,strcmp(headers,'Cr2O3'))./Cr2O3_mw).*2; %for Cr2O3
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

if any(strcmp(headers,'MgO'))
    MC(:,10)=data(:,strcmp(headers,'MgO'))./MgO_mw; % for MgO
end

MC(:,11)=data(:,strcmp(headers,'CaO'))./CaO_mw; % for CaO

if any(strcmp(headers,'SrO'))
    MC(:,12)=data(:,strcmp(headers,'SrO'))./SrO_mw; % for SrO
end

%calculates for Na2O if it is included in the analysis
if any(strcmp(headers,'Na2O'))
    MC(:,13)=(data(:,strcmp(headers,'Na2O'))./Na2O_mw).*2; %for Na2O
end

%calculates for K2O if it is included in the analysis
if any(strcmp(headers,'K2O'))
    MC(:,14)=(data(:,strcmp(headers,'K2O'))./K2O_mw).*2; %for K2O
end

%calculates for La2O3 if it is included in the analysis
if any(strcmp(headers,'La2O3'))
    MC(:,15)=(data(:,strcmp(headers,'La2O3'))./La2O3_mw).*2; %for La2O3
end

%calculates for Ce2O3 if it is included in the analysis
if any(strcmp(headers,'Ce2O3'))
    MC(:,16)=(data(:,strcmp(headers,'Ce2O3'))./Ce2O3_mw).*2; %for Ce2O3
end

%calculates for Pr2O3 if it is included in the analysis
if any(strcmp(headers,'Pr2O3'))
    MC(:,17)=(data(:,strcmp(headers,'Pr2O3'))./Pr2O3_mw).*2; %for Pr2O3
end

%calculates for Nd2O3 if it is included in the analysis
if any(strcmp(headers,'Nd2O3'))
    MC(:,18)=(data(:,strcmp(headers,'Nd2O3'))./Nd2O3_mw).*2; %for Nd2O3
end

%calculates for Sm2O3 if it is included in the analysis
if any(strcmp(headers,'Sm2O3'))
    MC(:,19)=(data(:,strcmp(headers,'Sm2O3'))./Sm2O3_mw).*2; %for Sm2O3
end

%calculates for Eu2O3 if it is included in the analysis
if any(strcmp(headers,'Eu2O3'))
    MC(:,20)=(data(:,strcmp(headers,'Eu2O3'))./Eu2O3_mw).*2; %for Eu2O3
end

%calculates for Gd2O3 if it is included in the analysis
if any(strcmp(headers,'Gd2O3'))
    MC(:,21)=(data(:,strcmp(headers,'Gd2O3'))./Gd2O3_mw).*2; %for Gd2O3
end

%calculates for Tb2O3 if it is included in the analysis
if any(strcmp(headers,'Tb2O3'))
    MC(:,22)=(data(:,strcmp(headers,'Tb2O3'))./Tb2O3_mw).*2; %for Tb2O3
end

%calculates for Dy2O3 if it is included in the analysis
if any(strcmp(headers,'Dy2O3'))
    MC(:,23)=(data(:,strcmp(headers,'Dy2O3'))./Dy2O3_mw).*2; %for Dy2O3
end

%calculates for Ho2O3 if it is included in the analysis
if any(strcmp(headers,'Ho2O3'))
    MC(:,24)=(data(:,strcmp(headers,'Ho2O3'))./Ho2O3_mw).*2; %for Ho2O3
end

%calculates for Er2O3 if it is included in the analysis
if any(strcmp(headers,'Er2O3'))
    MC(:,25)=(data(:,strcmp(headers,'Er2O3'))./Er2O3_mw).*2; %for Er2O3
end

%calculates for Tm2O3 if it is included in the analysis
if any(strcmp(headers,'Tm2O3'))
    MC(:,26)=(data(:,strcmp(headers,'Tm2O3'))./Tm2O3_mw).*2; %for Tm2O3
end

%calculates for Yb2O3 if it is included in the analysis
if any(strcmp(headers,'Yb2O3'))
    MC(:,27)=(data(:,strcmp(headers,'Yb2O3'))./Yb2O3_mw).*2; %for Yb2O3
end

%calculates for Lu2O3 if it is included in the analysis
if any(strcmp(headers,'Lu2O3'))
    MC(:,28)=(data(:,strcmp(headers,'Lu2O3'))./Lu2O3_mw).*2; %for Lu2O3
end

%calculates for Y2O3 if it is included in the analysis
if any(strcmp(headers,'Y2O3'))
    MC(:,29)=(data(:,strcmp(headers,'Y2O3'))./Y2O3_mw).*2; %for Y2O3
end

%calculates for PbO if it is included in the analysis
if any(strcmp(headers,'PbO'))
    MC(:,30)=(data(:,strcmp(headers,'PbO'))./PbO_mw); %for PbO
end

%calculates for UO2 if it is included in the analysis
if any(strcmp(headers,'UO2'))
    MC(:,31)=(data(:,strcmp(headers,'UO2'))./UO2_mw); %for UO2
end

%calculates for ThO2 if it is included in the analysis
if any(strcmp(headers,'ThO2'))
    MC(:,32)=(data(:,strcmp(headers,'ThO2'))./ThO2_mw); %for ThO2
end

%calculates for P2O5 if it is included in the analysis
if any(strcmp(headers,'P2O5'))
    MC(:,33)=(data(:,strcmp(headers,'P2O5'))./P2O5_mw).*2; %for P2O5
end

%calculates for F if it is included in the analysis
if any(strcmp(headers,'F'))
    MC(:,34)=data(:,strcmp(headers,'F'))./F_mw; %for F
end

%calculates for Cl if it is included in the analysis
if any(strcmp(headers,'Cl'))
    MC(:,35)=data(:,strcmp(headers,'Cl'))./Cl_mw; %for Cl
end

%calculates for H2O if it is included in the analysis
if any(strcmp(headers,'H2O'))
    MC(:,36)=(data(:,strcmp(headers,'H2O'))./H2O_mw)*2; %for H2O
end

%% Calculate Oxygen Units

O2(:,1)=MC(:,1).*2; %for SiO2
O2(:,2)=MC(:,2).*2; %for TiO2
O2(:,3)=MC(:,3).*2; %for ZrO2
O2(:,4)=MC(:,4).*(3/2); %for Al2O3
O2(:,5)=MC(:,5).*(3/2); %for Cr2O3
O2(:,6)=MC(:,6).*(3/2); %for Fe2O3
O2(:,7)=MC(:,7).*(3/2); %for Mn2O3
O2(:,8)=MC(:,8); %for FeO
O2(:,9)=MC(:,9); %for MnO
O2(:,10)=MC(:,10); %for MgO
O2(:,11)=MC(:,11); %for CaO
O2(:,12)=MC(:,12); %for SrO
O2(:,13)=MC(:,13).*(1/2); %for Na2O
O2(:,14)=MC(:,14).*(1/2); %for K2O
O2(:,15)=MC(:,15).*(3/2); %for La2O3
O2(:,16)=MC(:,16).*(3/2); %for Ce2O3
O2(:,17)=MC(:,17).*(3/2); %for Pr2O3
O2(:,18)=MC(:,18).*(3/2); %for Nd2O3
O2(:,19)=MC(:,19).*(3/2); %for Sm2O3
O2(:,20)=MC(:,20).*(3/2); %for Eu2O3
O2(:,21)=MC(:,21).*(3/2); %for Gd2O3
O2(:,22)=MC(:,22).*(3/2); %for Tb2O3
O2(:,23)=MC(:,23).*(3/2); %for Dy2O3
O2(:,24)=MC(:,24).*(3/2); %for Ho2O3
O2(:,25)=MC(:,25).*(3/2); %for Er2O3
O2(:,26)=MC(:,26).*(3/2); %for Tm2O3
O2(:,27)=MC(:,27).*(3/2); %for Yb2O3
O2(:,28)=MC(:,28).*(3/2); %for Lu2O3
O2(:,29)=MC(:,29).*(3/2); %for Y2O3
O2(:,30)=MC(:,30); %for PbO
O2(:,31)=MC(:,31).*2; %for UO2
O2(:,32)=MC(:,32).*2; %for ThO2
O2(:,33)=MC(:,33).*(5/2); %for ThO2
O2(:,34)=MC(:,34); %F
O2(:,35)=MC(:,35); %Cl
O2(:,36)=MC(:,36).*(1/2); %OH

O2_Sum=sum(O2(:,1:36),2)-0.5.*(O2(:,34)+O2(:,35)); %sum of O2, including F and Cl

O2_N=(13)./O2_Sum; %normalization factor

%normalized moles of anions
N_Ox=O2.*O2_N;

if not(any(strcmp(headers,'H2O')))

    %initial oxygens of OH
    Hin=1-(N_Ox(:,34)+N_Ox(:,35)); %H = 1 - (F+Cl)
    O_OH=(0.5.*Hin)./O2_N; %oxygen moles of H2O
    H2Oi=O_OH.*H2O_mw; %initial H2O wt %

    for z=1:50

        O2(:,1)=MC(:,1).*2; %for SiO2
        O2(:,2)=MC(:,2).*2; %for TiO2
        O2(:,3)=MC(:,3).*2; %for ZrO2
        O2(:,4)=MC(:,4).*(3/2); %for Al2O3
        O2(:,5)=MC(:,5).*(3/2); %for Cr2O3
        O2(:,6)=MC(:,6).*(3/2); %for Fe2O3
        O2(:,7)=MC(:,7).*(3/2); %for Mn2O3
        O2(:,8)=MC(:,8); %for FeO
        O2(:,9)=MC(:,9); %for MnO
        O2(:,10)=MC(:,10); %for MgO
        O2(:,11)=MC(:,11); %for CaO
        O2(:,12)=MC(:,12); %for SrO
        O2(:,13)=MC(:,13).*(1/2); %for Na2O
        O2(:,14)=MC(:,14).*(1/2); %for K2O
        O2(:,15)=MC(:,15).*(3/2); %for La2O3
        O2(:,16)=MC(:,16).*(3/2); %for Ce2O3
        O2(:,17)=MC(:,17).*(3/2); %for Pr2O3
        O2(:,18)=MC(:,18).*(3/2); %for Nd2O3
        O2(:,19)=MC(:,19).*(3/2); %for Sm2O3
        O2(:,20)=MC(:,20).*(3/2); %for Eu2O3
        O2(:,21)=MC(:,21).*(3/2); %for Gd2O3
        O2(:,22)=MC(:,22).*(3/2); %for Tb2O3
        O2(:,23)=MC(:,23).*(3/2); %for Dy2O3
        O2(:,24)=MC(:,24).*(3/2); %for Ho2O3
        O2(:,25)=MC(:,25).*(3/2); %for Er2O3
        O2(:,26)=MC(:,26).*(3/2); %for Tm2O3
        O2(:,27)=MC(:,27).*(3/2); %for Yb2O3
        O2(:,28)=MC(:,28).*(3/2); %for Lu2O3
        O2(:,29)=MC(:,29).*(3/2); %for Y2O3
        O2(:,30)=MC(:,30); %for PbO
        O2(:,31)=MC(:,31).*2; %for UO2
        O2(:,32)=MC(:,32).*2; %for ThO2
        O2(:,33)=MC(:,33).*(5/2); %for ThO2
        O2(:,34)=MC(:,34); %F
        O2(:,35)=MC(:,35); %Cl
        O2(:,36)=O_OH; %H2O

        O2_Sum=sum(O2(:,1:36),2)-0.5.*(O2(:,13)+O2(:,14)); %sum of O2, including F and Cl

        O2_N=(13)./O2_Sum; %normalization factor

        %normalized moles of anions
        N_Ox=O2.*O2_N;

        Hin=1-(N_Ox(:,34)+N_Ox(:,35)); %H = 1 - (F+Cl)
        O_OH=(0.5.*Hin)./O2_N; %oxygen moles of H2O
        H2Oi=O_OH.*H2O_mw; %initial H2O wt %

    end
end

%% atoms pfu
apfu(:,1)=N_Ox(:,1)./2; %Si
apfu(:,2)=N_Ox(:,2)./2; %Ti
apfu(:,3)=N_Ox(:,3)./2; %Zr
apfu(:,4)=N_Ox(:,4).*(2/3); %Al
apfu(:,5)=N_Ox(:,5).*(2/3); %Cr
apfu(:,6)=N_Ox(:,6).*(2/3); %Fe3
apfu(:,7)=N_Ox(:,7).*(2/3); %Mn3
apfu(:,8)=N_Ox(:,8); %Fe2
apfu(:,9)=N_Ox(:,9); %Mn2
apfu(:,10)=N_Ox(:,10); %Mg
apfu(:,11)=N_Ox(:,11); %Ca
apfu(:,12)=N_Ox(:,12); %Sr
apfu(:,13)=N_Ox(:,13).*2; %Na
apfu(:,14)=N_Ox(:,14).*2; %K
apfu(:,15)=N_Ox(:,15).*(2/3); %La
apfu(:,16)=N_Ox(:,16).*(2/3); %Ce
apfu(:,17)=N_Ox(:,17).*(2/3); %Pr
apfu(:,18)=N_Ox(:,18).*(2/3); %Nd
apfu(:,19)=N_Ox(:,19).*(2/3); %Sm
apfu(:,20)=N_Ox(:,20).*(2/3); %Eu
apfu(:,21)=N_Ox(:,21).*(2/3); %Gd
apfu(:,22)=N_Ox(:,22).*(2/3); %Tb
apfu(:,23)=N_Ox(:,23).*(2/3); %Dy
apfu(:,24)=N_Ox(:,24).*(2/3); %Ho
apfu(:,25)=N_Ox(:,25).*(2/3); %Er
apfu(:,26)=N_Ox(:,26).*(2/3); %Tm
apfu(:,27)=N_Ox(:,27).*(2/3); %Yb
apfu(:,28)=N_Ox(:,28).*(2/3); %Lu
apfu(:,29)=N_Ox(:,29).*(2/3); %Y
apfu(:,30)=N_Ox(:,30); %Pb
apfu(:,31)=N_Ox(:,31)./2; %U4
apfu(:,32)=N_Ox(:,32)./2; %Th4
apfu(:,33)=N_Ox(:,33).*(2/5); %P
apfu(:,34)=N_Ox(:,34); %F
apfu(:,35)=N_Ox(:,35); %Cl
apfu(:,36)=N_Ox(:,36).*2; %OH
apfu(:,37)=sum(apfu(:,1:1:33),2); %calculations the total

%% Structural Formula
StrctFrm=zeros(m,30);

%T SITE

StrctFrm(:,1)=apfu(:,1); %Si
StrctFrm(:,2)=apfu(:,33); %P
StrctFrm(:,3)=apfu(:,3); %Zr

%Al(T)
for c=1:m
    if 3-sum(StrctFrm(c,1:1:3),2)>0 %Is 3-Si-P-Zr > 0? If y, then some Al goes into T
        StrctFrm(c,4)=3-sum(StrctFrm(c,1:1:3),2);
    else
        StrctFrm(c,4)=0; %if Si=3, then no Al goes into T
    end
end

StrctFrm(:,5)=sum(StrctFrm(:,1:1:4),2); %T sum

%M Site

StrctFrm(:,6)=apfu(:,4)-StrctFrm(:,4); % Al(M) = Altotal - Al(T)
StrctFrm(:,7)=apfu(:,2); %Ti (M)
StrctFrm(:,8)=apfu(:,5); %Cr (M)
StrctFrm(:,9)=apfu(:,6); %Fe3 (M)
StrctFrm(:,10)=apfu(:,7); %Mn3 (M)


%setting up equipartitioning of Fe2, Mg, and Mn2 between A and M sites
M_rem = 3 - sum(StrctFrm(:,6:1:10),2); % remaining space on M after trivalent cations

XMn=apfu(:,9)./sum(apfu(:,8:1:10),2); %Mn / (Mg + Fe2 + Mn2)
XMg=apfu(:,10)./sum(apfu(:,8:1:10),2); %Mg / (Mg + Fe2 + Mn2)
XFe=apfu(:,8)./sum(apfu(:,8:1:10),2); %Fe2 / (Mg + Fe2 + Mn2)

%Remove division by zero
XMn(isnan(XMn))=0;
XMg(isnan(XMg))=0;
XFe(isnan(XFe))=0;

%assigning divalent cations
for c=1:m
    if M_rem(c,1) > 0
        if M_rem(c,1) >= sum(apfu(c,8:1:10),2)
            %if the remaining space on M is ≥ total Fe2, Mn2, and Mg
            %all Fe2, Mn2, and Mg are assigned to M
            StrctFrm(c,11)=apfu(c,8); %Fe2(M)
            StrctFrm(c,12)=apfu(c,9); %Mn2(M)
            StrctFrm(c,13)=apfu(c,10); %Mg(M)
        else
            %if there is extra Mg, Fe2, and Mn2, these elements will be
            %equipartitioned between M and A
            StrctFrm(c,11)=M_rem(c,1).*XFe(c,1); %Fe2(M)
            StrctFrm(c,12)=M_rem(c,1).*XMn(c,1); %Mn2(M)
            StrctFrm(c,13)=M_rem(c,1).*XMg(c,1); %Mg(M)
        end
    else
        %if 3+ and 4+ cations fill the M site, no Fe2, Mn2, and Mg are assigned to
        %M
    end
end

StrctFrm(:,14)=sum(StrctFrm(:,6:1:13),2); %M sum

%A site
StrctFrm(:,15)=apfu(:,8)-StrctFrm(:,11); %Fe2(A) = Fe2total - Fe2(M)
StrctFrm(:,16)=apfu(:,9)-StrctFrm(:,12); %Mn2(A) = Mn2total - Mn2(M)
StrctFrm(:,17)=apfu(:,10)-StrctFrm(:,13); %Mg(A) = Mgtotal - Mg(M)
StrctFrm(:,18)=apfu(:,11); %Ca(A)
StrctFrm(:,19)=apfu(:,12); %Sr(A)
StrctFrm(:,20)=apfu(:,13); %Na(A)
StrctFrm(:,21)=apfu(:,14); %K(A)
StrctFrm(:,22)=sum(apfu(:,15:1:29),2); %REE + Y (A)
StrctFrm(:,23)=apfu(:,30); %Pb(A)
StrctFrm(:,24)=apfu(:,31); %U4(A)
StrctFrm(:,25)=apfu(:,32); %Th(A)

StrctFrm(:,26)=sum(StrctFrm(:,15:1:25),2); %A sum

%halogens
StrctFrm(:,27)=apfu(:,34); %F (H)
StrctFrm(:,28)=apfu(:,35); %Cl (H)
StrctFrm(:,29)=apfu(:,36); %OH (H)
StrctFrm(:,30)=sum(StrctFrm(:,27:1:29),2); %H sum

%limit on significant digits (eliminates rounding noise)
StrctFrm(StrctFrm<1e-6) = 0;
apfu(apfu<1e-6) = 0;

StrctFrm=array2table(StrctFrm,'VariableNames',{'Si_T','P_T','Zr_T','Al_T','Sum_T', ...
    'Al_M','Ti_M','Cr_M','Fe3_M','Mn3_M','Fe2_M','Mn2_M','Mg_M','Sum_M','Fe2_A','Mn2_A', ...
    'Mg_A','Ca_A','Sr_A','Na_A','K_A','YREE_A','Pb_A','U_A','Th_A','Sum_A','F_H','Cl_H','OH_H','Sum_H'});
apfu=array2table(apfu,'VariableNames',{'Si','Ti','Zr','Al','Cr','Fe3','Mn3','Fe2','Mn2', ...
    'Mg','Ca','Sr','Na','K','La','Ce3','Pr','Nd','Sm','Eu3','Gd','Tb','Dy', ...
    'Ho','Er','Tm','Yb','Lu','Y','Pb','U4','Th','P','F','Cl','OH','Cation_Sum'});















