%% Spinel structural formula
% last modified 01.08.2024

function [StrctFrm, apfu,options_definition]=spinel_Fe3unknown(data,headers,options)
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

cat=3.0; %cations per formula unit
Opfu=4.0; %oxygens per formula unit

%% Molecular weights

SiO2_mw=60.083;
TiO2_mw=79.865;
V2O3_mw=149.880;
Al2O3_mw=101.961;
Cr2O3_mw=151.989;
Fe2O3_mw=159.6874;
Y2O3_mw=225.809;
NiO_mw=74.692;
ZnO_mw=81.381;
CoO_mw=74.932;
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

%% Calculate cations units

MC=zeros(m,12);

%adds a column of zeros if SiO2 is not included in the calculation
if any(strcmp(headers,'SiO2'))
    MC(:,1)=data(:,strcmp(headers,'SiO2'))./SiO2_mw; %for SiO2
end

%adds a column of zeros if TiO2 is not included in the calculation
if any(strcmp(headers,'TiO2'))
    MC(:,2)=data(:,strcmp(headers,'TiO2'))./TiO2_mw; %for TiO2
end 

MC(:,3)=(data(:,strcmp(headers,'Al2O3'))./Al2O3_mw).*2; %for Al2O3

%adds a column of zeros if V2O3 is not included in the calculation
if any(strcmp(headers,'V2O3'))
    MC(:,4)=(data(:,strcmp(headers,'V2O3'))./V2O3_mw).*2; %for V2O3
end 

MC(:,5)=(data(:,strcmp(headers,'Cr2O3'))./Cr2O3_mw).*2; %for Cr2O3

%adds a column of zeros if NiO is not included in the calculation
if any(strcmp(headers,'NiO'))
    MC(:,6)=data(:,strcmp(headers,'NiO'))./NiO_mw; %for NiO
end

%adds a column of zeros if ZnO is not included in the calculation
if any(strcmp(headers,'ZnO'))
    MC(:,7)=data(:,strcmp(headers,'ZnO'))./ZnO_mw; %for ZnO
end

%adds a column of zeros if CoO is not included in the calculation
if any(strcmp(headers,'CoO'))
    MC(:,8)=data(:,strcmp(headers,'CoO'))./CoO_mw; %for CoO
end

%calculates for FeO and/or Fe2O3
if any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) %if FeO and Fe2O3 are both included as inputs

    %If FeO and Fe2O3 are both given, Fe2O3 is converted to FeO and
    %combined with FeO
    MC(:,9)=(data(:,strcmp(headers,'Fe2O3')).*((2.*FeO_mw)./Fe2O3_mw)+data(:,strcmp(headers,'FeO')))./FeO_mw;
else
    if any(strcmp(headers,'FeO')) %if FeO is FeO total

        MC(:,9)=data(:,strcmp(headers,'FeO'))./FeO_mw;

    else

        %calculate moles of cations
        MC(:,9)=(data(:,strcmp(headers,'Fe2O3')).*((2.*FeO_mw)./Fe2O3_mw))./FeO_mw;

    end
end

%calculates for MnO and/or Mn2O3
if any(strcmp(headers,'MnO')) && any(strcmp(headers,'Mn2O3')) %if MnO and Mn2O3 are both included as inputs

    %If MnO and Mn2O3 are both given, Mn2O3 is converted to MnO and
    %combined with MnO
    MC(:,10)=(data(:,strcmp(headers,'Mn2O3')).*((2.*MnO_mw)./Mn2O3_mw)+data(:,strcmp(headers,'MnO')))./MnO_mw;
else
    if any(strcmp(headers,'MnO')) %if MnO is MnO total

        MC(:,10)=data(:,strcmp(headers,'MnO'))./MnO_mw;

    elseif any(strcmp(headers,'Mn2O3')) %if Mn2O3 total is given, converted to MnO total

        %calculate moles of cations
        MC(:,10)=(data(:,strcmp(headers,'Mn2O3')).*((2.*MnO_mw)./Mn2O3_mw))./MnO_mw;

    else
        % if MnO total or Mn2O3 total is not included, Mn = 0
    end
end

MC(:,11)=data(:,strcmp(headers,'MgO'))./MgO_mw; %for MgO

%% Calculate normalized cations units

MCnormfact=cat./sum(MC,2); %normalization factor
MCnorm=MCnormfact.*MC; %creates a matrix of normalized cation

%% Calculate Oxygen Units

O2(:,1)=MCnorm(:,1).*2; %for SiO2
O2(:,2)=MCnorm(:,2).*2; %for TiO2
O2(:,3)=MCnorm(:,3).*(3/2); %for Al2O3
O2(:,4)=MCnorm(:,4).*(3/2); %for V2O3
O2(:,5)=MCnorm(:,5).*(3/2); %for Cr2O3
O2(:,6)=MCnorm(:,6); %for NiO
O2(:,7)=MCnorm(:,7); %for ZnO
O2(:,8)=MCnorm(:,8); %for CoO
O2(:,9)=MCnorm(:,9); %for FeO
O2(:,10)=MCnorm(:,10); %for MnO
O2(:,11)=MCnorm(:,11); %for MgO

O2total=sum(O2,2); %O2 totals

%% Atoms PFU

apfu(:,1)=MCnorm(:,1); %for Si
apfu(:,2)=MCnorm(:,2); %for Ti
apfu(:,3)=MCnorm(:,3); %for Al
apfu(:,4)=MCnorm(:,4); %for V
apfu(:,5)=MCnorm(:,5); %for Cr
apfu(:,7)=MCnorm(:,6); %for Ni
apfu(:,8)=MCnorm(:,7); %for Zn
apfu(:,9)=MCnorm(:,8); %for Co
apfu(:,11)=MCnorm(:,10); %for Mn
apfu(:,12)=MCnorm(:,11); %for Mg

%calculation of Fe3+ from stoichiometry and charge balance
%the following if statement firsts checks if totalO2 = 4
%if so, then there is no Fe3+
%if totalO2 < 4, then we assume that the deficiency is caused by the
%assumption Fetotal = Fe2+
%in the nested if statement, if FeTotal > 2*(4-totalO2) then the amount
%of Fe3+ = 2*(4-totalO2), if false then, all Fe is Fe3+

for c=1:m
    if (Opfu-O2total(c,1)) >= 0
        if MCnorm(c,9) > 2*(Opfu-O2total(c,1))
            apfu(c,6)=2*(Opfu-O2total(c,1)); 
        else
            apfu(c,6)=MCnorm(c,9);
        end
    else
        apfu(c,6)=0;
    end
end

apfu(:,10)=MCnorm(:,9)-apfu(:,6); %the apfu of Fe2+ equals totalFe-Fe3+

apfu(:,13)=sum(apfu,2); %calculations the total, which should be 3

% Oxygen deficiency 
apfu(:,14)=Opfu-O2total; %must be greater than zero

XMg = apfu(:,12)./(apfu(:,12)+apfu(:,10)); %Mg/(Mg + Fe2)
XFe = 1-XMg; 

%% Structural Formula

StrctFrm=zeros(m,15);

%D 
StrctFrm(:,1)=apfu(:,3); %Al (D)
StrctFrm(:,2)=apfu(:,4); %V (D)
StrctFrm(:,3)=apfu(:,5); %Cr (D)
StrctFrm(:,4)=apfu(:,6); %Fe3 (D)

D_left = 2 - StrctFrm(:,1:1:5); %space remaining after Al + V + Cr + Fe3
SiTi = apfu(:,1) + apfu(:,2); %Si + Ti on A

%Fe2+ (D)
for c=1:m
    if SiTi(c,1) > 0 %Is Si + Ti > 0? If y, then some Fe2+ is assigned to D
        if SiTi(c,1) >= (apfu(c,10) + apfu(c,12)) %For low Fe2+ spinel, Si + Ti may be > Fe2+
            StrctFrm(c,5)=apfu(c,10); %Fe2+ is assigned to D but maintaining equipartitioning
        else
            StrctFrm(c,5)= SiTi(c,1) * XFe(c,1); %Fe2 in D = (Si + Ti) * XFe
        end
    else
        StrctFrm(c,5)=0; %if Si + Ti < 0, then no Fe2+ goes into D
    end
end

%Mg (D)
for c=1:m
    if SiTi(c,1) > 0 %Is Si + Ti > 0? If y, then some Mg is assigned to D
        if SiTi(c,1) >= (apfu(c,10) + apfu(c,12)) %For low Mg spinel, Si + Ti may be > Mg
            StrctFrm(c,6)=apfu(c,12); %Mg is assigned to D but maintaining equipartitioning
        else
            StrctFrm(c,6)= SiTi(c,1) * XMg(c,1); %Mg in D = (Si + Ti) * XMg
        end
    else
        StrctFrm(c,6)=0; %if Si + Ti < 0, then no Fe2+ goes into D
    end
end

StrctFrm(:,7)=sum(StrctFrm(:,1:1:6),2); %Sum of D

%A 
StrctFrm(:,8)=apfu(:,1); %Si (A)
StrctFrm(:,9)=apfu(:,2); %Ti (A)
StrctFrm(:,10)=apfu(:,7); %Ni (A)
StrctFrm(:,11)=apfu(:,8); %Zn (A)
StrctFrm(:,12)=apfu(:,9); %Co (A)
StrctFrm(:,13)=apfu(:,10)-StrctFrm(:,5); %Fe2 (A)
StrctFrm(:,14)=apfu(:,11); %Mn2 (A)
StrctFrm(:,15)=apfu(:,12)-StrctFrm(:,6); %Mg (A)
StrctFrm(:,16)=sum(StrctFrm(:,8:1:15),2); %sum of A

%% Endmembers 

XAl = apfu(:,3)./(StrctFrm(:,7)); %Al / sum of D cations
XV = apfu(:,4)./(StrctFrm(:,7)); %V / sum of D cations
XCr = apfu(:,5)./(StrctFrm(:,7)); %Cr / sum of D cations
XFe3 = apfu(:,6)./(StrctFrm(:,7)); %Fe3 / sum of D cations
XUlvGrp = (StrctFrm(:,5)+StrctFrm(:,6))./(StrctFrm(:,7)); %Fe2 in D / sum of D cations

XSi = apfu(:,1)./(apfu(:,1)+apfu(:,2)); %Si / (Si + Ti)
XTi = 1-XSi; 

%A site fractions
XA(:,1) = StrctFrm(:,13)./sum(StrctFrm(:,10:1:15),2); %Fe2/(Fe2 + Mg + Mn2 + Ni + Zn + Co)
XA(:,2) = StrctFrm(:,15)./sum(StrctFrm(:,10:1:15),2); %Mg/(Fe2 + Mg + Mn2 + Ni + Zn + Co)
XA(:,3) = StrctFrm(:,14)./sum(StrctFrm(:,10:1:15),2); %Mn2/(Fe2 + Mg + Mn2 + Ni + Zn + Co)
XA(:,4) = StrctFrm(:,10)./sum(StrctFrm(:,10:1:15),2); %Ni/(Fe2 + Mg + Mn2 + Ni + Zn + Co)
XA(:,5) = StrctFrm(:,11)./sum(StrctFrm(:,10:1:15),2); %Zn/(Fe2 + Mg + Mn2 + Ni + Zn + Co)
XA(:,6) = StrctFrm(:,12)./sum(StrctFrm(:,10:1:15),2); %Co/(Fe2 + Mg + Mn2 + Ni + Zn + Co)

%endmembers
Endmembers(:,1:1:6) = XCr .* XA; %Cr endmembers 
Endmembers(:,7:1:12) = XFe3 .* XA; %Fe3 endmembers 
Endmembers(:,13:1:18) = XAl .* XA; %Al endmembers 
Endmembers(:,19:1:24) = XV .* XA; %V endmembers 

%Si Ulvöspinels
Endmembers(:,25) = XUlvGrp .* XSi .* XFe; %ahrensite (SiFe2O4)
Endmembers(:,26) = XUlvGrp .* XSi .* XMg; %Ringwoodite (SiMg2O4)

%Ti Ulvöspinels
Endmembers(:,27) = XUlvGrp .* XTi .* XFe; %Ulvöspinel (TiFe2O4)
Endmembers(:,28) = XUlvGrp .* XTi .* XMg; %Qandilite (TiMg2O4)

%% Outputs

%limit on significant digits (eliminates rounding noise)
StrctFrm(StrctFrm<1e-6) = 0;
Endmembers(Endmembers<1e-3) = 0; %limit on endmember noise (cannot be less than a fraction of a percent)
apfu(apfu<1e-6) = 0;
apfu(:,14)=Opfu-O2total; %also adds O2 def to the apfu output

all=[StrctFrm Endmembers apfu(:,14)];
all(isnan(all))=0; %remove NaN

StrctFrm=array2table(all,'VariableNames',{'Al_B','V3_B','Cr_B','Fe3_B','Fe2_B','Mg_B' ...
    'B_Sum','Si_A','Ti_A','Ni_A','Zn_A','Co_A','Fe2_A','Mn2_A','Mg_A','A_Sum', ...
    'Xchr','Xmchr','XMnchr','XNichr','XZnchr','XCochr','Xmag','Xmfr','Xjcb','Xtrv' ...
    'Xfrk','XComag','Xhc','Xspl','Xglx','Xchi','Xghn','XCospl','Xcou','Xmcou','Xvou' ...
    'XNicou','XZncou','XCocou','Xahr','Xrwd','Xqnd','Xuspl','O2_deficiency'});
apfu=array2table(apfu,'VariableNames',{'Si','Ti','Al','V3','Cr','Fe3','Ni','Zn','Co','Fe2','Mn2','Mg','Cation_Sum','O2_deficiency'});

end
            
