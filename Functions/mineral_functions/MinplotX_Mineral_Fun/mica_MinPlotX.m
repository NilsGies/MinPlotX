%mica structural formula
%last modified 01.08.2024

function [StrctFrm, apfu,options_definition]=mica_MinPlotX(data,headers,options)
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

%% start calculation
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

%% Moles of oxides

MC=zeros(m,15);%create columns of zeros if optional data are not included

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

%calculates for MnO and/or Mn2O3
if any(strcmp(headers,'MnO')) && any(strcmp(headers,'Mn2O3')) %if MnO and Mn2O3 are both included as inputs 

     %If MnO and Mn2O3 are both given, Mn2O3 is converted to MnO and
     %combined with MnO
     MC(:,7)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw)+data(:,strcmp(headers,'MnO')))./MnO_mw;
else
    if any(strcmp(headers,'MnO')) %if MnO is MnO total

        MC(:,7)=data(:,strcmp(headers,'MnO'))./MnO_mw;

    elseif any(strcmp(headers,'Mn2O3'))%if Mn2O3 total is given, converted to MnO total

        %calculate moles of cations
        MC(:,7)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw))./MnO_mw;
    else
        %Mn2O3 and MnO are not included, Mn2 = 0
    end
end

MC(:,8)=data(:,strcmp(headers,'MgO'))./MgO_mw; %for MgO

%calculates for CaO if it is included in the analysis
if any(strcmp(headers,'CaO'))
    MC(:,9)=data(:,strcmp(headers,'CaO'))./CaO_mw; %for CaO
end

MC(:,10)=(data(:,strcmp(headers,'Na2O'))./Na2O_mw).*2; %for Na2O
MC(:,11)=(data(:,strcmp(headers,'K2O'))./K2O_mw).*2; %for K2O

%calculates for BaO if it is included in the analysis
if any(strcmp(headers,'BaO'))
    MC(:,12)=data(:,strcmp(headers,'BaO'))./BaO_mw; %for BaO
end

%calculates for F if it is included in the analysis
if any(strcmp(headers,'F'))
    MC(:,13)=data(:,strcmp(headers,'F'))./F_mw; %for F
end

%calculates for Cl if it is included in the analysis
if any(strcmp(headers,'Cl'))
    MC(:,14)=data(:,strcmp(headers,'Cl'))./Cl_mw; %for Cl
end

%calculates for H2O if it is included in the analysis
if any(strcmp(headers,'H2O'))
    MC(:,15)=(data(:,strcmp(headers,'H2O'))./H2O_mw).*2; %for H2O
end

%% Moles of O2 units
O2(:,1)=MC(:,1).*2; %SiO2
O2(:,2)=MC(:,2).*2; %TiO2
O2(:,3)=MC(:,3).*(3/2); %Al2O3
O2(:,4)=MC(:,4).*(3/2); %Cr2O3
O2(:,5)=MC(:,5).*(3/2); %Fe2O3
O2(:,6)=MC(:,6); %FeO
O2(:,7)=MC(:,7); %MnO
O2(:,8)=MC(:,8); %MgO
O2(:,9)=MC(:,9); %CaO
O2(:,10)=MC(:,10).*(1/2); %Na2O
O2(:,11)=MC(:,11).*(1/2); %K2O
O2(:,12)=MC(:,12); %BaO
O2(:,13)=MC(:,13); %F
O2(:,14)=MC(:,14); %Cl
O2(:,15)=MC(:,15).*(1/2); %H2O

O2_Sum=sum(O2(:,1:15),2)-0.5.*(O2(:,13)+O2(:,14)); %sum of O2, including F and Cl

O2_N=(12)./O2_Sum; %normalization factor

%normalized moles of anions
N_Ox=O2.*O2_N;

%% Iterate OH

if not(any(strcmp(headers,'H2O')))

    %initial oxygens of OH
    Hin=2-(N_Ox(:,13)+N_Ox(:,14)); %H = 2 - (F+Cl)
    O_OH=(0.5.*Hin)./O2_N; %oxygen moles of H2O

    for z=1:50
        O2(:,1)=MC(:,1).*2; %SiO2
        O2(:,2)=MC(:,2).*2; %TiO2
        O2(:,3)=MC(:,3).*(3/2); %Al2O3
        O2(:,4)=MC(:,4).*(3/2); %Cr2O3
        O2(:,5)=MC(:,5).*(3/2); %Fe2O3
        O2(:,6)=MC(:,6); %FeO
        O2(:,7)=MC(:,7); %MnO
        O2(:,8)=MC(:,8); %MgO
        O2(:,9)=MC(:,9); %CaO
        O2(:,10)=MC(:,10).*(1/2); %Na2O
        O2(:,11)=MC(:,11).*(1/2); %K2O
        O2(:,12)=MC(:,12); %BaO
        O2(:,13)=MC(:,13); %F
        O2(:,14)=MC(:,14); %Cl
        O2(:,15)=O_OH; %H2O

        O2_Sum=sum(O2(:,1:15),2)-0.5.*(O2(:,13)+O2(:,14)); %sum of O2, including F and Cl

        O2_N=(12)./O2_Sum; %normalization factor

        %normalized moles of anions
        N_Ox=O2.*O2_N;

        Hin=2-(N_Ox(:,13)+N_Ox(:,14)); %OH = 2 - (F+Cl)
        O_OH=(0.5.*Hin)./O2_N; %oxygen moles of H2O
    end
end

%% atoms pfu

apfu(:,1)=N_Ox(:,1)./2; %Si
apfu(:,2)=N_Ox(:,2)./2; %Ti
apfu(:,3)=N_Ox(:,3).*(2/3); %Al
apfu(:,4)=N_Ox(:,4).*(2/3); %Cr
apfu(:,5)=N_Ox(:,5).*(2/3); %Fe3+
apfu(:,6)=N_Ox(:,6); %Fe2+
apfu(:,7)=N_Ox(:,7); %Mn
apfu(:,8)=N_Ox(:,8); %Mg
apfu(:,9)=N_Ox(:,9); %Ca
apfu(:,10)=N_Ox(:,10).*2; %Na
apfu(:,11)=N_Ox(:,11).*2; %K
apfu(:,12)=N_Ox(:,12); %Ba
apfu(:,13)=N_Ox(:,13); %F
apfu(:,14)=N_Ox(:,14); %Cl
apfu(:,15)=N_Ox(:,15).*2; %OH
apfu(:,16)=sum(apfu(:,1:1:12),2); %calculations the total

%% Structural Formula

StrctFrm=zeros(m,20);

%T site
StrctFrm(:,1)=apfu(:,1); %Si (T)

StrctFrm(:,3)=apfu(:,5).*tetra_Fe3; %Fe3+ (T)

%Al (T)
for c=1:m
    if 4-StrctFrm(c,1)-StrctFrm(c,3) > apfu(c,3)
        StrctFrm(c,2)=apfu(c,3);
    else
        StrctFrm(c,2)=4-StrctFrm(c,1)-StrctFrm(c,3);
    end
end

StrctFrm(:,4)=StrctFrm(:,1)+StrctFrm(:,2)+StrctFrm(:,3); %Sum of T

StrctFrm(:,5)=apfu(:,3)-StrctFrm(:,2); %Altotal-Al(T)=Al(M)
StrctFrm(:,6)=apfu(:,2); %Ti (M)
StrctFrm(:,7)=apfu(:,4); %Cr (M)
StrctFrm(:,8)=apfu(:,5)-StrctFrm(:,3); %Fe3+ (M)
StrctFrm(:,9)=apfu(:,6); %Fe2 (M)
StrctFrm(:,10)=apfu(:,7); %Mn (M)
StrctFrm(:,11)=apfu(:,8); %Mg (M)
StrctFrm(:,12)=sum(StrctFrm(:,5:1:11),2); %sum (M)

StrctFrm(:,13)=apfu(:,9); %Ca (I)
StrctFrm(:,14)=apfu(:,10); %Na (I)
StrctFrm(:,15)=apfu(:,11); %K (I)
StrctFrm(:,16)=apfu(:,12); %Ba (I)

StrctFrm(:,17)=sum(StrctFrm(:,13:1:16),2)+StrctFrm(:,12)+StrctFrm(:,4); %cation sum

StrctFrm(:,18)=apfu(:,13); % F (A)
StrctFrm(:,19)=apfu(:,14); % Cl (A)
StrctFrm(:,20)=apfu(:,15); % OH (A)

%% Endmembers

%Dioctahedral Micas
%the sum of M=3 cations in TriOct micas, and M=2 cations in DiOct,

for c=1:m
    if StrctFrm(c,12)>2 %If M>2, there is a TriOct component
        XTriOct(c,:)=StrctFrm(c,12)-2; %Sum of M-2, scales from 0 to 1
    else
        if StrctFrm(c,12)>3
            XTriOct(c,:)=1;  %If M>3, there is a no DiOct Component
        else
            XTriOct(c,:)=0; %If M<3 and M<2 then there is no TriOct component
        end
    end
end

XDiOct=1-XTriOct; %fraction of dioctahedral mica

%calculate fractionn of Ms + Pg + Mrg + Pyl
%looks at octahedrally coordinated Al
for c=1:m
    if StrctFrm(c,5)>2
        XMsum(c,:)=1; %If viAl > 2, fix it as 1 (upper limit for muscovite, paragonite, margarite, pyrophyllite)
    else
        if StrctFrm(c,5) + StrctFrm(c,8) <1
            XMsum(c,:)=0;  %If viAl + viFe3 < 1
        else
            XMsum(c,:)=(StrctFrm(c,5)+StrctFrm(c,8))-1; %scales viAl between 0 and 1
        end
    end
end

%Calculate XMg
for c=1:m
    if StrctFrm(c,11)<=0
        XMg(c,:)=0; %Avoids division by 0 if no Mg is present
    else
        XMg(c,:)=StrctFrm(c,11)./(StrctFrm(c,11)+StrctFrm(c,9)); %calculates XMg
    end
end

%total Celadonite component
XCelTot=1-XMsum; %total amout (Fe Al-Celadonite + Mg Al-Celadonite)

%Fe3-Celadonite component 
for c=1:m
    if StrctFrm(c,8)+StrctFrm(c,5) > 0
        XFe3Cel(c,1)=StrctFrm(c,8)./(StrctFrm(c,8)+StrctFrm(c,5)); %Fe3/(Fe3 + viAl) on M 
    else 
        XFe3Cel(c,1)=0; %corrects for NaN if no viAl or viFe3
    end
end

XCel=XFe3Cel.*XCelTot; %fraction of Al-Celadonite (Mg-bearing endmember)
XAlCel=XCelTot-XCel; %fraction ofl-Celadonite

%fraction Ms, Pg, and Mrg
XMPM=(StrctFrm(:,13)+StrctFrm(:,14)+StrctFrm(:,15)).*XMsum; %Sum of Ca + Na + K multiplied by the total fraction of Ms, Pg, Mrg, and Prl

XPrl=XMsum-XMPM; %fraction of pyrophyllite

%Calculate XMrg: XCa * total proportion of Ms + Pg + Mrg
for c=1:m
    if (StrctFrm(c,13)+StrctFrm(c,14)+StrctFrm(c,15))<=0
        XMrg(c,:)=0; %If there is no Ca, Na, and K, then no XCa, avoids division by 0
    else
        XMrg(c,:)=(StrctFrm(c,13)./(StrctFrm(c,13)+StrctFrm(c,14)+StrctFrm(c,15))).*XMPM(c,:);
    end
end

%Calculate XPg: XNa * total proportion of Ms + Pg + Mrg
for c=1:m
    if (StrctFrm(c,13)+StrctFrm(c,14)+StrctFrm(c,15))<=0
        XPg(c,:)=0; %If there is no Ca, Na, and K, then no XNa, avoids division by 0
    else
        XPg(c,:)=(StrctFrm(c,14)./(StrctFrm(c,13)+StrctFrm(c,14)+StrctFrm(c,15))).*XMPM(c,:);
    end
end

%Calculate XMs: XK * total proportion of Ms + Pg + Mrg
for c=1:m
    if (StrctFrm(c,13)+StrctFrm(c,14)+StrctFrm(c,15))<=0
        XMs(c,:)=0; %If there is no Ca, Na, and K, then no XNa, avoids division by 0
    else
        XMs(c,:)=(StrctFrm(c,15)./(StrctFrm(c,13)+StrctFrm(c,14)+StrctFrm(c,15))).*XMPM(c,:);
    end
end

DiOct_Endmembers(:,1)=XCel.*XDiOct; %Celadonite
DiOct_Endmembers(:,2)=XAlCel.*XDiOct; %Al-Celadonite
DiOct_Endmembers(:,3)=XPrl.*XDiOct; %Pyrophyllite
DiOct_Endmembers(:,4)=XMrg.*XDiOct; %Margarite
DiOct_Endmembers(:,5)=XPg.*XDiOct; %Paragonite
DiOct_Endmembers(:,6)=XMs.*XDiOct; %Muscovite
DiOct_Endmembers(:,7)=XTriOct; %trioctahedral component

%Trioctahedral Micas

%pholgopite-annite solid solution
XPhlAnn=apfu(:,1)-2; %fraction of the Phl+Ann endmembers
%Phl and Ann has 3 Si, whereas Sid and East has 2
XPhl=XPhlAnn.*XMg; %unnormalized fraction of Phl
XAnn=XPhlAnn-XPhl; %unnormalized fraction of Ann

%siderophyllite-eastonite solid solution
XSidEast=1-XPhlAnn; %fraction of the Sid+Eastendmembers
XEast=XSidEast.*XMg;
XSid=XSidEast-XEast;

TriOct_Endmembers(:,1)=XPhl.*XTriOct; %phlogopite
TriOct_Endmembers(:,2)=XAnn.*XTriOct; %annite
TriOct_Endmembers(:,3)=XEast.*XTriOct; %eastonite
TriOct_Endmembers(:,4)=XSid.*XTriOct; %siderophyllite
TriOct_Endmembers(:,5)=XDiOct; %dioctohedral component

TriOct_Endmembers(TriOct_Endmembers<1e-3) = 0; %limit on endmember noise (cannot be less than a fraction of a percent)
DiOct_Endmembers(DiOct_Endmembers<1e-3) = 0; %limit on endmember noise (cannot be less than a fraction of a percent)

%% Outputs


all=[StrctFrm DiOct_Endmembers TriOct_Endmembers];

%limit on significant digits (eliminates rounding noise)
all(all<1e-6) = 0;
apfu(apfu<1e-6) = 0;
all(isnan(all))=0; %remove NaN

StrctFrm=array2table(all,'VariableNames',{'Si_T','Al_T','Fe3_T','Sum_T','Al_M','Ti_M','Cr_M','Fe3_M','Fe2_M','Mn_M','Mg_M','M_Sum','Ca_I','Na_I','K_I','Ba_I','Total_Cations','F_A','Cl_A','OH_A','Xcel','XAlcel','Xprl','Xmrg','Xpg','Xms','XTriOct','Xphl','Xann','Xeas','Xsid','XDiOct'});
apfu=array2table(apfu,'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn2','Mg','Ca','Na','K','Ba','F','Cl','OH','Cation_Sum'});


end


