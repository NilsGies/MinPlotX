%%Garnet Structural Formula
% last modified 01.08.2024

function [StrctFrm, apfu]=garnet_Fe3known(data,headers,options)

[m,~]=size(data); %finds the x and y size of the input data matrix

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
Mn2O3_mw=157.873;
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

MC(:,1)=data(:,strcmp(headers,'SiO2'))./SiO2_mw; %for SiO2

%adds a column of zeros if Ti is not included in the calculation
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
    MC(:,7)=Fe2_wtper(:,1)./FeO_mw;

elseif any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) %if FeO and Fe2O3 are both included as inputs

     MC(:,4)=(data(:,strcmp(headers,'Fe2O3'))./Fe2O3_mw).*2;
     MC(:,7)=data(:,strcmp(headers,'FeO'))./FeO_mw;
     
else
    if any(strcmp(headers,'FeO')) %if FeO is FeO total

        %first calculate the weight percent of Fe2O3 and FeO from FeO total
        Fe3_wtper(:,1)=(data(:,strcmp(headers,'FeO')).*Fe_rat)./((2*FeO_mw)./Fe2O3_mw); %calculates Fe2O3 from Fe3+ ratio and total FeO
        Fe2_wtper(:,1)=data(:,strcmp(headers,'FeO'))-(data(:,strcmp(headers,'FeO')).*Fe_rat); %FeO = Total FeO - Fe3+ in Total FeO

        %calculate moles of cations
        MC(:,4)=(Fe3_wtper(:,1)./Fe2O3_mw).*2;
        MC(:,7)=Fe2_wtper(:,1)./FeO_mw;

    else %if Fe2O3 is Fe2O3 total

        %first calculate the weight percent of Fe2O3 and FeO from Fe2O3 total
        Fe3_wtper(:,1)=data(:,strcmp(headers,'Fe2O3'))-(data(:,strcmp(headers,'Fe2O3')).*(1-Fe_rat)); %Fe2O3 = Total Fe2O3 - Fe2+ in Total Fe2O3 
        Fe2_wtper(:,1)=(data(:,strcmp(headers,'Fe2O3')).*(1-Fe_rat)).*((2*FeO_mw)./Fe2O3_mw); %calculates FeO from Fe3+ ratio and total Fe2O3

        %calculate moles of cations
        MC(:,4)=(Fe3_wtper(:,1)./Fe2O3_mw).*2;
        MC(:,7)=Fe2_wtper(:,1)./FeO_mw;

    end
end

%adds a column of zeros if Cr is not included in the calculation
if any(strcmp(headers,'Cr2O3'))
    MC(:,5)=(data(:,strcmp(headers,'Cr2O3'))./Cr2O3_mw).*2; %for Cr2O3
end

%adds a column of zeros if Y is not included in the calculation
if any(strcmp(headers,'Y2O3'))
    MC(:,6)=(data(:,strcmp(headers,'Y2O3'))./Y2O3_mw).*2; %for Y2O3
end


%calculates for MnO and/or Mn2O3
if any(strcmp(headers,'MnO')) && any(strcmp(headers,'Mn2O3')) %if MnO and Mn2O3 are both included as inputs 

     %If MnO and Mn2O3 are both given, Mn2O3 is converted to MnO and
     %combined with MnO
     MC(:,8)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw)+data(:,strcmp(headers,'MnO')))./MnO_mw;
else
    if any(strcmp(headers,'MnO')) %if MnO is MnO total

        MC(:,8)=data(:,strcmp(headers,'MnO'))./MnO_mw;

    else  %if Mn2O3 total is given, converted to MnO total

        %calculate moles of cations
        MC(:,8)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw))./MnO_mw;
    end
end

MC(:,9)=data(:,strcmp(headers,'MgO'))./MgO_mw; %for MgO
MC(:,10)=data(:,strcmp(headers,'CaO'))./CaO_mw; %for CaO

%adds a column of zeros if Na is not included in the calculation
if any(strcmp(headers,'Na2O'))
    MC(:,11)=(data(:,strcmp(headers,'Na2O'))./Na2O_mw).*2; %for Na2O
end

%% Calculate Oxygen Units

O2(:,1)=MC(:,1).*2; %for SiO2
O2(:,2)=MC(:,2).*2; %for TiO2
O2(:,3)=MC(:,3).*(3/2); %for Al2O3
O2(:,4)=MC(:,4).*(3/2); %for Fe2O3
O2(:,5)=MC(:,5).*(3/2); %for Cr2O3
O2(:,6)=MC(:,6).*(3/2); %for Y2O3
O2(:,7)=MC(:,7); %for FeO
O2(:,8)=MC(:,8); %for MnO
O2(:,9)=MC(:,9); %for MgO
O2(:,10)=MC(:,10); %for CaO
O2(:,11)=MC(:,11)./2; %for Na2O


O2total=sum(O2,2); %O2 totals
MCnormfact=Opfu./sum(O2,2); %normalization factor

%% Atoms pfu

apfu=zeros(m,13); %create columns of zeros if optional data are not included
%apfu=MCnormfact.*MC; %creates a matrix of normalized cations
apfu(:,1:11)=MCnormfact.*MC; %creates a matrix of normalized cations
apfu(:,12)=sum(apfu,2); %calculations the total, which should be close to 8
O2_def=apfu(:,13); %Makes O2_def a column of zeros


%% structural formula calculation

StrctFrm=zeros(m,20);

%T SITE
%Si
for c=1:m
    if apfu(c,1)<3.000
        StrctFrm(c,1)=apfu(c,1); %If Si < 3, then Si(T) = the measured Si content
    else
        StrctFrm(c,1)=3; %If Si is in excess, then Si(T) = 3
    end
end

%Al(T)
for c=1:m
    if 3-StrctFrm(c,1)>0 %Is 3-Si > 0? If y, then some Al goes into T
        if 3-StrctFrm(c,1)>apfu(c,3) %For low Al Grt, 3-Si may be > Al
            StrctFrm(c,2)=apfu(c,3); %All Al goes into T
        else
            StrctFrm(c,2)=3-StrctFrm(c,1); %if there isn't enough space in T for all Al, the rest will go to Y
        end
    else
        StrctFrm(c,2)=0; %if Si=3, then no Al goes into T
    end
end

%Fe3+(T)
for c=1:m
    if 3-StrctFrm(c,1)-StrctFrm(c,2)>0 %Is 3-(Si+Al) > 0? If y, then some Fe3+ goes into T
        if 3-StrctFrm(c,1)-StrctFrm(c,2)>apfu(c,4) %For low Fe3+ grt, 3-(Si+Al) may be > Fe3+
            StrctFrm(c,3)=apfu(c,4); %All Fe3+ goes into T
        else
            StrctFrm(c,3)=3-StrctFrm(c,1)-StrctFrm(c,2); %if there isn't enough space in T for all Fe3+, the rest will go to M1
        end
    else
        StrctFrm(c,3)=0; %if Si+Al=2, then no Fe3+ goes into T
    end
end

%Sum of T site
StrctFrm(:,4)=StrctFrm(:,1)+StrctFrm(:,2)+StrctFrm(:,3);

%Y SITE

%Si (Y)
for c=1:m
    if apfu(c,1)<3.000
        StrctFrm(c,5)=0; %If Si < 3, then there is no Si on the octahedral site
    else
        StrctFrm(c,5)=apfu(c,1)-3; %If Si is in excess, then some Si is assigned to the octahedral site
    end
end

%Al (Y)
StrctFrm(:,6)=apfu(:,3)-StrctFrm(:,2); %Al(M1) = Total Al - Al(T)

%Ti (Y)
StrctFrm(:,7)=apfu(:,2);

%Cr (Y)
StrctFrm(:,8)=apfu(:,5);

%Fe3+ (Y)
StrctFrm(:,9)=apfu(:,4)-StrctFrm(:,3); %Fe3+(M1) = Total Fe3+ - Fe3+(T)

%Mg (Y)
for c=1:m
    if (StrctFrm(c,5)+StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8)+StrctFrm(c,9))<2.000 %Mg is only considered if the octahedral site is not yet filled
        if (2-(StrctFrm(c,5)+StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8)+StrctFrm(c,9))) > apfu(c,9)
            StrctFrm(c,10)=apfu(c,9); %all Mg goes into octahedral site
        else
            StrctFrm(c,10)=2-(StrctFrm(c,5)+StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8)+StrctFrm(c,9)); % only some Mg goes into the octahedral site
        end
    else
        StrctFrm(c,10)=0; % no Mg goes into the octahedral site (it's already filled)
    end
end

%Fe2+ (Y)
for c=1:m
    if (StrctFrm(c,5)+StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8)+StrctFrm(c,9)+StrctFrm(c,10))<2.000 %Fe2+ is only considered if the octahedral site is not yet filled
        if (2-(StrctFrm(c,5)+StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8)+StrctFrm(c,9)+StrctFrm(c,10))) > apfu(c,7)
            StrctFrm(c,11)=apfu(c,7); %all Fe2+ goes into octahedral site
        else
            StrctFrm(c,11)=2-(StrctFrm(c,5)+StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8)+StrctFrm(c,9)+StrctFrm(c,10)); % only some Fe2+ goes into the octahedral site
        end
    else
        StrctFrm(c,11)=0; % no Fe2+ goes into the octahedral site (it's already filled)
    end
end

%Mn (Y)
for c=1:m
    if (StrctFrm(c,5)+StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8)+StrctFrm(c,9)+StrctFrm(c,10)+StrctFrm(c,11))<2.000 %Mn is only considered if the octahedral site is not yet filled
        if (2-(StrctFrm(c,5)+StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8)+StrctFrm(c,9)+StrctFrm(c,10)+StrctFrm(c,11))) > apfu(c,8)
            StrctFrm(c,12)=apfu(c,8); %all Mn goes into octahedral site
        else
            StrctFrm(c,12)=2-(StrctFrm(c,5)+StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8)+StrctFrm(c,9)+StrctFrm(c,10)+StrctFrm(c,11)); % only some Mn goes into the octahedral site
        end
    else
        StrctFrm(c,12)=0; % no Mn goes into the octahedral site (it's already filled)
    end
end


%Y sum
StrctFrm(:,13)=StrctFrm(:,12)+StrctFrm(:,11)+StrctFrm(:,10)+StrctFrm(:,9)+StrctFrm(:,8)+StrctFrm(:,7)+StrctFrm(:,6)+StrctFrm(:,5);

%X SITE

%Y (X)
StrctFrm(:,14)=apfu(:,6);

%Mg (X)
StrctFrm(:,15)=apfu(:,9)-StrctFrm(:,10);

%Fe (X)
StrctFrm(:,16)=apfu(:,7)-StrctFrm(:,11);

%Mn (X)
StrctFrm(:,17)=apfu(:,8)-StrctFrm(:,12);

%Ca (X)
StrctFrm(:,18)=apfu(:,10);

%Na (X)
StrctFrm(:,19)=apfu(:,11);

% X Sum
StrctFrm(:,20)=StrctFrm(:,19)+StrctFrm(:,18)+StrctFrm(:,17)+StrctFrm(:,16)+StrctFrm(:,15)+StrctFrm(:,14);

%% end member calculations

A(:,1)=apfu(:,10); %Ca
A(:,2)=apfu(:,9); %Mg
A(:,3)=apfu(:,4)+apfu(:,7); %Fetotal
A(:,4)=apfu(:,5); %Cr
A(:,5)=apfu(:,8); %Mn
A(:,6)=apfu(:,3); %Al
AT=transpose(A); %transpose of A

M=[0 0 0 3 3 3; 0 3 0 0 0 0; 3 0 0 0 2 0; 0 0 0 0 0 2; 0 0 3 0 0 0; 2 2 2 2 0 0];

X=zeros(6,m);
for c=1:m
    X(:,c)=inv(M)*AT(:,c); %calculates endmembers
end

Xtot=sum(X);%sum of endmembers
Xnorm=X./sum(X); %normalizes the endmembers to 1
XnormT=transpose(Xnorm); %transposes back

Endmembers(:,1)=XnormT(:,1); % XAlm
Endmembers(:,2)=XnormT(:,2); % XPrp
Endmembers(:,3)=XnormT(:,3); % XSps
Endmembers(:,4)=XnormT(:,4); % XGrs
Endmembers(:,5)=XnormT(:,5); % XAdr
Endmembers(:,6)=XnormT(:,6); % XUv

Endmembers(Endmembers<1e-3) = 0; %limit on endmember noise (cannot be less than a fraction of a percent)

%Final outputs

all=[StrctFrm Endmembers O2_def];

%limit on significant digits (eliminates rounding noise)
all(all<1e-6) = 0;
apfu(apfu<1e-6) = 0;
all(isnan(all))=0; %remove NaN

StrctFrm=array2table(all,'VariableNames',{'Si_T','Al_T','Fe3_T','T_Sum','Si_Oct','Al_Oct','Ti_Oct','Cr_Oct','Fe3_Oct','Mg_Oct','Fe2_Oct','Mn2_Oct','Oct_Sum','Y_Dod','Mg_Dod','Fe_Dod','Mn2_Dod','Ca_Dod','Na_Dod','Dod_Sum','Xalm','Xprp','Xsps','Xgrs','Xadr','Xuv','O2_deficiency'});
apfu=array2table(apfu,'VariableNames',{'Si','Ti','Al','Fe3','Cr','Y','Fe2','Mn2','Mg','Ca','Na','Cation_Sum','O2_deficiency'});
apfu.Fe3_FeTot=apfu.Fe3./(apfu.Fe3+apfu.Fe2);

end
