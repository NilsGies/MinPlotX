%%amphibole Structural Formula
% last modified 21.05.2024

function [StrctFrm, apfu,options_definition]=amphibole_Fe3known(data,headers,options)
%% define empty output variables
StrctFrm=[];
apfu=[];

%% options definition
options_definition=struct();

%Parameter 1
options_definition.Fe3_ratio.question='Enter Fe3_ratio';
options_definition.Fe3_ratio.Value=0;
options_definition.Fe3_ratio.limits=[0 1];

%Parameter 2
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
%create columns of zeros if optional data are not included
MC=zeros(m,14);
WT=zeros(m,14);

MC(:,1)=data(:,strcmp(headers,'SiO2'))./SiO2_mw; %for SiO2
WT(:,1)=data(:,strcmp(headers,'SiO2')); % SiO2 (wt %)

%calculates for TiO2 if it is included in the analysis
if any(strcmp(headers,'TiO2'))
    MC(:,2)=data(:,strcmp(headers,'TiO2'))./TiO2_mw; %for TiO2
    WT(:,2)=data(:,strcmp(headers,'TiO2')); % TiO2 (wt %)
end

MC(:,3)=(data(:,strcmp(headers,'Al2O3'))./Al2O3_mw).*2; %for Al2O3
WT(:,3)=data(:,strcmp(headers,'Al2O3')); % Al2O3 (wt %)

%calculates for Cr2O3 if it is included in the analysis
if any(strcmp(headers,'Cr2O3'))
    MC(:,4)=(data(:,strcmp(headers,'Cr2O3'))./Cr2O3_mw).*2; %for Cr2O3
    WT(:,4)=data(:,strcmp(headers,'Cr2O3')); % Cr2O3 (wt %)
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

    WT(:,5)=Fe3_wtper; % Fe2O3 (wt %)
    WT(:,6)=Fe2_wtper; % FeO (wt %)

elseif any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) %if FeO and Fe2O3 are both included as inputs

    MC(:,5)=(data(:,strcmp(headers,'Fe2O3'))./Fe2O3_mw).*2;
    WT(:,5)=data(:,strcmp(headers,'Fe2O3')); % Fe2O3 (wt %)
    MC(:,6)=data(:,strcmp(headers,'FeO'))./FeO_mw;
    WT(:,6)=data(:,strcmp(headers,'FeO')); % FeO (wt %)
else
    if any(strcmp(headers,'FeO')) %if FeO is FeO total

        %first calculate the weight percent of Fe2O3 and FeO from FeO total
        Fe3_wtper(:,1)=(data(:,strcmp(headers,'FeO')).*Fe_rat)./((2*FeO_mw)./Fe2O3_mw); %calculates Fe2O3 from Fe3+ ratio and total FeO
        Fe2_wtper(:,1)=data(:,strcmp(headers,'FeO'))-(data(:,strcmp(headers,'FeO')).*Fe_rat); %FeO = Total FeO - Fe3+ in Total FeO

        %calculate moles of cations
        MC(:,5)=(Fe3_wtper(:,1)./Fe2O3_mw).*2;
        MC(:,6)=Fe2_wtper(:,1)./FeO_mw;

        WT(:,5)=Fe3_wtper(:,1); % Fe2O3 (wt %)
        WT(:,6)=Fe2_wtper(:,1); % FeO (wt %)

    else %if Fe2O3 is Fe2O3 total

        %first calculate the weight percent of Fe2O3 and FeO from Fe2O3 total
        Fe3_wtper(:,1)=data(:,strcmp(headers,'Fe2O3'))-(data(:,strcmp(headers,'Fe2O3')).*(1-Fe_rat)); %Fe2O3 = Total Fe2O3 - Fe2+ in Total Fe2O3 
        Fe2_wtper(:,1)=(data(:,strcmp(headers,'Fe2O3')).*(1-Fe_rat)).*((2*FeO_mw)./Fe2O3_mw); %calculates FeO from Fe3+ ratio and total Fe2O3

        %calculate moles of cations
        MC(:,5)=(Fe3_wtper(:,1)./Fe2O3_mw).*2;
        MC(:,6)=Fe2_wtper(:,1)./FeO_mw;

        WT(:,5)=Fe3_wtper(:,1); % Fe2O3 (wt %)
        WT(:,6)=Fe2_wtper(:,1); % FeO (wt %)

    end
end

%calculates for MnO and/or Mn2O3
if any(strcmp(headers,'MnO')) && any(strcmp(headers,'Mn2O3')) %if MnO and Mn2O3 are both included as inputs

    %If MnO and Mn2O3 are both given, Mn2O3 is converted to MnO and
    %combined with MnO
    MC(:,7)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw)+data(:,strcmp(headers,'MnO')))./MnO_mw;
    WT(:,7)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw)+data(:,strcmp(headers,'MnO')));
else
    if any(strcmp(headers,'MnO')) %if MnO is MnO total

        MC(:,7)=data(:,strcmp(headers,'MnO'))./MnO_mw;
        WT(:,7)=data(:,strcmp(headers,'MnO'));

    elseif any(strcmp(headers,'Mn2O3')) %if Mn2O3 total is given, converted to MnO total

        %calculate moles of cations
        MC(:,7)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw))./MnO_mw;
        WT(:,7)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw));

    else
        %Mn2O3 and MnO are not included, Mn2 = 0
    end
end

MC(:,8)=data(:,strcmp(headers,'MgO'))./MgO_mw; %for MgO
WT(:,8)=data(:,strcmp(headers,'MgO')); % MgO (wt %)
MC(:,9)=data(:,strcmp(headers,'CaO'))./CaO_mw; %for CaO
WT(:,9)=data(:,strcmp(headers,'CaO')); % CaO (wt %)

%calculates for Na2O if it is included in the analysis
if any(strcmp(headers,'Na2O'))
    MC(:,10)=(data(:,strcmp(headers,'Na2O'))./Na2O_mw).*2; %for Na2O
    WT(:,10)=data(:,strcmp(headers,'Na2O')); %Na2O (wt %)
end

%calculates for K2O if it is included in the analysis
if any(strcmp(headers,'K2O'))
    MC(:,11)=(data(:,strcmp(headers,'K2O'))./K2O_mw).*2; %for K2O
    WT(:,11)=data(:,strcmp(headers,'K2O')); %K2O (wt %)
end

%calculates for H2O if it is included in the analysis
if any(strcmp(headers,'H2O'))
    MC(:,12)=(data(:,strcmp(headers,'H2O'))./H2O_mw).*2; %for H2O
    WT(:,12)=data(:,strcmp(headers,'H2O')); %H2O (wt %)
end

%calculates for F if it is included in the analysis
if any(strcmp(headers,'F'))
    MC(:,13)=data(:,strcmp(headers,'F'))./F_mw; %for F
    WT(:,13)=data(:,strcmp(headers,'F')); %F (wt %)
end

%calculates for Cl if it is included in the analysis
if any(strcmp(headers,'Cl'))
    MC(:,14)=data(:,strcmp(headers,'Cl'))./Cl_mw; %for Cl
    WT(:,14)=data(:,strcmp(headers,'Cl')); %Cl (wt %)
end



%% Oxygen Units and initial apfu normalization

O2_I(:,1)=MC(:,1).*2; %for SiO2
O2_I(:,2)=MC(:,2).*2; %for TiO2
O2_I(:,3)=MC(:,3).*(3/2); %for Al2O3
O2_I(:,4)=MC(:,4).*(3/2); %for Cr2O3
O2_I(:,5)=MC(:,5).*(3/2); %for Fe2O3
O2_I(:,6)=MC(:,6); %for FeO
O2_I(:,7)=MC(:,7); %for MnO
O2_I(:,8)=MC(:,8); %for MgO
O2_I(:,9)=MC(:,9); %for CaO
O2_I(:,10)=MC(:,10).*(1/2); %for Na2O
O2_I(:,11)=MC(:,11).*(1/2); %for K2O
O2_I(:,12)=MC(:,12).*(1/2); %H2O
O2_I(:,13)=MC(:,13); %F
O2_I(:,14)=MC(:,14); %Cl


%correction of O2 sum for F and Cl
Anion_Sum=sum(O2_I,2)-(0.5*O2_I(:,13))-(0.5*O2_I(:,14)); %anion sum, O2 equivalents minus 1/2 F and 1/2 Cl
Anhydrous_Sum=Anion_Sum-O2_I(:,12)-(0.5*O2_I(:,13))-(0.5*O2_I(:,14)); %sum not including F or OH

%normalizes all cations + F + Cl by the anion_sum and calculated to
%24 oxygen equivalents:

apfu_N24 = (MC./Anion_Sum).*24;

%correction of O2 sum for Ti-O substitution
Ti_adj(:,1)=zeros(m,1);

for c=1:m
    %if the moles of Ti > 8-(Si + Al) & 8-(Si + Al) > 0, then some Ti goes into the T site
    if apfu_N24(c,2) > 8-(apfu_N24(c,1)+apfu_N24(c,3)) && 8-(apfu_N24(c,1)+apfu_N24(c,3)) > 0
        %The Ti in T is subtracted from the Ti that is used to estimate OH
        Ti_adj(c,1) = apfu_N24(c,2) - (8-(apfu_N24(c,1)+apfu_N24(c,3)));
    else
        Ti_adj(c,1) = apfu_N24(c,2); %all Ti is used to estimate OH
    end
end

%correct for OH=2-2Ti: decides whether to use the O2 equivalent
%from the Ti-O estimate or 24-(F+Cl)
O2_max(:,1)=zeros(m,1);

for c=1:m
    if (23+Ti_adj(c,1))>24
        O2_max(c,1)=24-0.5*(apfu_N24(c,13)+apfu_N24(c,4));
    else
        O2_max(c,1)=23+Ti_adj(c,1);
    end
end

apfu_I=zeros(m,14);

%normalization of initial apfu
for c=1:m

    if O2_I(c,12) > 0 %if OH is included, the initial OH value is used 

        %selected procedure for apfu normalization if H2O is measured

        O2_min = 24 - 0.5.*(apfu_N24(c,13)+apfu_N24(c,4)); %calculates minimum number of O2 equivalent
        apfu_I(c,:) = (MC(c,:)./Anhydrous_Sum(c,1)).*O2_min; %normalization

    else

        %some analyses may use the Ti correct and others not
        if options.TiOH_correction.Value==true

            %normalization for Ti-OH correction
            apfu_I(c,:) = (MC(c,:)./Anhydrous_Sum(c,1)).*O2_max(c,1); %normalization

        else
            %if OH is not included and the OH=2-2Ti correction is not selected
            apfu_I(c,:)  = (MC(c,:)./Anhydrous_Sum(c,1)).*23; %normalization
        end
    end
end

%% Initial OH estimate

OH_prelim(:,1)=zeros(m,1);

for c=1:m
    %calculation of preliminary OH content
    if O2_I(c,12) > 0

        %if H2O is measured, OH new = OH old

    else
        if (apfu_I(c,12)+apfu_I(c,13)+apfu_I(c,14)) > 2 % is OH + F + Cl > 2?

            if (apfu_I(c,13)+apfu_I(c,14)) < 2 %is F + Cl < 2?

                OH_prelim(c,1) = 2 - (apfu_I(c,13)+apfu_I(c,14)); %OH = 2-(F+Cl)

            else

                OH_prelim(c,1) = 0; %otherwise the W site occupancy is exceeded and OH = 0

            end

        else

            if (apfu_I(c,12)+apfu_I(c,13)+apfu_I(c,14)) < 2 % is OH + F + Cl < 2?

                %is there space on W site after Cl and F are taken into
                %account?

                if 2-(apfu_I(c,13)+apfu_I(c,14))>0  %2-(F + Cl) > 0?

                    OH_prelim(c,1) = 2 - (apfu_I(c,13)+apfu_I(c,14)); %OH = 2-(F+Cl)
                else
                    OH_prelim(c,1) = 0; %otherwise the W site occupancy is exceeded and OH = 0
                end
            else
                OH_prelim(c,1) = apfu_I(c,12); %OH_preliminary estimate is the initial estimate,
                %may occur if OH is accidentally calculated already for some reason
            end
        end
    end
end

%calculation of preliminary OH content
Ticorr(:,1)=zeros(m,1);

for c=1:m
    %calculation of preliminary OH content
    if O2_I(c,12) > 0

        %nothing happens, H2O is already measured

    else

        %recalculate the Ti correction from normalized initial apfu

        %if Ti > (8-(Si+Al)) & (8-(Si+Al))>0
        if apfu_I(c,2)>(8-(apfu_I(c,1)+apfu_I(c,3))) && (8-(apfu_I(c,1)+apfu_I(c,3))) > 0
            Ticorr(c,1)=apfu_I(c,2)-(8-(apfu_I(c,1)+apfu_I(c,3))); %Ti in M = Titotal - Ti in T
        else
            Ticorr(c,1)=apfu_I(c,2); %otherwise all Ti is in M
        end


        % Determine appropriate OH correction
        if options.TiOH_correction.Value==true
            if (apfu_I(c,12) + apfu_I(c,13) + apfu_I(c,14)) > 2 %If (OH + F + Cl) > 2
                if ((2.*apfu_I(c,12))./((apfu_I(c,12)+ apfu_I(c,13) + apfu_I(c,14))-2.*apfu_I(c,2)))>0 %if ((2*OH)/((OH+F+Cl)-2*Ti)) >0
                    apfu_I(c,12)=((2.*apfu_I(c,12))./((apfu_I(c,12)+ apfu_I(c,13) + apfu_I(c,14))-2.*apfu_I(c,2)));
                else
                    apfu_I(c,12)=0; %W is filled by F, Cl, and O
                end
            else
                if (apfu_I(c,12) + apfu_I(c,13) + apfu_I(c,14)) < 2 %If (OH + F + Cl) < 2
                    if ((2-apfu_I(c,13)-apfu_I(c,14))-2.*Ticorr(c,1))>0 %if 2-F-Cl-2*Ti > 0
                        apfu_I(c,12)=((2-apfu_I(c,13)-apfu_I(c,14))-2.*Ticorr(c,1)); %OH=(2-F-Cl-2*Ti)
                    else
                        apfu_I(c,12)=0; %W is filled by F, Cl, and O
                    end
                else
                    apfu_I(c,12)=apfu_I(c,12); %OH content is not adjusted for Ti
                end
            end
        else
            %OH is estimated from 2-(F+Cl) without Ti corection
            apfu_I(c,12)=OH_prelim(c,1);
        end
    end
end

%Moles of anions
apfu_an(:,1)=apfu_I(:,1).*2; %SiO2
apfu_an(:,2)=apfu_I(:,2).*2; %TiO2
apfu_an(:,3)=apfu_I(:,3).*1.5; %Al2O3
apfu_an(:,4)=apfu_I(:,4).*1.5; %Cr2O3
apfu_an(:,5)=apfu_I(:,5).*1.5; %Fe2O3
apfu_an(:,6)=apfu_I(:,6); %FeO
apfu_an(:,7)=apfu_I(:,7); %MnO
apfu_an(:,8)=apfu_I(:,8); %MgO
apfu_an(:,9)=apfu_I(:,9); %CaO
apfu_an(:,10)=apfu_I(:,10).*0.5; %Na2O
apfu_an(:,11)=apfu_I(:,11).*0.5; %K2O
apfu_an(:,12)=apfu_I(:,12).*0.5; %H2O

%for F
for c=1:m
    if (apfu_I(c,13)+apfu_I(c,14))>2 %F + Cl cannot > 2
        apfu_an(c,13)=2.*(apfu_I(c,13)./(apfu_I(c,13)+apfu_I(c,14))); %scales F
    else
        apfu_an(c,13)=apfu_I(c,13); %Does not scale F
    end
end

%for Cl
for c=1:m
    if (apfu_I(c,13)+apfu_I(c,14))>2 %F + Cl cannot > 2
        apfu_an(c,14)=2.*(apfu_I(c,14)./(apfu_I(c,13)+apfu_I(c,14))); %scales Cl
    else
        apfu_an(c,14)=apfu_I(c,14); %Does not scale Cl
    end
end

%the anion sum may be > 24
anion_sum=sum(apfu_an(:,1:12),2)+0.5*apfu_an(:,13)+0.5*apfu_an(:,14);
anion_norm=24./anion_sum;

apfu_n=apfu_I.*anion_norm; %cations normalized to 24 anions
apfu_ann=apfu_an.*anion_norm; %anions normalized to 24 anions

%molar mass (g/mol) = Si*SiO2_mw + Ti*TiO2_mw*0.5 + Al*Al2O3_mw*0.5 +
%Cr*Cr2O3_mw*0.5 + Fe3*Fe2O3_mw*0.5 + Fe2*FeO_mw + Mn2*MnO_mw + Mg*MgO_mw +
%Ca*CaO_mw + Na*Na2O_mw*0.5 + K*K2O_mw*0.5 + H*H2O_mw*0.5 + F*F_mw +
%Cl*Cl_mw - Cl*15.9994.*0.5  - F*15.9994.*0.5 


molar_mass=apfu_n(:,1).*SiO2_mw+apfu_n(:,2).*TiO2_mw+(apfu_n(:,3)./2).*Al2O3_mw+(apfu_n(:,4)./2).*Cr2O3_mw...
+(apfu_n(:,5)./2).*Fe2O3_mw+apfu_n(:,6).*FeO_mw+apfu_n(:,7).*MnO_mw+apfu_n(:,8).*MgO_mw+apfu_n(:,9).*CaO_mw...
+(apfu_n(:,10)./2).*Na2O_mw+(apfu_n(:,11)./2).*K2O_mw+(apfu_n(:,12)./2).*H2O_mw+apfu_ann(:,13).*F_mw...
+apfu_ann(:,14).*Cl_mw-apfu_ann(:,13).*15.9994.*0.5-apfu_ann(:,14).*15.9994.*0.5;

%new wt% of oxides with or without estimated H2O
D2=WT;

for c=1:m
    if O2_I(c,12) > 0 %is OH in the inital input > 0
        D2(c,12)=WT(c,12); %keeps original H2O content if measured
    else
        %uses estimated initial OH content
        D2(c,12)=apfu_n(c,12).*H2O_mw.*0.5.*(1./molar_mass(c,1))*100;
    end
end

%% Itterates the calculation 10 times
for z=1:10

    MC2=zeros(m,14);
    MC2(:,1)=D2(:,1)./SiO2_mw; %for SiO2
    MC2(:,2)=D2(:,2)./TiO2_mw; %for TiO2
    MC2(:,3)=(D2(:,3)./Al2O3_mw).*2; %for Al2O3
    MC2(:,4)=(D2(:,4)./Cr2O3_mw).*2; %for Cr2O3
    MC2(:,5)=(D2(:,5)./Fe2O3_mw).*2;
    MC2(:,6)=D2(:,6)./FeO_mw;
    MC2(:,7)=D2(:,7)./MnO_mw;
    MC2(:,8)=D2(:,8)./MgO_mw; %for MgO
    MC2(:,9)=D2(:,9)./CaO_mw; %for CaO
    MC2(:,10)=(D2(:,10)./Na2O_mw).*2; %for Na2O
    MC2(:,11)=(D2(:,11)./K2O_mw).*2; %for K2O
    MC2(:,12)=(D2(:,12)./H2O_mw).*2; %for H2O
    MC2(:,13)=D2(:,13)./F_mw; %for F
    MC2(:,14)=D2(:,14)./Cl_mw; %for Cl

    %% Oxygen Units and initial apfu normalization

    O2(:,1)=MC2(:,1).*2; %for SiO2
    O2(:,2)=MC2(:,2).*2; %for TiO2
    O2(:,3)=MC2(:,3).*(3/2); %for Al2O3
    O2(:,4)=MC2(:,4).*(3/2); %for Cr2O3
    O2(:,5)=MC2(:,5).*(3/2); %for Fe2O3
    O2(:,6)=MC2(:,6); %for FeO
    O2(:,7)=MC2(:,7); %for MnO
    O2(:,8)=MC2(:,8); %for MgO
    O2(:,9)=MC2(:,9); %for CaO
    O2(:,10)=MC2(:,10).*(1/2); %for Na2O
    O2(:,11)=MC2(:,11).*(1/2); %for K2O
    O2(:,12)=MC2(:,12).*(1/2); %H2O
    O2(:,13)=MC2(:,13); %F
    O2(:,14)=MC2(:,14); %Cl

    %correction of O2 sum for F and Cl
    Anion_Sum2=sum(O2,2)-(0.5*O2(:,13))-(0.5*O2(:,14)); %anion sum, O2 equivalents minus 1/2 F and 1/2 Cl
    Anhydrous_Sum2=Anion_Sum2-O2(:,12)-(0.5*O2(:,13))-(0.5*O2(:,14)); %sum not including F or OH

    %normalizes all cations + F + Cl by the anion_sum and calculated to
    %24 oxygen equivalents:

    apfu_N24_2 = (MC2./Anion_Sum2).*24;

    %correction of O2 sum for Ti-O substitution
    Ti_adj2(:,1)=zeros(m,1);

    for c=1:m
        %if the moles of Ti > 8-(Si + Al) & 8-(Si + Al) > 0, then some Ti goes into the T site
        if apfu_N24_2(c,2) > 8-(apfu_N24_2(c,1)+apfu_N24_2(c,3)) && 8-(apfu_N24_2(c,1)+apfu_N24_2(c,3)) > 0
            %The Ti in T is subtracted from the Ti that is used to estimate OH
            Ti_adj2(c,1) = apfu_N24_2(c,2) - (8-(apfu_N24_2(c,1)+apfu_N24_2(c,3)));
        else
            Ti_adj2(c,1) = apfu_N24_2(c,2); %all Ti is used to estimate OH
        end
    end

    %correct for OH=2-2Ti: decides whether to use the O2 equivalent
    %from the Ti-O estimate or 24-(F+Cl)
    O2_max2(:,1)=zeros(m,1);

    for c=1:m
        if (23+Ti_adj2(c,1))>24
            O2_max2(c,1)=24-0.5*(apfu_N24_2(c,13)+apfu_N24_2(c,4));
        else
            O2_max2(c,1)=23+Ti_adj2(c,1);
        end
    end

    apfu=zeros(m,14);

    %normalization of initial apfu
    for c=1:m

        if O2_I(c,12) > 0 %if OH is included, the initial OH value is used

            %selected procedure for apfu normalization if H2O is measured

            O2_min2 = 24 - 0.5.*(apfu_N24_2(c,13)+apfu_N24_2(c,4)); %calculates minimum number of O2 equivalent
            apfu(c,:) = (MC2(c,:)./Anhydrous_Sum2(c,1)).*O2_min2; %normalization

        else

            %some analyses may use the Ti correct and others not
            if options.TiOH_correction.Value==true

                %normalization for Ti-OH correction
                apfu(c,:) = (MC2(c,:)./Anhydrous_Sum2(c,1)).*O2_max2(c,1); %normalization

            else
                %if OH is not included and the OH=2-2Ti correction is not selected
                apfu(c,:)  = (MC2(c,:)./Anhydrous_Sum2(c,1)).*23; %normalization
            end
        end
    end

    % Initial OH estimate
    OH_prelim2(:,1)=zeros(m,1);

    for c=1:m
        %calculation of preliminary OH content
        if O2_I(c,12) > 0

            %if H2O is measured, OH new = OH old

        else
            if (apfu(c,12)+apfu(c,13)+apfu(c,14)) > 2 % is OH + F + Cl > 2?

                if (apfu(c,13)+apfu(c,14)) < 2 %is F + Cl < 2?

                    OH_prelim2(c,1) = 2 - (apfu(c,13)+apfu(c,14)); %OH = 2-(F+Cl)

                else

                    OH_prelim2(c,1) = 0; %otherwise the W site occupancy is exceeded and OH = 0

                end

            else

                if (apfu(c,12)+apfu(c,13)+apfu(c,14)) < 2 % is OH + F + Cl < 2?

                    %is there space on W site after Cl and F are taken into
                    %account?

                    if 2-(apfu(c,13)+apfu(c,14))>0  %2-(F + Cl) > 0?

                        OH_prelim2(c,1) = 2 - (apfu(c,13)+apfu(c,14)); %OH = 2-(F+Cl)
                    else
                        OH_prelim2(c,1) = 0; %otherwise the W site occupancy is exceeded and OH = 0
                    end
                else
                    OH_prelim2(c,1) = apfu(c,12); %OH_preliminary estimate is the initial estimate,
                    %may occur if OH is accidentally calculated already for some reason
                end
            end
        end
    end

    %calculation of preliminary OH content
    Ticorr2(:,1)=zeros(m,1);

    for c=1:m
        %calculation of preliminary OH content
        if O2_I(c,12) > 0

            %nothing happens, H2O is already measured

        else

            %recalculate the Ti correction from normalized initial apfu

            %if Ti > (8-(Si+Al)) & (8-(Si+Al))>0
            if apfu(c,2)>(8-(apfu(c,1)+apfu(c,3))) && (8-(apfu(c,1)+apfu(c,3))) > 0
                Ticorr2(c,1)=apfu(c,2)-(8-(apfu(c,1)+apfu(c,3))); %Ti in M = Titotal - Ti in T
            else
                Ticorr2(c,1)=apfu(c,2); %otherwise all Ti is in M
            end


            % Determine appropriate OH correction
            if options.TiOH_correction.Value==true
                if (apfu(c,12) + apfu(c,13) + apfu(c,14)) > 2 %If (OH + F + Cl) > 2
                    if ((2.*apfu(c,12))./((apfu(c,12)+ apfu(c,13) + apfu(c,14))-2.*apfu(c,2)))>0 %if ((2*OH)/((OH+F+Cl)-2*Ti)) >0
                        apfu(c,12)=((2.*apfu(c,12))./((apfu(c,12)+ apfu(c,13) + apfu(c,14))-2.*apfu(c,2)));
                    else
                        apfu(c,12)=0; %W is filled by F, Cl, and O
                    end
                else
                    if (apfu(c,12) + apfu(c,13) + apfu(c,14)) < 2 %If (OH + F + Cl) < 2
                        if ((2-apfu(c,13)-apfu(c,14))-2.*Ticorr2(c,1))>0 %if 2-F-Cl-2*Ti > 0
                            apfu(c,12)=((2-apfu(c,13)-apfu(c,14))-2.*Ticorr2(c,1)); %OH=(2-F-Cl-2*Ti)
                        else
                            apfu(c,12)=0; %W is filled by F, Cl, and O
                        end
                    else
                        apfu(c,12)=apfu(c,12); %OH content is not adjusted for Ti
                    end
                end
            else
                %OH is estimated from 2-(F+Cl) without Ti corection
                apfu(c,12)=OH_prelim2(c,1);
            end
        end
    end

    %Moles of anions
    apfu_an2(:,1)=apfu(:,1).*2; %SiO2
    apfu_an2(:,2)=apfu(:,2).*2; %TiO2
    apfu_an2(:,3)=apfu(:,3).*1.5; %Al2O3
    apfu_an2(:,4)=apfu(:,4).*1.5; %Cr2O3
    apfu_an2(:,5)=apfu(:,5).*1.5; %Fe2O3
    apfu_an2(:,6)=apfu(:,6); %FeO
    apfu_an2(:,7)=apfu(:,7); %MnO
    apfu_an2(:,8)=apfu(:,8); %MgO
    apfu_an2(:,9)=apfu(:,9); %CaO
    apfu_an2(:,10)=apfu(:,10).*0.5; %Na2O
    apfu_an2(:,11)=apfu(:,11).*0.5; %K2O
    apfu_an2(:,12)=apfu(:,12).*0.5; %H2O

    %for F
    for c=1:m
        if (apfu(c,13)+apfu(c,14))>2 %F + Cl cannot > 2
            apfu_an2(c,13)=2.*(apfu(c,13)./(apfu(c,13)+apfu(c,14))); %scales F
        else
            apfu_an2(c,13)=apfu(c,13); %Does not scale F
        end
    end

    %for Cl
    for c=1:m
        if (apfu(c,13)+apfu(c,14))>2 %F + Cl cannot > 2
            apfu_an2(c,14)=2.*(apfu(c,14)./(apfu(c,13)+apfu(c,14))); %scales Cl
        else
            apfu_an2(c,14)=apfu(c,14); %Does not scale Cl
        end
    end

    %the anion sum may be > 24
    anion_sum2=sum(apfu_an2(:,1:12),2)+0.5*apfu_an2(:,13)+0.5*apfu_an2(:,14);
    anion_norm2=24./anion_sum2;

    apfu_n2=apfu.*anion_norm2; %cations normalized to 24 anions
    apfu_ann2=apfu_an2.*anion_norm2; %anions normalized to 24 anions

    %molar mass (g/mol) = Si*SiO2_mw + Ti*TiO2_mw*0.5 + Al*Al2O3_mw*0.5 +
    %Cr*Cr2O3_mw*0.5 + Fe3*Fe2O3_mw*0.5 + Fe2*FeO_mw + Mn2*MnO_mw + Mg*MgO_mw +
    %Ca*CaO_mw + Na*Na2O_mw*0.5 + K*K2O_mw*0.5 + H*H2O_mw*0.5 + F*F_mw +
    %Cl*Cl_mw - Cl*15.9994.*0.5  - F*15.9994.*0.5


    molar_mass2=apfu_n2(:,1).*SiO2_mw+apfu_n2(:,2).*TiO2_mw+(apfu_n2(:,3)./2).*Al2O3_mw+(apfu_n2(:,4)./2).*Cr2O3_mw...
        +(apfu_n2(:,5)./2).*Fe2O3_mw+apfu_n2(:,6).*FeO_mw+apfu_n2(:,7).*MnO_mw+apfu_n2(:,8).*MgO_mw+apfu_n2(:,9).*CaO_mw...
        +(apfu_n2(:,10)./2).*Na2O_mw+(apfu_n2(:,11)./2).*K2O_mw+(apfu_n2(:,12)./2).*H2O_mw+apfu_ann2(:,13).*F_mw...
        +apfu_ann2(:,14).*Cl_mw-apfu_ann2(:,13).*15.9994.*0.5-apfu_ann2(:,14).*15.9994.*0.5;

    %new wt% of oxides
    D2=WT;

    for c=1:m
        if O2_I(c,12) > 0 %is OH in the inital input > 0
            D2(c,12)=WT(c,12); %keeps original H2O content if measured
        else
            %uses estimated initial OH content
            D2(c,12)=apfu_n2(c,12).*H2O_mw.*0.5.*(1./molar_mass2(c,1))*100;
        end
    end
end

%% structural Formula  & outputs

StrctFrm=StrctFrm_amph_fun(apfu_n2);

%W cations
W_site(:,1)=apfu_n2(:,12); %OH
W_site(:,2)=apfu_n2(:,13); %F
W_site(:,3)=apfu_n2(:,14); %Cl
W_site(:,4)=apfu_n2(:,12)+apfu_n2(:,13)+apfu_n2(:,14); %W sum

all=[StrctFrm W_site];

%limit on significant digits (eliminates rounding noise)
all(all<1e-6) = 0;
apfu_n2(apfu_n2<1e-6) = 0;

StrctFrm=array2table(all,'VariableNames',{'Si_T','Al_T','Sum_T','Al_C','Ti_C','Cr_C','Fe3_C','Mg_C','Fe2_C','Mn_C','Sum_C','Mg_B','Fe2_B','Mn_B','Ca_B','Na_B','Sum_B','Ca_A','Na_A','K_A','Cation_Sum','OH_W','F_W','Cl_W','W_sum'});
apfu=array2table(apfu_n2,'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K','OH','F','Cl'});


end










