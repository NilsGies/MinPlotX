%%Calculates amphibole formula following Hawthorne et al. (2012) and
%Leake et al., (1997)
%last modified 11.07.2024

function [StrctFrm, apfu,options_definition, apfu_FeT, Fe3_limits, Fe3_class, apfu_SiT, apfu_AfullT, apfu_NaAT, apfu_Fe2O3T, apfu_SiAlT, apfu_KAT, apfu_CaNaBT]=amphibole_MinPlotX(data,headers,options)

%% define empty output variables
StrctFrm=[];
apfu=[];
Fe3_limits=[];
Fe3_class=[];
apfu_SiT=[];
apfu_AfullT=[];
apfu_NaAT=[];
apfu_Fe2O3T=[];
apfu_SiAlT=[];
apfu_KAT=[];
apfu_CaNaBT=[];
apfu_FeT=[];
%% options definition
%Parameter 1
options_definition.want_ferric_Fe3.question='Calculate Fe3+ from stoichiometry?';
options_definition.want_ferric_Fe3.options={true,false}; % optional
options_definition.want_ferric_Fe3.Value=true; %default_value

%Parameter 2
options_definition.TiOH_correction.question='Do you want to correct for OH=2-2Ti?';
options_definition.TiOH_correction.Value=true; %default_value
options_definition.TiOH_correction.options={true,false}; % optional

%Parameter 3
options_definition.UseKnownFe3.question='Is Fe3+/FeTotal Ratio known?';
options_definition.UseKnownFe3.Value=false;

%Parameter 4
options_definition.Fe3_ratio.question='Fe3+/FeTotal Ratio?';
options_definition.Fe3_ratio.Value=0; %default_value
options_definition.Fe3_ratio.limits=[0 1]; % optional

%Parameter 5
options_definition.RecalcTotalFe.question='Recalculate FeO and Fe2O3?';
options_definition.RecalcTotalFe.Value=false;


if any(strcmp(headers,'H2O'))
    options.H2Oknown.Value='y';
else
    options.H2Oknown.Value='n';
end

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
H2O_mw=18.015;

W=[SiO2_mw,TiO2_mw,Al2O3_mw,Cr2O3_mw,Fe2O3_mw,FeO_mw,MnO_mw,MgO_mw,CaO_mw,Na2O_mw,K2O_mw,H2O_mw,F_mw,Cl_mw];

%% makes a new array, D, where the oxides are in the correct order and
%filling in zeros where optional oxides are not included 

D=zeros(m,14); %create columns of zeros if optional data are not included

D(:,1)=data(:,strcmp(headers,'SiO2')); %for SiO2

%calculates for TiO2 if it is included in the analysis
if any(strcmp(headers,'TiO2'))
    D(:,2)=data(:,strcmp(headers,'TiO2')); %for TiO2
end

D(:,3)=data(:,strcmp(headers,'Al2O3')); %for Al2O3

%calculates for Cr2O3 if it is included in the analysis
if any(strcmp(headers,'Cr2O3'))
    D(:,4)=data(:,strcmp(headers,'Cr2O3')); %for Cr2O3
end 

%Fe2O3 (column 5) is blank for D

%converts FeO and/or Fe2O3 to FeO total 
if any(strcmp(headers,'FeO')) && any(strcmp(headers,'Fe2O3')) %if FeO and Fe2O3 are both included as inputs 

     %If FeO and Fe2O3 are both given, Fe2O3 is converted to FeO and
     %combined with FeO
     D(:,6)=(data(:,strcmp(headers,'Fe2O3')).*((2*FeO_mw)./Fe2O3_mw)+data(:,strcmp(headers,'FeO')));
else
    if any(strcmp(headers,'FeO')) %if FeO is FeO total

        D(:,6)=data(:,strcmp(headers,'FeO'));

    else %if Fe2O3 total is given, converted to FeO total

        %calculate moles of cations
        D(:,6)=(data(:,strcmp(headers,'Fe2O3')).*((2*FeO_mw)./Fe2O3_mw));

    end
end

%calculates for MnO and/or Mn2O3
if any(strcmp(headers,'MnO')) && any(strcmp(headers,'Mn2O3')) %if MnO and Mn2O3 are both included as inputs 

     %If MnO and Mn2O3 are both given, Mn2O3 is converted to MnO and
     %combined with MnO
     D(:,7)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3_mw)+data(:,strcmp(headers,'MnO')));
else
    if any(strcmp(headers,'MnO')) %if MnO is MnO total

        D(:,7)=data(:,strcmp(headers,'MnO'));

    elseif any(strcmp(headers,'Mn2O3'))%if Mn2O3 total is given, converted to MnO total

        %calculate moles of cations
        D(:,7)=(data(:,strcmp(headers,'Mn2O3')).*((2*MnO_mw)./Mn2O3));
    else
        %Mn2O3 and MnO are not included, Mn2 = 0
    end
end


D(:,8)=data(:,strcmp(headers,'MgO')); %for MgO
D(:,9)=data(:,strcmp(headers,'CaO')); %for CaO
D(:,10)=data(:,strcmp(headers,'Na2O')); %for Na2O
D(:,11)=data(:,strcmp(headers,'K2O')); %for K2O

%calculates for H2O if it is included in the analysis
if any(strcmp(headers,'H2O'))
    D(:,12)=data(:,strcmp(headers,'H2O')); %for H2O
end 

%calculates for F if it is included in the analysis
if any(strcmp(headers,'F'))
    D(:,13)=data(:,strcmp(headers,'F')); %for F
end 

%calculates for Cl if it is included in the analysis
if any(strcmp(headers,'Cl'))
    D(:,14)=data(:,strcmp(headers,'Cl')); %for Cl
end 


%% Initial Fe3+ estimation

%calculate the initial apfu of cations and anions 

[O2_Ti,apfu_I,MC]=Amph_O2Ti(D,W,options.TiOH_correction.Value,options.H2Oknown.Value); %calls the initial amphibole apfu calculation

%Lower Limits on Fe3+:

%Atoms per formula unit (apfu) assuming that all Fe is FeO
apfu_Fe=Amph_Fe(D,W,O2_Ti);

if isfield(options,'want_ferric_Fe3') && options.want_ferric_Fe3.Value==true && (isfield(options,'UseKnownFe3') && options.UseKnownFe3.Value==false )
    %apfu assuming 8 Si cations
    apfu_Si=Amph_Si(apfu_Fe,O2_Ti);
    
    %apfu assuming 16 cations (no vacancies on A)
    apfu_Afull=Amph_Afull(apfu_Fe,O2_Ti);
    
    %apfu assuming Na in A site only
    apfu_NaA=Amph_NaA(apfu_Fe,O2_Ti);
    
    apfu_low=zeros(m,11);
    % Fe3+ minimum estimate: highest of 3 normalization factors
    for c=1:m
        %if the 3 minima criteria have lower normalization factors than the
        %four maximum criteria, then Fe3+ cannot be estimated and FeTotal=Fe2+
        if find(apfu_Fe(c,15:17)>apfu_Fe(c,12:14))
            apfu_low(c,:)=apfu_Fe(c,1:1:11);
            low_check{c,1}='Fe3+ cannot be estimated';
        else %otherwise,a Fe3+ lower estimate can be made
            %if the normalization factors are all greater than 1, then all Fe is assumed to be FeO
            if (apfu_Fe(c,12) > 1) && (apfu_Fe(c,13) > 1) && (apfu_Fe(c,14) > 1)
                apfu_low(c,:)=apfu_Fe(c,1:1:11);
                low_check{c,1}={'All Fe is FeO'};
            else
                %If (8/Si) < (16/sum(all cations)) &  < (15/sum(cations Si to Ca))
                %, then the Si normalized formula is correct
                if (apfu_Fe(c,12) < apfu_Fe(c,13)) && (apfu_Fe(c,12) < apfu_Fe(c,14))
                    apfu_low(c,:)=apfu_Si(c,1:1:11);
                    low_check{c,1}={'Criterion 1-1: 8 Si cations on T'};
                else %otherwise, consider the reamining two options
                    
                    %If 16/sum(all cations) < (8/Si) & < (15/sum(cations Si to Ca))
                    %, then the formula normalized to 16 cations is correct
                    if (apfu_Fe(c,13) < apfu_Fe(c,12)) && (apfu_Fe(c,13) < apfu_Fe(c,14))
                        apfu_low(c,:)=apfu_Afull(c,1:1:11);
                        low_check{c,1}={'Criterion 1-2: 16 total cations (no vac. on A)'};
                    else
                        %otherwise, 15/sum(cations Si to Ca) is the only remaining
                        %option for a normalization factor
                        apfu_low(c,:)=apfu_NaA(c,1:1:11);
                        low_check{c,1}={'Criterion 1-3: Na on A site only'};
                        
                    end
                end
            end
        end
    end
    
    
    % Upper Limits on Fe3+
    %apfu assuming all Fe is Fe2O3
    apfu_Fe2O3=Amph_Fe2O3(D,W,O2_Ti);
    O2_Nfact=apfu_Fe2O3(:,19)./apfu_Fe(:,19); %Normalization factor for all Fe as Fe2O3 relative to all Fe as FeO
    
    %apfu assuming Al + Si in T site only
    apfu_SiAlT=Amph_SiAlT(apfu_Fe,O2_Ti);
    
    %apfu assuming K only on the A site
    apfu_KA=Amph_KA(apfu_Fe,O2_Ti);
    
    %apfu assuming fully occupied T & B sites
    apfu_CaNaB=Amph_CaNaB(apfu_Fe,O2_Ti);
    
    %apfu assuming tetrahedral sites completely filled by 3+ and 4+ cations
    apfu_10Fe3=Amph_10Fe3(apfu_Fe,O2_Ti);
    
    
    %Fe3+ maximum estimate: highest of 4 normalization factors
    apfu_hi=zeros(m,11);
    for c=1:m
        
        %if the 3 minima criteria have lower normalization factors than the
        %four maximum criteria, then Fe3+ cannot be estimated and FeTotal=Fe2+
        if find(apfu_Fe(c,15:17)>apfu_Fe(c,12:14))
            apfu_hi(c,:)=apfu_Fe(c,1:1:11);
            hi_check{c,1}={'Fe3+ cannot be estimated'};
        else
            %If the all Fe3+ normalization factor is grater than the other 
            if (O2_Nfact(c,1) > apfu_Fe(c,15)) && (O2_Nfact(c,1) > apfu_Fe(c,16)) && (O2_Nfact(c,1) > apfu_Fe(c,17)) && (O2_Nfact(c,1) > apfu_Fe(c,18))
                apfu_hi(c,:)=apfu_Fe2O3(c,1:1:11);
                hi_check{c,1}={'Criterion 2-5: All Fe is Fe2O3'};
            else
                if (apfu_Fe(c,15) > O2_Nfact(c,1)) && (apfu_Fe(c,15) > apfu_Fe(c,16)) && (apfu_Fe(c,15) > apfu_Fe(c,17)) && (apfu_Fe(c,15) > apfu_Fe(c,18))
                    apfu_hi(c,:)=apfu_SiAlT(c,1:1:11);
                    hi_check{c,1}={'Criterion 2-1: Al + Si normalized to 8 in T site'};
                else
                    if (apfu_Fe(c,16) > O2_Nfact(c,1)) && (apfu_Fe(c,16) > apfu_Fe(c,15)) && (apfu_Fe(c,16) > apfu_Fe(c,17)) && (apfu_Fe(c,16) > apfu_Fe(c,18))
                        apfu_hi(c,:)=apfu_KA(c,1:1:11);
                        hi_check{c,1}={'Criterion 2-2: K only on the A site, fully occupied T, B, & C sites'};
                    else
                        if (apfu_Fe(c,17) > O2_Nfact(c,1)) && (apfu_Fe(c,17) > apfu_Fe(c,15)) && (apfu_Fe(c,17) > apfu_Fe(c,16)) && (apfu_Fe(c,17) > apfu_Fe(c,18))
                            apfu_hi(c,:)=apfu_CaNaB(c,1:1:11);
                            hi_check{c,1}={'Criterion 2-3: Fully occupied T & B sites'};
                        else
                            apfu_hi(c,:)=apfu_10Fe3(c,1:1:11);
                            hi_check{c,1}={'Criterion 2-4: tetrahedral sites completely filled by 3+ and 4+ cations'};
                        end
                    end
                end
            end
        end
        
    end
    
    % calculate the formula with Fe3+/Fetotal determination from high and low
    % limits

    for c=1:m
        %checks if Fe is present
        if apfu_I(c,5) + apfu_I(c,6) > 0     
            apfu(c,:)=(apfu_hi(c,:)+apfu_low(c,:))./2; 
        else
         %if no Fe is present, then the Fe2 only apfu is automatically
         %chosen to prevent incorrect normalizations being chosen for some
         %compositions
            apfu(c,:)=apfu_Fe(c,1:11);
        end
    end
    
    %Adds OH, scales for subsequent itterations
    for c=1:m
        if apfu_I(c,12) > 0
            apfu(c,12) = apfu_I(c,12); %OH from the normalized moles cations
        else
            apfu(c,12) = MC(c,12).*(MC(c,1)./apfu(c,1)); %OH is scaled based on the change in Si content
        end
    end
    
    %Adds F
    for c=1:m
        if apfu_I(c,13) > 0
            apfu(c,13) = apfu_I(c,13); %F from the normalized moles cations
        else
            apfu(c,13) = MC(c,13).*(MC(c,1)./apfu(c,1)); %F is scaled based on the change in Si content
        end
    end
    
    %Adds Cl
    for c=1:m
        if apfu_I(c,14) > 0
            apfu(c,14) = apfu_I(c,14); %Cl from the normalized moles cations
        else
            apfu(c,14) = MC(c,14).*(MC(c,1)./apfu(c,1)); %Cl is scaled based on the change in Si content
        end
    end
    

    %Checks if H2O is measured (given as a column in the input data)
    if options.H2Oknown.Value == 'n'

        for c=1:m
            if (apfu(c,12)+apfu(c,13)+apfu(c,14)) > 2 %if (OH_initial + F + Cl) > 2
                if (apfu(c,13)+apfu(c,14)) < 2 % if (F + Cl) < 2
                    OH_1(c,1)=2-(apfu(c,13)+apfu(c,14)); %initial OH guess is 2-(F+Cl)
                else
                    OH_1(c,1)=0; %W site is filled by Cl and F
                end
            else
                if (apfu(c,12)+apfu(c,13)+apfu(c,14)) < 2 %if (OH_initial + F + Cl) < 2
                    if (2-apfu(c,13)-apfu(c,14)) > 0 %if (2-F-Cl) < 2
                        OH_1(c,1)=(2-apfu(c,13)-apfu(c,14)); % OH = 2-F-Cl
                    else
                        OH_1(c,1)=0; %W site is filled by Cl and F
                    end
                else
                    OH_1(c,1)=apfu(c,12);
                end
            end
        end

        %OH estimation, Step 2: Adjustment for OH=2-2Ti

        %some Ti goes into the T site instead of the M site
        for c=1:m
            %if Ti > (8-(Si+Al)) & (8-(Si+Al))>0
            if apfu(c,2)>(8-(apfu(c,1)+apfu(c,3))) && (8-(apfu(c,1)+apfu(c,3))) > 0
                Ticorr(c,1)=apfu(c,2)-(8-(apfu(c,1)+apfu(c,3))); %Ti in M = Titotal - Ti in T
            else
                Ticorr(c,1)=apfu(c,2); %otherwise all Ti is in M
            end
        end

        % Determine appropriate OH correction
        for c=1:m
            if options.TiOH_correction.Value==true
                if (apfu(c,12) + apfu(c,13) + apfu(c,14)) > 2 %If (OH + F + Cl) > 2
                    if ((2.*apfu(c,12))./((apfu(c,12)+ apfu(c,13) + apfu(c,14))-2.*apfu(c,2)))>0 %if ((2*OH)/((OH+F+Cl)-2*Ti)) >0
                        apfu(c,12)=((2.*apfu(c,12))./((apfu(c,12)+ apfu(c,13) + apfu(c,14))-2.*apfu(c,2)));
                    else
                        apfu(c,12)=0; %W is filled by F, Cl, and O
                    end
                else
                    if (apfu(c,12) + apfu(c,13) + apfu(c,14)) < 2 %If (OH + F + Cl) < 2
                        if ((2-apfu(c,13)-apfu(c,14))-2.*Ticorr(c,1))>0 %if 2-F-Cl-2*Ti > 0
                            apfu(c,12)=((2-apfu(c,13)-apfu(c,14))-2.*Ticorr(c,1)); %OH=(2-F-Cl-2*Ti)
                        else
                            apfu(c,12)=0; %W is filled by F, Cl, and O
                        end
                    else
                        apfu(c,12)=apfu(c,12); %OH content is not adjusted for Ti
                    end
                end
            else
                apfu(c,12)=OH_1(c,1);
            end
        end
        
    else
        %If H2O is in the input data nothing happens
    end
    
    %Moles of anions
    apfu_an(:,1)=apfu(:,1).*2; %SiO2
    apfu_an(:,2)=apfu(:,2).*2; %TiO2
    apfu_an(:,3)=apfu(:,3).*1.5; %Al2O3
    apfu_an(:,4)=apfu(:,4).*1.5; %Cr2O3
    apfu_an(:,5)=apfu(:,5).*1.5; %Fe2O3
    apfu_an(:,6)=apfu(:,6); %FeO
    apfu_an(:,7)=apfu(:,7); %MnO
    apfu_an(:,8)=apfu(:,8); %MgO
    apfu_an(:,9)=apfu(:,9); %CaO
    apfu_an(:,10)=apfu(:,10).*0.5; %Na2O
    apfu_an(:,11)=apfu(:,11).*0.5; %K2O
    apfu_an(:,12)=apfu(:,12).*0.5; %H2O
    
    %for F
    for c=1:m
        if (apfu(c,13)+apfu(c,14))>2 %F + Cl cannot > 2
            apfu_an(c,13)=2.*(apfu(c,13)./(apfu(c,13)+apfu(c,14))); %scales F
        else
            apfu_an(c,13)=apfu(c,13); %Does not scale F
        end
    end
    
    %for Cl
    for c=1:m
        if (apfu(c,13)+apfu(c,14))>2 %F + Cl cannot > 2
            apfu_an(c,14)=2.*(apfu(c,14)./(apfu(c,13)+apfu(c,14))); %scales Cl
        else
            apfu_an(c,14)=apfu(c,14); %Does not scale Cl
        end
    end
    
    
    %the anion sum may be > 24
    anion_sum=sum(apfu_an(:,1:12),2)+0.5*apfu_an(:,13)+0.5*apfu_an(:,14);
    anion_norm=24./anion_sum;
    
    apfu_n=apfu.*anion_norm; %cations normalized to 24 anions
    apfu_ann=apfu_an.*anion_norm; %anions normalized to 24 anions
    Fe3_ratio(:,1)=apfu(:,5)./(apfu(:,5)+apfu(:,6));
    Fe3_ratio(isnan(Fe3_ratio),1)=0;

    molar_mass=apfu_n(:,1).*SiO2_mw+apfu_n(:,2).*TiO2_mw+(apfu_n(:,3)./2).*Al2O3_mw+(apfu_n(:,4)./2).*Cr2O3_mw+(apfu_n(:,5)./2).*Fe2O3_mw+apfu_n(:,6).*FeO_mw+apfu_n(:,7).*MnO_mw+apfu_n(:,8).*MgO_mw+apfu_n(:,9).*CaO_mw+(apfu_n(:,10)./2).*Na2O_mw+(apfu_n(:,11)./2).*K2O_mw+(apfu_n(:,12)./2).*H2O_mw+apfu_ann(:,13).*F_mw+apfu_ann(:,14).*Cl_mw-apfu_ann(:,13).*15.9994.*0.5-apfu_ann(:,14).*15.9994.*0.5;
    
    %new wt% of oxides
    D2=D;
    D2(:,5)=Fe3_ratio(:,1).*(D(:,5).*(2*FeO_mw/Fe2O3_mw)+D(:,6))./(2*FeO_mw/Fe2O3_mw);
    D2(:,6)=(1-Fe3_ratio(:,1)).*(D(:,5).*(2*FeO_mw/Fe2O3_mw)+D(:,6));

    %Checks if H2O is measured (given as a column in the input data)
    if options.H2Oknown.Value == 'n'
        D2(:,12)=apfu_n(:,12).*H2O_mw.*0.5.*(1./molar_mass)*100;
    end
    
    D2(D2<1e-3) = 0;


    %% Itterates the calculation 10 times
    for z=1:25
        %calculate the initial apfu of cations and anions
        
        [O2_Ti2,apfu_I2,MC2]=Amph_O2Ti(D2,W,options.TiOH_correction,options.H2Oknown.Value); %calls the initial amphibole apfu calculation
        
        %Lower Limits on Fe3+:
        
        %Atoms per formula unit (apfu) assuming that all Fe is FeO
        apfu_Fe2=Amph_Fe(D2,W,O2_Ti2);
        
        %apfu assuming 8 Si cations
        apfu_Si2=Amph_Si(apfu_Fe2,O2_Ti2);
        
        %apfu assuming 16 cations (no vacancies on A)
        apfu_Afull2=Amph_Afull(apfu_Fe2,O2_Ti2);
        
        %apfu assuming Na in A site only
        apfu_NaA2=Amph_NaA(apfu_Fe2,O2_Ti2);
        
        apfu_low2=zeros(m,11);
        % Fe3+ minimum estimate: highest of 3 normalization factors
        for c=1:m
            %if the 3 minima criteria have lower normalization factors than the
            %four maximum criteria, then Fe3+ cannot be estimated and FeTotal=Fe2+
            if find(apfu_Fe(c,15:17)>apfu_Fe(c,12:14))
                apfu_low2(c,:)=apfu_Fe2(c,1:1:11);
                low_check2{c,1}='Fe3+ cannot be estimated';
            else %otherwise,a Fe3+ lower estimate can be made
                %if the normalization factors are all greater than 1, then all Fe is assumed to be FeO
                if (apfu_Fe2(c,12) > 1) && (apfu_Fe2(c,13) > 1) && (apfu_Fe2(c,14) > 1)
                    apfu_low2(c,:)=apfu_Fe2(c,1:1:11);
                    low_check2{c,1}={'All Fe is FeO'};
                else
                    %If (8/Si) < (16/sum(all cations)) &  < (15/sum(cations Si to Ca))
                    %, then the Si normalized formula is correct
                    if (apfu_Fe2(c,12) < apfu_Fe2(c,13)) && (apfu_Fe2(c,12) < apfu_Fe2(c,14))
                        apfu_low2(c,:)=apfu_Si2(c,1:1:11);
                        low_check2{c,1}={'Criterion 1-1: 8 Si cations on T'};
                    else %otherwise, consider the reamining two options
                        
                        %If 16/sum(all cations) < (8/Si) & < (15/sum(cations Si to Ca))
                        %, then the formula normalized to 16 cations is correct
                        if (apfu_Fe2(c,13) < apfu_Fe2(c,12)) && (apfu_Fe2(c,13) < apfu_Fe2(c,14))
                            apfu_low2(c,:)=apfu_Afull2(c,1:1:11);
                            low_check2{c,1}={'Criterion 1-2: 16 total cations (no vac. on A)'};
                        else
                            %otherwise, 15/sum(cations Si to Ca) is the only remaining
                            %option for a normalization factor
                            apfu_low2(c,:)=apfu_NaA2(c,1:1:11);
                            low_check2{c,1}={'Criterion 1-3: Na on A site only'};
                            
                        end
                    end
                end
            end
        end
        
        
        % Upper Limits on Fe3+
        %apfu assuming all Fe is Fe2O3
        apfu_Fe2O32=Amph_Fe2O3(D2,W,O2_Ti2);
        O2_Nfact2=apfu_Fe2O32(:,19)./apfu_Fe2(:,19); %Normalization factor for all Fe as Fe2O3 relative to all Fe as FeO
        
        %apfu assuming Al + Si in T site only
        apfu_SiAlT2=Amph_SiAlT(apfu_Fe2,O2_Ti2);
        
        %apfu assuming K only on the A site
        apfu_KA2=Amph_KA(apfu_Fe2,O2_Ti2);
        
        %apfu assuming fully occupied T & B sites
        apfu_CaNaB2=Amph_CaNaB(apfu_Fe2,O2_Ti2);
        
        %apfu assuming tetrahedral sites completely filled by 3+ and 4+ cations
        apfu_10Fe32=Amph_10Fe3(apfu_Fe2,O2_Ti2);
        
        
        %Fe3+ maximum estimate: highest of 4 normalization factors
        apfu_hi2=zeros(m,11);
        for c=1:m
            
            %if the 3 minima criteria have lower normalization factors than the
            %four maximum criteria, then Fe3+ cannot be estimated and FeTotal=Fe2+
            if find(apfu_Fe(c,15:17)>apfu_Fe(c,12:14))
                apfu_hi2(c,:)=apfu_Fe2(c,1:1:11);
                hi_check2{c,1}={'Fe3+ cannot be estimated'};
            else
                %If the all Fe3+ normalization factor is grater than the other th
                if (O2_Nfact2(c,1) > apfu_Fe2(c,15)) && (O2_Nfact2(c,1) > apfu_Fe2(c,16)) && (O2_Nfact2(c,1) > apfu_Fe2(c,17)) && (O2_Nfact2(c,1) > apfu_Fe2(c,18))
                    apfu_hi2(c,:)=apfu_Fe2O32(c,1:1:11);
                    hi_check2{c,1}={'Criterion 2-5: All Fe is Fe2O3'};
                else
                    if (apfu_Fe2(c,15) > O2_Nfact2(c,1)) && (apfu_Fe2(c,15) > apfu_Fe2(c,16)) && (apfu_Fe2(c,15) > apfu_Fe2(c,17)) && (apfu_Fe2(c,15) > apfu_Fe2(c,18))
                        apfu_hi2(c,:)=apfu_SiAlT2(c,1:1:11);
                        hi_check2{c,1}={'Criterion 2-1: Al + Si normalized to 8 in T site'};
                    else
                        if (apfu_Fe2(c,16) > O2_Nfact2(c,1)) && (apfu_Fe2(c,16) > apfu_Fe2(c,15)) && (apfu_Fe2(c,16) > apfu_Fe2(c,17)) && (apfu_Fe2(c,16) > apfu_Fe2(c,18))
                            apfu_hi2(c,:)=apfu_KA2(c,1:1:11);
                            hi_check2{c,1}={'Criterion 2-2: K only on the A site, fully occupied T, B, & C sites'};
                        else
                            if (apfu_Fe2(c,17) > O2_Nfact2(c,1)) && (apfu_Fe2(c,17) > apfu_Fe2(c,15)) && (apfu_Fe2(c,17) > apfu_Fe2(c,16)) && (apfu_Fe2(c,17) > apfu_Fe2(c,18))
                                apfu_hi2(c,:)=apfu_CaNaB2(c,1:1:11);
                                hi_check2{c,1}={'Criterion 2-3: Fully occupied T & B sites'};
                            else
                                apfu_hi2(c,:)=apfu_10Fe32(c,1:1:11);
                                hi_check2{c,1}={'Criterion 2-4: tetrahedral sites completely filled by 3+ and 4+ cations'};
                            end
                        end
                    end
                end
            end
        end
        
        % calculate the formula with Fe3+/Fetotal determination from high and low
        % limits

        for c=1:m
            %checks if Fe is present
            if apfu_I(c,5) + apfu_I(c,6) > 0
                apfu2(c,1:11)=(apfu_hi2(c,:)+apfu_low2(c,:))./2;
            else
                %if no Fe is present, then the Fe2 only apfu is automatically
                %chosen to prevent incorrect normalizations being chosen for some
                %compositions
                apfu2(c,1:11)=apfu_Fe2(c,1:11);
            end
        end
        
        %Adds OH
        for c=1:m
            if apfu_I2(c,12) > 0
                apfu2(c,12) = apfu_I2(c,12); %OH from the normalized moles cations
            else
                apfu2(c,12) = MC2(c,12).*(MC2(c,1)./apfu2(c,1)); %OH is scaled based on the change in Si content
            end
        end
        
        %Adds F
        for c=1:m
            if apfu_I2(c,13) > 0
                apfu2(c,13) = apfu_I2(c,13); %F from the normalized moles cations
            else
                apfu2(c,13) = MC2(c,13).*(MC2(c,1)./apfu2(c,1)); %F is scaled based on the change in Si content
            end
        end
        
        %Adds Cl
        for c=1:m
            if apfu_I2(c,14) > 0
                apfu2(c,14) = apfu_I2(c,14); %Cl from the normalized moles cations
            else
                apfu2(c,14) = MC2(c,14).*(MC2(c,1)./apfu2(c,1)); %Cl is scaled based on the change in Si content
            end
        end

        %Checks if H2O is measured (given as a column in the input data)
        if options.H2Oknown.Value == 'n'

            %initial OH estimate if H2O is not measured
            for c=1:m
                if (apfu2(c,12)+apfu2(c,13)+apfu2(c,14)) > 2 %if (OH_initial + F + Cl) > 2
                    if (apfu2(c,13)+apfu2(c,14)) < 2 % if (F + Cl) < 2
                        OH_2(c,1)=2-(apfu2(c,13)+apfu2(c,14)); %initial OH guess is 2-(F+Cl)
                    else
                        OH_2(c,1)=0; %W site is filled by Cl and F
                    end
                else
                    if (apfu2(c,12)+apfu2(c,13)+apfu2(c,14)) < 2 %if (OH_initial + F + Cl) < 2
                        if (2-apfu2(c,13)-apfu2(c,14)) > 0 %if (2-F-Cl) < 2
                            OH_2(c,1)=(2-apfu2(c,13)-apfu2(c,14)); % OH = 2-F-Cl
                        else
                            OH_2(c,1)=0; %W site is filled by Cl and F
                        end
                    else
                        OH_2(c,1)=apfu2(c,12);
                    end
                end
            end

            %OH estimation, Step 2: Adjustment for OH=2-2Ti

            %some Ti goes into the T site instead of the M site
            for c=1:m
                %if Ti > (8-(Si+Al)) & (8-(Si+Al))>0
                if apfu2(c,2)>(8-(apfu2(c,1)+apfu2(c,3))) && (8-(apfu2(c,1)+apfu2(c,3))) > 0
                    Ticorr2(c,1)=apfu2(c,2)-(8-(apfu2(c,1)+apfu2(c,3))); %Ti in M = Titotal - Ti in T
                else
                    Ticorr2(c,1)=apfu2(c,2); %otherwise all Ti is in M
                end
            end

            %Determine appropriate OH correction
            for c=1:m
                if  options.TiOH_correction.Value==true
                    if (apfu2(c,12) + apfu2(c,13) + apfu2(c,14)) > 2 %If (OH + F + Cl) > 2
                        if ((2.*apfu2(c,12))./((apfu2(c,12)+ apfu2(c,13) + apfu2(c,14))-2.*apfu2(c,2)))>0 %if ((2*OH)/((OH+F+Cl)-2*Ti)) >0
                            apfu2(c,12)=((2.*apfu2(c,12))./((apfu2(c,12)+ apfu2(c,13) + apfu2(c,14))-2.*apfu2(c,2)));
                        else
                            apfu2(c,12)=0; %W is filled by F, Cl, and O
                        end
                    else
                        if (apfu2(c,12) + apfu2(c,13) + apfu2(c,14)) < 2 %If (OH + F + Cl) < 2
                            if ((2-apfu2(c,13)-apfu2(c,14))-2.*Ticorr2(c,1))>0 %if 2-F-Cl-2*Ti > 0
                                apfu2(c,12)=((2-apfu2(c,13)-apfu2(c,14))-2.*Ticorr2(c,1)); %OH=(2-F-Cl-2*Ti)
                            else
                                apfu2(c,12)=0; %W is filled by F, Cl, and O
                            end
                        else
                            apfu2(c,12)=apfu2(c,12); %OH content is not adjusted for Ti
                        end
                    end
                else
                    apfu2(c,12)=OH_2(c,1);
                end
            end

        else
            %if H2O is measured nothing happens
        end

        %Moles of anions
        apfu_an2(:,1)=apfu2(:,1)*2; %SiO2
        apfu_an2(:,2)=apfu2(:,2)*2; %TiO2
        apfu_an2(:,3)=apfu2(:,3)*1.5; %Al2O3
        apfu_an2(:,4)=apfu2(:,4)*1.5; %Cr2O3
        apfu_an2(:,5)=apfu2(:,5)*1.5; %Fe2O3
        apfu_an2(:,6)=apfu2(:,6); %FeO
        apfu_an2(:,7)=apfu2(:,7); %MnO
        apfu_an2(:,8)=apfu2(:,8); %MgO
        apfu_an2(:,9)=apfu2(:,9); %CaO
        apfu_an2(:,10)=apfu2(:,10)*0.5; %Na2O
        apfu_an2(:,11)=apfu2(:,11)*0.5; %K2O
        apfu_an2(:,12)=apfu2(:,12)*0.5; %H2O
        
        %for F
        for c=1:m
            if (apfu2(c,13)+apfu2(c,14))>2 %F + Cl cannot > 2
                apfu_an2(c,13)=2.*(apfu2(c,13)./(apfu2(c,13)+apfu2(c,14))); %scales F
            else
                apfu_an2(c,13)=apfu2(c,13); %Does not scale F
            end
        end
        
        %for Cl
        for c=1:m
            if (apfu2(c,13)+apfu2(c,14))>2 %F + Cl cannot > 2
                apfu_an2(c,14)=2.*(apfu2(c,14)./(apfu2(c,13)+apfu2(c,14))); %scales Cl
            else
                apfu_an2(c,14)=apfu2(c,14); %Does not scale Cl
            end
        end
        
        %the anionn sum may be > 24
        anion_sum2=sum(apfu_an2(:,1:12),2)+0.5*apfu_an2(:,13)+0.5*apfu_an2(:,14);
        anion_norm2=24./anion_sum2;
        
        apfu_n2=apfu2.*anion_norm2; %cations normalized to 24 anions
        apfu_ann2=apfu_an2.*anion_norm2; %anions normalized to 24 anions
        Fe3_ratio2(:,1)=apfu_n2(:,5)./(apfu_n2(:,5)+apfu_n2(:,6));
        Fe3_ratio2(isnan(Fe3_ratio2),1)=0;

        
        molar_mass2=apfu_n2(:,1).*SiO2_mw+apfu_n2(:,2).*TiO2_mw+(apfu_n2(:,3)./2).*Al2O3_mw+(apfu_n2(:,4)./2).*Cr2O3_mw+(apfu_n2(:,5)./2).*Fe2O3_mw+apfu_n2(:,6).*FeO_mw+apfu_n2(:,7).*MnO_mw+apfu_n2(:,8).*MgO_mw+apfu_n2(:,9).*CaO_mw+(apfu_n2(:,10)./2).*Na2O_mw+(apfu_n2(:,11)./2).*K2O_mw+(apfu_n2(:,12)./2).*H2O_mw+apfu_ann2(:,13).*F_mw+apfu_ann2(:,14).*Cl_mw-apfu_ann2(:,13).*15.9994.*0.5-apfu_ann2(:,14).*15.9994.*0.5;
        
        %new wt% of oxides
        D2=D;
        D2(:,5)=Fe3_ratio(:,1).*(D(:,5).*(2*FeO_mw/Fe2O3_mw)+D(:,6))./(2*FeO_mw/Fe2O3_mw);
        D2(:,6)=(1-Fe3_ratio(:,1)).*(D(:,5).*(2*FeO_mw/Fe2O3_mw)+D(:,6));

        %Checks if H2O is measured (given as a column in the input data)
        if options.H2Oknown.Value == 'n'
            D2(:,12)=apfu_n(:,12).*H2O_mw.*0.5.*(1./molar_mass)*100;
        end

        D2(D2<1e-3) = 0;

    end
    
    %% Calculate structural formula
    StrctFrm=StrctFrm_amph_fun(apfu_n2);
    
    %% Data Plotting
    
        
        if isfield(options,'plot_state') && options.plot_state==true
            
            Amp_Plot(:,1)=Fe3_ratio2(:,1); %Fe3+/Fetotal
            Amp_Plot(:,2)=(apfu_n2(:,5)+apfu_n2(:,6)); %Fetotal
            Amp_Plot(:,3)=StrctFrm(:,15)./(StrctFrm(:,15)+StrctFrm(:,16)); %Ca/(Ca+Na) in B
            Amp_Plot(:,4)=StrctFrm(:,19)+StrctFrm(:,20)+2.*StrctFrm(:,18); %A(Na + K + 2Ca), Ca=0
            Amp_Plot(:,5)=StrctFrm(:,4)+StrctFrm(:,7)+2.*StrctFrm(:,5); %C(Al + Fe3+ +2Ti)
            Amp_Plot(:,6)=StrctFrm(:,8)./(StrctFrm(:,8)+StrctFrm(:,9)+StrctFrm(:,13)); %XMg
            Amp_Plot(:,7)=StrctFrm(:,9)./(StrctFrm(:,9)+StrctFrm(:,8)+StrctFrm(:,10)); %Fe2+/(Fe2+ + Mg + Mn) in C
            Amp_Plot(:,8)=StrctFrm(:,7)./(StrctFrm(:,7)+StrctFrm(:,4)+StrctFrm(:,5)); %Fe3+/(Fe3+ + Al + Ti) in C
            Amp_Plot(:,9)=StrctFrm(:,1); %Si (T)
            
              %assigns the variable the appropriate symbol marker
        symb=options.plot_settings.symbol;
        %assigns the variable the appropriate fill color
        fil=options.plot_settings.Color;

        symbsize=options.plot_settings.symbol_size;
            
            if  options.want_classplot.Value==true
                
                %Plots for Ca amphiboles
                figure('Name','Calcium Amphiboles');
                xlim([0 2])
                hold on
                ylim([0 1])
                hold on
                plot([0 2],[0.5 0.5],'k','linewidth',0.5)
                hold on
                plot([0.5 0.5],[0.0 1.0],'k','linewidth',0.5)
                hold on
                plot([1.5 1.5],[0.0 1.0],'k','linewidth',0.5)
                hold on
                text(0.25,0.25,'Tremolite','FontSize',12,'HorizontalAlignment','center')
                text(0.25,0.75,'Edenite','FontSize',12,'HorizontalAlignment','center')
                text(1,0.25,'Magnesio-Hornblende','FontSize',12,'HorizontalAlignment','center')
                text(1,0.75,'Pargasite','FontSize',12,'HorizontalAlignment','center')
                text(1.75,0.75,'Sadanagaite','FontSize',12,'HorizontalAlignment','center')
                text(1.75,0.25,'Tschermakite','FontSize',12,'HorizontalAlignment','center')
                xlabel('^{C}(Al + Fe^{3+} + 2Ti) apfu')
                ylabel('^{A}(Na + K + 2Ca) apfu')
                box on
                hold on
                
                for c=1:m
                    if Amp_Plot(c,3) >= 0.75
                        scatter(Amp_Plot(c,5),Amp_Plot(c,4),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
                        hold on
                    end
                end
                
                figure('Name','Calcium Amphiboles 2');
                xlim([5.5 8])
                hold on
                ylim([0 1])
                hold on
                plot([5.5 8],[0.5 0.5],'k','linewidth',0.5)
                hold on
                plot([6.5 6.5],[0.0 1.0],'k','linewidth',0.5)
                hold on
                plot([7.5 7.5],[0.0 1.0],'k','linewidth',0.5)
                hold on
                plot([7.5 8],[0.9 0.9],'k','linewidth',0.5)
                hold on
                text(7.0,0.25,'Ferro-hornblende','FontSize',12,'HorizontalAlignment','center')
                text(6.0,0.25,'Ferro-tschermakite','FontSize',12,'HorizontalAlignment','center')
                text(7.0,0.75,'Magnesio-Hornblende','FontSize',12,'HorizontalAlignment','center')
                text(6.0,0.75,'Tschermakite','FontSize',12,'HorizontalAlignment','center')
                text(7.75,0.25,'Ferro-actinolite','FontSize',12,'HorizontalAlignment','center')
                text(7.75,0.75,'Actinolite','FontSize',12,'HorizontalAlignment','center')
                text(7.75,0.95,'Tremolite','FontSize',12,'HorizontalAlignment','center')
                xlabel('Si apfu')
                ylabel('X_{Mg}')
                box on
                hold on
                
                for c=1:m
                    if Amp_Plot(c,3) >= 0.75
                        scatter(Amp_Plot(c,9),Amp_Plot(c,6),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
                        hold on
                    end
                end
                
                
                %Na-Ca amphiboles
                figure('Name','Sodium-Calcium Amphiboles');
                xlim([0 2])
                hold on
                ylim([0 1])
                hold on
                plot([0 2],[0.5 0.5],'k','linewidth',0.5)
                hold on
                plot([0.5 0],[0 0.5],'k','linewidth',0.5)
                hold on
                plot([0.5 0.5],[0.5 1.0],'k','linewidth',0.5)
                hold on
                plot([1.5 1.5],[0.0 1.0],'k','linewidth',0.5)
                hold on
                text(0.25,0.75,'Richterite','FontSize',12,'HorizontalAlignment','center')
                text(1,0.25,'Winchite','FontSize',12,'HorizontalAlignment','center')
                text(1,0.75,'Katophorite','FontSize',12,'HorizontalAlignment','center')
                text(1.75,0.75,'Taramite','FontSize',12,'HorizontalAlignment','center')
                text(1.75,0.25,'Barroisite','FontSize',12,'HorizontalAlignment','center')
                xlabel('^{C}(Al + Fe^{3+} + 2Ti) apfu')
                ylabel('^{A}(Na + K + 2Ca) apfu')
                box on
                hold on
                
                for c=1:m
                    if Amp_Plot(c,3) < 0.75 && Amp_Plot(c,3) > 0.25
                        scatter(Amp_Plot(c,5),Amp_Plot(c,4),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
                        hold on
                    end
                end
                
                %Na amphiboles
                
                figure('Name','Sodium Amphiboles');
                xlim([0 2])
                hold on
                ylim([0 1])
                hold on
                plot([1 2],[0.5 0.5],'k','linewidth',0.5)
                hold on
                plot([1.5 1.5],[0.5 1.0],'k','linewidth',0.5)
                hold on
                plot([0.5 1.5],[1 0],'k','linewidth',0.5)
                hold on
                text(1.15,0.75,'Eckermannite','FontSize',12,'HorizontalAlignment','center')
                text(1.75,0.75,'Nyboite','FontSize',12,'HorizontalAlignment','center')
                text(1.75,0.25,'Glaucophane','FontSize',12,'HorizontalAlignment','center')
                xlabel('^{C}(Al + Fe^{3+} + 2Ti) apfu')
                ylabel('^{A}(Na + K + 2Ca) apfu')
                box on
                hold on
                
                for c=1:m
                    if Amp_Plot(c,3) <= 0.25
                        scatter(Amp_Plot(c,5),Amp_Plot(c,4),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
                        hold on
                    end
                end
                
                figure('Name','Sodium Amphiboles 2');
                xlim([0 1])
                hold on
                ylim([0 1])
                hold on
                plot([0 1],[0.5 0.5],'k','linewidth',0.5)
                hold on
                plot([0.5 0.5],[0 1],'k','linewidth',0.5)
                hold on
                text(0.25,0.25,'Glaucophane','FontSize',12,'HorizontalAlignment','center')
                text(0.25,0.75,'Ferro-Glaucophane','FontSize',12,'HorizontalAlignment','center')
                text(0.75,0.25,'Magnesio-Riebeckite','FontSize',12,'HorizontalAlignment','center')
                text(0.75,0.75,'Riebeckite','FontSize',12,'HorizontalAlignment','center')
                xlabel('Fe^{3+}/(Fe^{3+} + Al + Ti)')
                ylabel('Fe^{2+}/(Fe^{2+} + Mg + Mn)')
                box on
                hold on
                
                for c=1:m
                    if Amp_Plot(c,3) <= 0.25
                        scatter(Amp_Plot(c,8),Amp_Plot(c,7),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
                        hold on
                    end
                end
            end
            
 
            if options.want_FePlot.Value
                
                figure('Name','Fe3+/Fetotal vs Fetotal');
                xlim([0 5])
                hold on
                ylim([0 1])
                hold on
                xlabel('\SigmaFe (apfu)')
                ylabel('Fe^{3+}/\SigmaFe')
                box on
                hold on
                scatter(Amp_Plot(:,2),Amp_Plot(:,1),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
            end
        end
    
    
    %% Calculates Outputs
    
    %automatic Fe3+ estimation
    
    %W cations
    W_site(:,1)=apfu_n2(:,12); %OH
    W_site(:,2)=apfu_n2(:,13); %F
    W_site(:,3)=apfu_n2(:,14); %Cl
    W_site(:,4)=apfu_n2(:,12)+apfu_n2(:,13)+apfu_n2(:,14); %W sum
    
    all=[StrctFrm W_site];

    %limit on significant digits (eliminates rounding noise)
    all(all<1e-5) = 0;
    apfu_n2(apfu_n2<1e-5) = 0;

    
    StrctFrm=array2table(all,'VariableNames',{'Si_T','Al_T','Sum_T','Al_C','Ti_C','Cr_C','Fe3_C','Mg_C','Fe2_C','Mn_C','Sum_C','Mg_B','Fe2_B','Mn_B','Ca_B','Na_B','Sum_B','Ca_A','Na_A','K_A','Cation_Sum','OH_W','F_W','Cl_W','W_sum'});
    apfu=array2table(apfu_n2,'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K','OH','F','Cl'});
    
    limits=[low_check hi_check];
    
    Fe3_limits=array2table(limits,'VariableNames',{'Low_Fe3_limits','High_Fe3_limits'});
    
    amph_class=[apfu_Fe(:,12:18) O2_Nfact];
    
    Fe3_class=array2table(amph_class, 'VariableNames',{'Crit1_1','Crit1_2','Crit1_3','Crit2_1','Crit2_2','Crit2_3','Crit2_4','Crit2_5'});
    
    %Fe3+ low estimates:
    
    %Fe2+ Only formula
    apfu_FeT=array2table(apfu_Fe(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});
    
    %Normalized to 8Si
    apfu_SiT=array2table(apfu_Si2(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});
    
    %Normalized to [Si + Al + Ti + Fe3+ + Fe2+ + Mn + Mg + Ca + Na + K] = 16
    apfu_AfullT=array2table(apfu_Afull2(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});
    
    %Normalized to [Si + Al + Ti + Fe3+ + Fe2+ + Mn + Mg + Ca] = 15
    apfu_NaAT=array2table(apfu_NaA2(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});
    
    %Fe3+ high estimates:
    %Fe3+ only
    
    apfu_Fe2O3T=array2table(apfu_Fe2O32(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});
    
    %Normalized to [Si + Al] = 8 cations
    apfu_SiAlT=array2table(apfu_SiAlT2(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});
    
    %Normalized to [Si + Al + Ti + Fe + Mn + Mg + Ca + Na] = 15 cations
    apfu_KAT=array2table(apfu_KA2(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});
    
    %Normalized to [Si + Al + Ti + Fe + Mn + Mg + Ca + Na] = 15 cations
    apfu_CaNaBT=array2table(apfu_CaNaB2(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});

else

    if isfield(options,'want_ferric_Fe3') && options.want_ferric_Fe3.Value==true
        [StrctFrm, apfu,options_definition]=amphibole_Fe3known(data,headers,options);
    else
        options.UseKnownFe3.Value=true;
        options.Fe3_ratio.Value=0;
        [StrctFrm, apfu,options_definition]=amphibole_Fe3known(data,headers,options);

    end
    %blank outputs not needed if Fe3 is known or FeO only is assumed

    low_check=zeros(m,1);
    hi_check=zeros(m,1);
   
    limits=[low_check hi_check];
    
    Fe3_limits=array2table(limits,'VariableNames',{'Low_Fe3_limits','High_Fe3_limits'});

    apfu_Fe=zeros(m,18);
    apfu_Si2=zeros(m,11);
    apfu_Afull2=zeros(m,11);
    apfu_NaA2=zeros(m,11);
    apfu_Fe2O32=zeros(m,11);
    apfu_SiAlT2=zeros(m,11);
    apfu_KA2=zeros(m,11);
    apfu_CaNaB2=zeros(m,11);

    O2_Nfact=zeros(m,1);

    amph_class=[apfu_Fe(:,12:18) O2_Nfact];
    
    Fe3_class=array2table(amph_class, 'VariableNames',{'Crit1_1','Crit1_2','Crit1_3','Crit2_1','Crit2_2','Crit2_3','Crit2_4','Crit2_5'});
    apfu_FeT=array2table(apfu_Fe(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});
    apfu_SiT=array2table(apfu_Si2(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});
    apfu_AfullT=array2table(apfu_Afull2(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});
    apfu_NaAT=array2table(apfu_NaA2(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});   
    apfu_Fe2O3T=array2table(apfu_Fe2O32(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});
    apfu_SiAlT=array2table(apfu_SiAlT2(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});
    apfu_KAT=array2table(apfu_KA2(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});
    apfu_CaNaBT=array2table(apfu_CaNaB2(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});
end

end
