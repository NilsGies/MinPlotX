function [StrctFrm,apfu,options_definition]=MinPlotX_fun(data,headers,MineralName,SampleName,options)
apfu=table();
StrctFrm=table();

%SampleName=app.MinPlotXData.SampleNames{app.active_entry};
% data=app.MinPlotXData.Sample(app.active_entry).Input.data; %data converted to array
% headers=app.MinPlotXData.Sample(app.active_entry).Input.headers; %saves a list of the oxides in the header
% MineralName = app.MineralDropDown.Value;
% options.save_state
% options.PathName=app.ExportDirEditField.Value;

%% garnet #check if ok
if strcmp(MineralName, 'garnet_MinPlotX')
    [StrctFrm,apfu,options_definition]=garnet_MinPlotX(data,headers,options);

    %% pyroxene #check if ok
elseif strcmp(MineralName, 'pyroxene_MinPlotX') || strcmp(MineralName, 'clinopyroxene_MinPlotX')|| strcmp(MineralName, 'orthopyroxene_MinPlotX')
    [StrctFrm,apfu,options_definition]=pyroxene_MinPlotX(data,headers,options);

    %% olivine # add plot options
elseif strcmp(MineralName, 'olivine_MinPlotX')
    [StrctFrm,apfu,options_definition]=olivine_MinPlotX(data,headers,options);

    %% feldspar #check if ok
elseif strcmp(MineralName, 'feldspar_MinPlotX')
    [StrctFrm,apfu,options_definition]=feldspar_MinPlotX(data,headers,options); %calls feldspar calculation

    %% mica #check if ok
elseif strcmp(MineralName, 'mica_MinPlotX')
    [StrctFrm,apfu,options_definition]=mica_MinPlotX(data,headers,options); %saves only apfu

    %% amphibole %% ready for JW testing
elseif strcmp(MineralName, 'amphibole_MinPlotX')
        [StrctFrm, apfu, options_definition, apfu_FeT, Fe3_limits, Fe3_class, apfu_SiT, apfu_AfullT, apfu_NaAT, apfu_Fe2O3T, apfu_SiAlT, apfu_KAT, apfu_CaNaBT]=amphibole_MinPlotX(data,headers,options); %#ok %outputs the AFPU only

    %% epidote #check if ok
elseif strcmp(MineralName, 'epidote_MinPlotX')
    [StrctFrm,apfu,options_definition]=epidote_MinPlotX(data,headers,options); %calls epidote calculation
%% allanite #check if ok
elseif strcmp(MineralName, 'allanite_MinPlotX')
    [StrctFrm,apfu,options_definition]=allanite_MinPlotX(data,headers,options); %calls epidote calculation

    %% chlorite #check if ok
elseif strcmp(MineralName, 'chlorite_MinPlotX')
    [StrctFrm,apfu,options_definition]=chlorite_MinPlotX(data,headers,options); %saves only apfu
    %% serpentine #check if ok
elseif strcmp(MineralName, 'serpentine_MinPlotX')
    [StrctFrm,apfu,options_definition]=serpentine_MinPlotX(data,headers,options); %saves only apfu

    %% titanite #check if ok
elseif strcmp(MineralName, 'titanite_MinPlotX')
    [StrctFrm,apfu,options_definition]=titanite_MinPlotX(data,headers,options); %saves only apfu

    %% talc #check if ok
elseif strcmp(MineralName, 'talc_MinPlotX')
    [StrctFrm,apfu,options_definition]=talc_MinPlotX(data,headers,options); %saves only apfu

    %% chloritoid % test
elseif strcmp(MineralName, 'chloritoid_MinPlotX')
    [StrctFrm,apfu,options_definition]=chloritoid_MinPlotX(data,headers,options); %saves only apfu

    %% staurolite % test
elseif strcmp(MineralName, 'staurolite_MinPlotX')
    [StrctFrm,apfu,options_definition]=staurolite_MinPlotX(data,headers,options); %saves only apfu but has StrctFrm in function

    %% spinel 
elseif strcmp(MineralName, 'oxyspinel_MinPlotX') || strcmp(MineralName, 'spinel_MinPlotX')
    [StrctFrm,apfu,options_definition]=spinel_MinPlotX(data,headers,options); %calls spinel calculation

    %% cordierite % test
elseif strcmp(MineralName, 'cordierite_MinPlotX')
    [StrctFrm,apfu,options_definition]=cordierite_MinPlotX(data,headers,options); 

    %% sulfide %% test
elseif strcmp(MineralName,'sulfide_MinPlotX')
    [StrctFrm, apfu,options_definition]=sulfide_MinPlotX(data,headers,options);

    %% ilmenite #check if ok
elseif strcmp(MineralName, 'ilmenite_MinPlotX')
    [StrctFrm,apfu,options_definition]=ilmenite_MinPlotX(data,headers,options); %calls ilmenite calculation

    %% apatite % test
elseif strcmp(MineralName, 'apatite_MinPlotX')
    [StrctFrm,apfu,options_definition]=apatite_MinPlotX(data,headers,options); %outputs the structural formula with cation assignment

    %% scapolite % test
elseif strcmp(MineralName, 'scapolite_MinPlotX')
    [StrctFrm, apfu,options_definition]=scapolite_MinPlotX(data,headers,options);

    %% lawsonite % test
elseif strcmp(MineralName, 'lawsonite_MinPlotX')
    [StrctFrm, apfu,options_definition]=lawsonite_MinPlotX(data,headers,options);

elseif strcmp(MineralName,'MultiMineralData')

elseif strcmp(MineralName,'unknown_MinPlotX')

  [StrctFrm, apfu, options_definition]=unknown_MinPlotX(data,headers,options);

%    ElOxDataDef=ReadDefFiles;

%     if options.cation.Value==true
%         [OutputData,OutputVariables] = SF_CatNorm(data,headers,options.n_cations.Value,ElOxDataDef);
%     elseif options.oxygen.Value==true
%         [OutputData,OutputVariables] = SF_OxNorm(data,headers,options.n_oxygen.Value,ElOxDataDef);
%     elseif options.anion.Value==true
%         [OutputData,OutputVariables] = SF_OxNorm(data,headers,options.n_anions.Value,ElOxDataDef); % ath change function
%     end
% 
%     apfu=array2table(OutputData);
%     if numel(unique(OutputVariables))==numel(OutputVariables)
%         apfu.Properties.VariableNames=OutputVariables; % ath check for Fe Fe / Mn Mn
%     else
%         id=find(strcmp(OutputVariables,'Fe'));
%         OutputVariables(id(1))={'Fe2'};
%         OutputVariables(id(2))={'Fe3'};
%         apfu.Properties.VariableNames=OutputVariables; % ath check for Fe Fe / Mn Mn
%     end
%                StrctFrm=apfu;

elseif endsWith(MineralName,'_XMapT')
    try
        ElOxDataDef=ReadDefFiles;
        ExtFct=char(MineralName);

        [OutputData,OutputVariables] = feval('Core_Call_SF',data,headers,ExtFct,ElOxDataDef); %#ok faster
        %[OutputData,OutputVariables] = Core_Call_SF(data,headers,ExtFct,ElOxDataDef);

        apfu=array2table(OutputData);
        apfu.Properties.VariableNames=OutputVariables;
               StrctFrm=apfu;
 catch ME
        disp('Calculation Failed in MinPlotX_fun - XMapT Fun')
        disp(ME.message)
    end
else
    try
        try
            [StrctFrm, apfu,options_definition] = feval(char(MineralName),data,headers,options);
        catch
            ElOxDataDef=ReadDefFiles;
            ExtFct=char(MineralName);
            [OutputData,OutputVariables] = feval(Core_Call_SF(data,headers,ExtFct,ElOxDataDef));
            %  [data, header] = feval(char(MineralName),data,headers,ElOxDataDef);
            apfu=array2table(OutputData);
            apfu.Properties.VariableNames=OutputVariables;
            StrctFrm=apfu;
        end
    catch ME
        disp('Calculation Failed in MinPlotX_fun - External Fun')
        disp(ME.message)
    end
end

%% Fe3_FeTot

if  istable(apfu)&& istable(StrctFrm)&& any(strcmp('Fe3',apfu.Properties.VariableNames))&& any(strcmp('Fe2',apfu.Properties.VariableNames))
    StrctFrm.Fe3_FeTot=apfu.Fe3./(apfu.Fe3+apfu.Fe2);
end

%% XMg

if   istable(apfu)&& istable(StrctFrm)&&any(strcmp('Mg',apfu.Properties.VariableNames))&& any(strcmp('Fe2',apfu.Properties.VariableNames))
    StrctFrm.XMg(:)=0;
    StrctFrm.XMg(apfu.Mg+apfu.Fe2>0)=apfu.Mg(apfu.Mg+apfu.Fe2>0)./(apfu.Mg(apfu.Mg+apfu.Fe2>0)+apfu.Fe2(apfu.Mg+apfu.Fe2>0));
end

if isfield(options,'save_state') && options.save_state==true
    if istable(StrctFrm)
        writetable(StrctFrm,fullfile(options.PathName,[SampleName,'_structuralformula.csv']),'Delimiter','\t');
    end

    if istable(apfu)
        writetable(apfu,fullfile(options.PathName,[SampleName,'_apfu.csv']),'Delimiter','\t');
    end
end

if isempty(apfu)
    apfu=table();
end
if isempty(StrctFrm)
    StrctFrm=table();
end



end