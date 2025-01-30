function [ax2plot,options_definition,legend_str] = plot_XY_data_fua(T,ax2plot,options)

[~,options_definition]=plot_XY_fun([],[]);

if not(exist('T','var'))
    T=[];
end
if not(exist('ax2plot','var'))
    ax2plot=[];
end
if not(exist('options','var'))
    options=[];
end

%% Plot Options definition
if  (isfield(options,'default') && options.default==true) || not(isfield(options,'plot_settings'))|| isempty(options.plot_settings)  || not(isfield(options.plot_settings,'MarkerFaceColor')) || isempty(options.plot_settings.MarkerFaceColor)
    fil=[0 0.4470 0.7410];
    options_definition.plot_settings.MarkerFaceColor=fil;
else
    fil=options.plot_settings.MarkerFaceColor;
end

if  (isfield(options,'default') && options.default==true)|| not(isfield(options,'plot_settings')) || isempty(options.plot_settings)  || not(isfield(options.plot_settings,'MarkerFaceAlpha')) || isempty(options.plot_settings.MarkerFaceAlpha) || not(isfield(options.plot_settings.MarkerFaceAlpha,'Value')) || isempty(options.plot_settings.MarkerFaceAlpha.Value)
    mfa=3/8;
    options_definition.plot_settings.MarkerFaceAlpha=mfa;
else
    mfa=options.plot_settings.MarkerFaceAlpha.Value;
end

if  (isfield(options,'default') && options.default==true) || not(isfield(options,'plot_settings'))|| isempty(options.plot_settings)  || not(isfield(options.plot_settings,'symbol')) || isempty(options.plot_settings.symbol)
    symb='o';
    options_definition.plot_settings.symbol=symb;
else
    symb=options.plot_settings.symbol;
end

if  (isfield(options,'default') && options.default==true) || not(isfield(options,'plot_settings'))|| isempty(options.plot_settings)  || not(isfield(options.plot_settings,'MarkerEdgeColor')) || isempty(options.plot_settings.MarkerEdgeColor)
    if strcmp(symb,'.')
        mec=fil;
    else
        mec=[0 0 0];
    end
    options_definition.plot_settings.MarkerEdgeColor=mec;
else
    mec=options.plot_settings.MarkerEdgeColor;

end

if  (isfield(options,'default') && options.default==true) || not(isfield(options,'plot_settings'))|| isempty(options.plot_settings)  || not(isfield(options.plot_settings,'MarkerEdgeAlpha')) || isempty(options.plot_settings.MarkerEdgeAlpha)  || not(isfield(options.plot_settings.MarkerEdgeAlpha,'Value')) || isempty(options.plot_settings.MarkerEdgeAlpha.Value)
    mea=4/8;
    options_definition.plot_settings.MarkerEdgeAlpha=mea;
else
    mea=options.plot_settings.MarkerEdgeAlpha.Value;
end

if  (isfield(options,'default') && options.default==true) || not(isfield(options,'plot_settings'))|| isempty(options.plot_settings)  || not(isfield(options.plot_settings,'MarkerLineWidth')) || isempty(options.plot_settings.MarkerLineWidth)  || not(isfield(options.plot_settings.MarkerLineWidth,'Value')) || isempty(options.plot_settings.MarkerLineWidth.Value)
    mlw=0.5;
    options_definition.plot_settings.MarkerLineWidth.Value=mlw;
else
    mlw=options.plot_settings.MarkerLineWidth.Value;
end

if  (isfield(options,'default') && options.default==true) || not(isfield(options,'plot_settings'))|| isempty(options.plot_settings)  || not(isfield(options.plot_settings,'symbol_size')) || isempty(options.plot_settings.symbol_size)
    symbsize=100;
    options_definition.plot_settings.symbol_size=symbsize;
else
    symbsize= options.plot_settings.symbol_size;
end

legend_str={};

if not(exist('options','var')) || isempty(options)
    if not(exist('ax2plot','var'))
        ax2plot=[];
    end
    return
end

if (isfield(options,'default') && options.default==true)
    type=options.type.Value;
    options=options_definition;
    options.type.Value=type;
end

if not(exist('options','var')) || isempty(options)
    if not(exist('ax2plot','var'))
        ax2plot=[];
    end
    return
end

if (isfield(options,'default') && options.default==true)
    type=options.type.Value;
    options=options_definition;
    options.type.Value=type;
end

if not(exist('ax2plot','var'))
    options.plot_tern.Value=true;
    [ax2plot]=plot_XY_fun([],options);
end

switch options.type.Value

    case 'cpx_QJ_XY'
        Y1=T.apfu_Ca+T.apfu_Mg+T.apfu_Fe2;
        X1=T.apfu_Na.*2;

    case 'cpx_catscaes'
        if ismember('StrctFrm_Xks',T.Properties.VariableNames) && ismember('StrctFrm_XCaes',T.Properties.VariableNames)
            Y1=T.StrctFrm_Xks;
            X1=T.StrctFrm_XCaes;
        end

    case 'cpx_Kcpxcaes'
        if ismember('StrctFrm_XKkos',T.Properties.VariableNames) && ismember('StrctFrm_XKjd',T.Properties.VariableNames)
            Y1=T.StrctFrm_XKkos + T.StrctFrm_XKjd;
            X1=T.StrctFrm_XCaes;
        end

    case 'ol_NiO_XY'
        X1=T.StrctFrm_XFo.*100;
        if any(strcmp(T.Properties.VariableNames,'NiO'))
            Y1=  T.NiO;
        else
            Y1=  T.SiO2.*0;
        end
    case 'amph_XY_Ca1'
        if ismember('StrctFrm_Ca_B',T.Properties.VariableNames)

            Amp_Ca_Na_B_cond=T.StrctFrm_Ca_B./(T.StrctFrm_Ca_B+T.StrctFrm_Na_B); %Ca/(Ca+Na) in B
            Amp_A_Na_K_2Ca=T.StrctFrm_Na_A+T.StrctFrm_K_A+2.*T.StrctFrm_Ca_A; %A(Na + K + 2Ca), Ca=0
            Amp_C_Al_Fe3_2Ti=T.StrctFrm_Al_C+T.StrctFrm_Fe3_C+2.*T.StrctFrm_Ti_C; %C(Al + Fe3+ +2Ti)

            condition=Amp_Ca_Na_B_cond >= 0.75;

            X1=Amp_C_Al_Fe3_2Ti(condition);
            Y1=Amp_A_Na_K_2Ca(condition);
        end
    case 'amph_XY_Ca2'

        if ismember('StrctFrm_Mg_C',T.Properties.VariableNames)
            Amp_Plot(:,6)=T.StrctFrm_Mg_C./(T.StrctFrm_Mg_C+T.StrctFrm_Fe2_C+T.StrctFrm_Fe2_B); %XMg
            Amp_Plot(:,9)=T.StrctFrm_Si_T; %Si (T)
            Amp_Ca_Na_B_cond=T.StrctFrm_Ca_B./(T.StrctFrm_Ca_B+T.StrctFrm_Na_B); %Ca/(Ca+Na) in B
            condition=Amp_Ca_Na_B_cond >= 0.75;

            X1=Amp_Plot(condition,9);
            Y1=Amp_Plot(condition,6);
        end
    case 'amph_XY_NaCa'
        if ismember('StrctFrm_Ca_B',T.Properties.VariableNames)

            Amp_Ca_Na_B_cond=T.StrctFrm_Ca_B./(T.StrctFrm_Ca_B+T.StrctFrm_Na_B); %Ca/(Ca+Na) in B
            Amp_A_Na_K_2Ca=T.StrctFrm_Na_A+T.StrctFrm_K_A+2.*T.StrctFrm_Ca_A; %A(Na + K + 2Ca), Ca=0
            Amp_C_Al_Fe3_2Ti=T.StrctFrm_Al_C+T.StrctFrm_Fe3_C+2.*T.StrctFrm_Ti_C; %C(Al + Fe3+ +2Ti)

            condition=Amp_Ca_Na_B_cond < 0.75 & Amp_Ca_Na_B_cond > 0.25;

            X1=Amp_C_Al_Fe3_2Ti(condition);
            Y1=Amp_A_Na_K_2Ca(condition);
        end
    case 'amph_XY_Na1'
        if ismember('StrctFrm_Ca_B',T.Properties.VariableNames)

            Amp_Ca_Na_B_cond=T.StrctFrm_Ca_B./(T.StrctFrm_Ca_B+T.StrctFrm_Na_B); %Ca/(Ca+Na) in B
            Amp_A_Na_K_2Ca=T.StrctFrm_Na_A+T.StrctFrm_K_A+2.*T.StrctFrm_Ca_A; %A(Na + K + 2Ca), Ca=0
            Amp_C_Al_Fe3_2Ti=T.StrctFrm_Al_C+T.StrctFrm_Fe3_C+2.*T.StrctFrm_Ti_C; %C(Al + Fe3+ +2Ti)

            condition=Amp_Ca_Na_B_cond <= 0.25;

            X1=Amp_C_Al_Fe3_2Ti(condition);
            Y1=Amp_A_Na_K_2Ca(condition);
        end
    case 'amph_XY_Na2'
        if ismember('StrctFrm_Ca_B',T.Properties.VariableNames)
            Amp_Ca_Na_B_cond=T.StrctFrm_Ca_B./(T.StrctFrm_Ca_B+T.StrctFrm_Na_B); %Ca/(Ca+Na) in B
            Amp_C_Fe2_Mg_Mn=T.StrctFrm_Fe2_C./(T.StrctFrm_Fe2_C+T.StrctFrm_Mg_C+T.StrctFrm_Mn_C); %Fe2+/(Fe2+ + Mg + Mn) in C
            Amp_C_Fe3_Al_Ti=T.StrctFrm_Fe3_C./(T.StrctFrm_Fe3_C+T.StrctFrm_Al_C+T.StrctFrm_Ti_C); %Fe3+/(Fe3+ + Al + Ti) in C

            condition=Amp_Ca_Na_B_cond <= 0.25;
            X1=Amp_C_Fe3_Al_Ti(condition);
            Y1=Amp_C_Fe2_Mg_Mn(condition);
        end
    case 'amph_XY_Fe'
        FeTotal=(T.apfu_Fe3+T.apfu_Fe2); %Fetotal
        Fe3_ratio2=T.apfu_Fe3./(T.apfu_Fe3+T.apfu_Fe2);
        X1=FeTotal;
        Y1=Fe3_ratio2;
    case 'mica_XY_ph'

        %only plots Dioctahedral, 50 % rule
        condition=T.StrctFrm_XDiOct >= 0.50;

        X1=T.apfu_Al(condition);
        Y1=T.apfu_Si(condition);

    case 'mica_XY_SiNa'
        X1=T.apfu_Na;
        Y1=T.apfu_Si;

    case 'mica_XY_NaKFeMg'
        X1=T.apfu_Mg+T.apfu_Fe2;
        Y1=T.apfu_Na./(T.apfu_Na+T.apfu_K);

    case 'mica_XY_bio'

        %only plots trioctahedral, 50 % rule
        condition=T.StrctFrm_XTriOct >= 0.50;

        XMg=zeros(size(T,1),1);
        XMg(T.StrctFrm_Mg_M>0)=T.StrctFrm_Mg_M(T.StrctFrm_Mg_M>0)./(T.StrctFrm_Mg_M(T.StrctFrm_Mg_M>0)+T.StrctFrm_Fe2_M(T.StrctFrm_Mg_M>0));
        X1=XMg(condition);
        Y1=T.StrctFrm_Al_M(condition);

    case 'chl_XY_simp'
        X1=T.StrctFrm_Al_M+T.StrctFrm_Fe3_M+T.StrctFrm_Cr_M;
        Y1=T.StrctFrm_XMg;
    case 'chl_XY_comp'
        X1=T.StrctFrm_Al_M+T.StrctFrm_Fe3_M+T.StrctFrm_Cr_M;;
        Y1=T.StrctFrm_XMg;
    case 'chl_XY_R2Si'
        Y1=T.apfu_Si;
        X1=T.apfu_Fe2+T.apfu_Mn2+T.apfu_Mg+T.apfu_Ni;
    case 'chl_XY_XFeSi'
        X1=T.apfu_Si;
        Y1=T.apfu_Fe3+T.apfu_Fe2;

    case 'aln_XY'
        if ismember('StrctFrm_Al_M',T.Properties.VariableNames)

            X1=T.StrctFrm_Al_M;
            if ismember('StrctFrm_YREE_A',T.Properties.VariableNames)
                Y1=T.StrctFrm_YREE_A+T.StrctFrm_Th_A+T.StrctFrm_U_A;
            else
                Y1=zeros(size(X1));
            end
        end
    case 'ttn_XY_Ti'
        X1=T.StrctFrm_Al_Oct+T.StrctFrm_Fe3_Oct;
        Y1=T.apfu_Ti;

    case 'ttn_XY_F'
        X1=T.StrctFrm_Al_Oct+T.StrctFrm_Fe3_Oct;
        Y1=T.apfu_F;
    case 'crd_XY_NaAlSi'
        Y1=T.apfu_Al+T.apfu_Na;
        X1=T.apfu_Si;
    case 'crd_XY_AlMeNaSi'
        Y1=T.apfu_Al+T.apfu_Fe2+T.apfu_Mn2+T.apfu_Mg;
        X1=T.apfu_Si+T.apfu_Na;
    case 'crd_XY_MevsSiAl'
        Y1=T.apfu_Fe2+T.apfu_Mn2+T.apfu_Mg;
        X1=T.apfu_Si+T.apfu_Al;
    case 'lws_XY_FeTiCrAl'
        Y1=T.apfu_Fe3+T.apfu_Ti+T.apfu_Cr;
        X1=T.apfu_Al;
    case 'custom_XY'
        if isfield(options,'X1')&& not(isempty(options.X1))
            X1=T.(char(options.X1));
        else
            X1=(1:size(T,1))';
        end
        if isfield(options,'Y1')&& not(isempty(options.Y1))
            Y1=T.(char(options.Y1));
        else
            Y1=(1:size(T,1))';
        end

end

if isfield(options,'C1')&& not(isempty(options.C1)) && ismember(options.C1,T.Properties.VariableNames)
    C1=T.(char(options.C1));
elseif isfield(options,'C1')&& not(isempty(options.C1)) && strcmp(options.C1,'1:n')
    C1=1:size(T,1);
else
    C1=1:size(T,1);
    options.ColorData.cbar_label='Colormap Error - Check Selected Variable';
end

if exist('condition','var')
    C1=C1(condition);
end

if  exist('X1','var') &&  exist('Y1','var') &&any(not(any((isnan([X1 Y1])),2)))
    if isfield(options,'ColorData') && isfield(options.ColorData,'Value') && options.ColorData.Value && isfield(options.ColorData,'colormap') && not(isempty(options.ColorData.colormap))
        colormap(ax2plot,options.ColorData.colormap)

        if  isfield(options.ColorData,'cbar_label')
            cbar_str=options.ColorData.cbar_label;
        else
            cbar_str='';
        end

        scatter(ax2plot,X1(:),Y1(:),symbsize,C1(:),symb,'filled','MarkerFaceAlpha',mfa,'MarkerEdgeAlpha',mea,'MarkerEdgeColor',mec,'LineWidth',mlw)
        cbar= colorbar;
        cbar.Label.String=cbar_str;

    else
        if isfield(options,'plot_settings') && isfield(options.plot_settings,'plot_mean_error')&& isfield(options.plot_settings.plot_mean_error,'Value') && options.plot_settings.plot_mean_error.Value==true
            x_plot_error= std(X1)./sqrt(length(X1));
            y_plot_error=std(Y1)./sqrt(length(Y1));
            errorbar(ax2plot,mean(X1),mean(Y1),y_plot_error,y_plot_error,x_plot_error,x_plot_error,'')
            scatter(ax2plot,mean(X1),mean(Y1),symbsize,symb,'filled','MarkerFaceAlpha',mfa,'MarkerEdgeAlpha',mea,'MarkerEdgeColor',mec,'MarkerFaceColor',fil,'LineWidth',mlw)
        else
            scatter(ax2plot,X1(:),Y1(:),symbsize,symb,'filled','MarkerFaceAlpha',mfa,'MarkerEdgeAlpha',mea,'MarkerEdgeColor',mec,'MarkerFaceColor',fil,'LineWidth',mlw)
        end
    end


    if  isfield(options,'legend') && isfield(options.legend,'String')
        legend_str=options.legend.String(any(not(isnan(X1)) & not(isnan(Y1))));
    end
end

end

