function [ax2plot,options_definition,legend_str] = plot_XYZ_data_fun(T,ax2plot,options)

[~,options_definition]=plot_XYZ_fun([],[]);

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

if not(exist('ax2plot','var')) || isempty(ax2plot)
    options.plot_tern.Value=true;
    [ax2plot]=plot_XYZ_fun([],options);
end

switch options.type.Value


    case 'ol_CaO_Al2O3_NiO_XYZ'
        if any(strcmp(T.Properties.VariableNames,'CaO'))
            X1=T.CaO;
        else
            X1=  zeros(size(T.SiO2));
        end

        if any(strcmp(T.Properties.VariableNames,'Al2O3'))
            Y1=  T.Al2O3;
        else
            Y1=  zeros(size(T.SiO2));
        end

        if any(strcmp(T.Properties.VariableNames,'NiO'))
            Z1=  T.NiO;
        else
            Z1=  zeros(size(T.SiO2));
        end

    case 'grt_CaO_Al2O3_MgO_XYZ'
        %         X1=T.CaO;
%         Y1=  T.Al2O3;
%         Z1=  T.MgO;
         X1=T.apfu_Ca;
        Y1=  T.apfu_Al;
        Z1=  T.apfu_Mg;
    case 'mica_xMg_Al_Si'

        %only plots trioctahedral, 50 % rule
        condition=T.StrctFrm_XTriOct >= 0.50;

        XMg=zeros(size(T,1),1);
        XMg(T.StrctFrm_Mg_M>0)=T.StrctFrm_Mg_M(T.StrctFrm_Mg_M>0)./(T.StrctFrm_Mg_M(T.StrctFrm_Mg_M>0)+T.StrctFrm_Fe2_M(T.StrctFrm_Mg_M>0));
        X1=XMg(condition);
        Y1=T.StrctFrm_Al_M(condition);
        Z1=T.apfu_Si(condition);

    case 'custom_XYZ'
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

        if isfield(options,'Z1')&& not(isempty(options.Z1))
            Z1=T.(char(options.Z1));
        else
            Z1=zeros(size(Y1));
        end
end

if isfield(options,'C1')&& not(isempty(options.C1)) && ismember(options.C1,T.Properties.VariableNames)
    C1=T.(char(options.C1));
elseif isfield(options,'C1')&& not(isempty(options.C1)) && strcmp(options.C1,'1:n')
    C1=(1:size(T,1))';
else
    C1=(1:size(T,1))';
    options.ColorData.cbar_label='Colormap Error - Check Selected Variable';
end

if exist('condition','var')
    C1=C1(condition);
end

if  exist('X1','var') &&  exist('Y1','var') &&  exist('Z1','var') &&any(not(any((isnan([X1 Y1 Z1])),2)))

    if isfield(options,'ColorData') && isfield(options.ColorData,'Value') && options.ColorData.Value && isfield(options.ColorData,'colormap') && not(isempty(options.ColorData.colormap))
        colormap(ax2plot,options.ColorData.colormap)

        if  isfield(options.ColorData,'cbar_label')
            cbar_str=options.ColorData.cbar_label;
        else
            cbar_str='';
        end

        sc_obj=  scatter3(ax2plot,X1(:),Y1(:),Z1(:),symbsize,C1(:),symb,'filled','MarkerFaceAlpha',mfa,'MarkerEdgeAlpha',mea,'MarkerEdgeColor',mec,'LineWidth',mlw);
        cbar= colorbar;
        cbar.Label.String=cbar_str;

     
    else
        scatter3(ax2plot,X1(:),Y1(:),Z1(:),symbsize,symb,'filled','MarkerFaceAlpha',mfa,'MarkerEdgeAlpha',mea,'MarkerEdgeColor',mec,'MarkerFaceColor',fil,'LineWidth',mlw);
 
    end
                

    if  isfield(options,'legend') && isfield(options.legend,'String')
        legend_str=options.legend.String(any(not(isnan(X1)) & not(isnan(Y1))));
    end

end

end

