function [ax2plot,options_definition,legend_str] = plot_multi_XY_data_fun(T,ax2plot,options)

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
    [ax2plot]=plot_multi_XY_fun([],options);
end

plot_data=true;

switch options.type.Value

    case 'amph_multiplot'
        plot_data=false;
        %%%%% 1
        options_sub=options;
        options_sub.type.Value='amph_XY_Ca1';
                        hold(ax2plot(1),'on');
        [ax2plot(1),~,legend_str_1]=plot_XY_data_fun(T,ax2plot(1),options_sub);
        %%%%% 2
        options_sub.type.Value='amph_XY_Ca2';
        hold(ax2plot(2),'on');
        [ax2plot(2),~,legend_str_2]=plot_XY_data_fun(T,ax2plot(2),options_sub);
        %%%%% 3
        options_sub.type.Value='amph_XY_Na1';
        hold(ax2plot(3),'on');
        [ax2plot(3),~,legend_str_3]=plot_XY_data_fun(T,ax2plot(3),options_sub);
        %%%%% 4
        options_sub.type.Value='amph_XY_NaCa';
        hold(ax2plot(4),'on');
        [ax2plot(4),~,legend_str_4]=plot_XY_data_fun(T,ax2plot(4),options_sub);
        %%%%% 6
        options_sub.type.Value='amph_XY_Na2';
        hold(ax2plot(5),'on');
        [ax2plot(5),~,legend_str_5]=plot_XY_data_fun(T,ax2plot(5),options_sub);
        %%%%% 6
        hold(ax2plot(6),'on');
        options_sub.type.Value='amph_XY_Fe';
        [ax2plot(6),~,legend_str_6]=plot_XY_data_fun(T,ax2plot(6),options_sub);

        legend_str=unique([legend_str_1, legend_str_2, legend_str_3, legend_str_4, legend_str_5, legend_str_6]);

    case 'ol_NiO_Mn_Ca_Mg_XXXY'
        X1(1:3)=  {100.*T.apfu_Mg./(T.apfu_Mg+T.apfu_Fe2)};
       if any(strcmp(T.Properties.VariableNames,'NiO'))
        Y1(1)=  {T.NiO};
       else
        Y1(1)=  {T.SiO2.*0}; 
       end
        Y1(2)=  {T.MnO};
        Y1(3)=  {T.CaO};

    case 'ol_NiO_CaO_Mg_XXXY'
        X1(1:2)=  {100.*T.apfu_Mg./(T.apfu_Mg+T.apfu_Fe2)};
if any(strcmp(T.Properties.VariableNames,'NiO'))
        Y1(1)=  {T.NiO};
       else
        Y1(1)=  {T.SiO2.*0}; 
end
Y1(2)=  {T.CaO};

        % X1=T.StrctFrm_XFo.*100};
    case 'ol_large_XY'

        X1(1:4)=  {100.*T.apfu_Mg./(T.apfu_Mg+T.apfu_Fe2)};
        X1(5:7)=  {T.SiO2};
        X1(8)=   {T.Al2O3};

        Y1(1)=  {T.SiO2};
        Y1(2)=  {T.TiO2};
        Y1(3)=  {T.Cr2O3};
        Y1(4)=  {T.Al2O3};
        Y1(5)=  {T.Na2O};
        Y1(6)=  {T.CaO};
        Y1(7:8)=  {T.K2O};
    case 'ol_test_XY'
        Y1(1)=  {T.FeO./T.MgO};
        X1(1)=  {T.SiO2./T.MnO};

        Y1(2)=   {T.FeO./T.MgO};
        X1(2)=  {T.SiO2./T.NiO};

        Y1(3)=  {T.NiO./T.MnO};
        X1(3)=  {T.MgO./T.FeO};

        Y1(4)=  {T.NiO./T.Cr2O3};
        X1(4)=  {T.MgO./T.FeO};
    case 'ol_logtest1_XY'
        Y1(1)=  {log(T.FeO./T.MgO)};
        X1(1)=  {log(T.SiO2./T.MnO)};

        Y1(2)=   {log(T.FeO./T.MgO)};
        X1(2)=  {log(T.SiO2./T.NiO)};

        Y1(3)=  {log(T.NiO./T.MnO)};
        X1(3)=  {log(T.MgO./T.FeO)};

        Y1(4)=  {log(T.NiO./T.Cr2O3)};
        X1(4)=  {log(T.MgO./T.FeO)};
  %  case 'ol_H2O_large'
% 
% 
%         X1(1:6)= {T.SIMS_H2O};
% 
% 
%         Y1(1)=  {T.StrctFrm_XFo}; %XFo
%         Y1(2)=  {T.SIMS_Ca_corr};
%         Y1(3)=  {T.Cr2O3.*(1/1.4616)*10000};
%         
%         
%         Y1(4)=  {T.SIMS_Al_corr};
%         Y1(5)=  {T.SIMS_Ti_corr};
%         Y1(6)=  {T.SIMS_Li_corr};
case 'ol_H2O_large'


        X1(1:6)= {T.SIMS_H2O};


        Y1(1)=  {T.StrctFrm_XFo}; %XFo
        Y1(2)=  {T.SIMS_Ca_corr};
        Y1(3)=  {T.Cr2O3.*(1/1.4616)*10000};
        
        
        Y1(4)=  {T.SIMS_Al_corr};
        Y1(5)=  {T.SIMS_Ti_corr};
        Y1(6)=  {T.SIMS_Li_corr};

    case 'clinopyroxene_xMg_NaO2_Cr2O3_Al2O3'
        X1(1:3)=  {100.*T.apfu_Mg./(T.apfu_Mg+T.apfu_Fe2)};
      Y1(1)=  {T.Na2O};
      Y1(2)=  {T.Cr2O3};
      Y1(3)=  {T.Al2O3};
    case 'clinopyroxene_xMg_Endmember_HPT'
        X1(1:10)=  {100.*T.apfu_Mg./(T.apfu_Mg+T.apfu_Fe2)};
        % if ismember(T.Properties.VariableNames,))
      Y1(1)=  {T.StrctFrm_Xjd};
      Y1(2)=  {T.StrctFrm_Xaeg};
      Y1(3)=  {T.StrctFrm_Xdihd};
      Y1(4)=  {T.StrctFrm_XCats};
      Y1(5)=  {T.StrctFrm_Xkos};
      Y1(6)=  {T.StrctFrm_XKkos};
      Y1(7)=  {T.StrctFrm_XKjd};
      Y1(8)=  {T.StrctFrm_XTicpx};
      Y1(9)=  {T.StrctFrm_XCaes};
      Y1(10)=  {T.StrctFrm_Xopx};
case 'clinopyroxene_xMg_Endmember'
        X1(1:7)=  {100.*T.apfu_Mg./(T.apfu_Mg+T.apfu_Fe2)};
        % if ismember(T.Properties.VariableNames,))
      Y1(1)=  {T.StrctFrm_Xwo};
      Y1(2)=  {T.StrctFrm_Xfs};
      Y1(3)=  {T.StrctFrm_Xen};
      Y1(4)=  {T.StrctFrm_Xjd};
      Y1(5)=  {T.StrctFrm_Xaeg};
      Y1(6)=  {T.StrctFrm_Xkos};
      Y1(7)=  {T.StrctFrm_Xquad};

end
try
    if plot_data==true
        if  exist('X1','var') &&  exist('Y1','var') && numel(X1)==numel(Y1)
            legend_test=false(numel(X1),1);
            for n=1:numel(X1)
                hold(ax2plot(n),'on');
                scatter(ax2plot(n),X1{n},Y1{n},symbsize,symb,'filled','MarkerFaceAlpha',mfa,'MarkerEdgeAlpha',mea,'MarkerEdgeColor',mec,'MarkerFaceColor',fil,'LineWidth',mlw)
                hold(ax2plot(n),'off');
                if sum(not(isnan(X1{n})))==numel(X1{n}) && sum(not(isnan(Y1{n})))==numel(Y1{n})
                    legend_test(n)=true;
                end
            end

            if  isfield(options,'legend') && isfield(options.legend,'String')
                legend_str=options.legend.String(any(legend_test));
            end
        end
    end
catch
    0
end

end

