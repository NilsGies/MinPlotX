function [ax2plot,options_definition,legend_str] = plot_tern_data_fun(T,ax2plot,options)

[~,options_definition]=plot_tern_fun([],[]);
if not(exist('T','var'))
    T=[];
end
if not(exist('ax2plot','var'))
    ax2plot=[];
end
if not(exist('options','var'))
    options=[];
end

legend_str={};
%% Plot Options definition
if  (isfield(options,'default') && options.default==true) || not(isfield(options,'plot_settings')) || isempty(options.plot_settings)  || not(isfield(options.plot_settings,'MarkerFaceColor')) || isempty(options.plot_settings.MarkerFaceColor)
    fil=[0 0.4470 0.7410];
    options_definition.plot_settings.MarkerFaceColor=fil;
else
    fil=options.plot_settings.MarkerFaceColor;
end

if  (isfield(options,'default') && options.default==true) || not(isfield(options,'plot_settings')) || isempty(options.plot_settings)  || not(isfield(options.plot_settings,'MarkerFaceAlpha')) || isempty(options.plot_settings.MarkerFaceAlpha)  ||  not(isfield(options.plot_settings.MarkerFaceAlpha,'Value')) || isempty(options.plot_settings.MarkerFaceAlpha.Value)
    mfa=3/8;
    options_definition.plot_settings.MarkerFaceAlpha=mfa;
else
    mfa=options.plot_settings.MarkerFaceAlpha.Value;
end

if  (isfield(options,'default') && options.default==true) || not(isfield(options,'plot_settings')) || isempty(options.plot_settings)  || not(isfield(options.plot_settings,'symbol')) || isempty(options.plot_settings.symbol)
    symb='o';
    options_definition.plot_settings.symbol=symb;
else
    symb=options.plot_settings.symbol;
end

if  (isfield(options,'default') && options.default==true) || not(isfield(options,'plot_settings')) || isempty(options.plot_settings)  || not(isfield(options.plot_settings,'MarkerEdgeColor')) || isempty(options.plot_settings.MarkerEdgeColor)
    if strcmp(symb,'.')
        mec=fil;
    else
        mec=[0 0 0];
    end
    options_definition.plot_settings.MarkerEdgeColor=mec;
else
    mec=options.plot_settings.MarkerEdgeColor;

end

if  (isfield(options,'default') && options.default==true) || not(isfield(options,'plot_settings')) || isempty(options.plot_settings)  || not(isfield(options.plot_settings,'MarkerEdgeAlpha')) || isempty(options.plot_settings.MarkerEdgeAlpha)  || not(isfield(options.plot_settings.MarkerEdgeAlpha,'Value')) || isempty(options.plot_settings.MarkerEdgeAlpha.Value) 
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

if  (isfield(options,'default') && options.default==true) || not(isfield(options,'plot_settings')) || isempty(options.plot_settings)  || not(isfield(options.plot_settings,'symbol_size')) || isempty(options.plot_settings.symbol_size)
    symbsize=100;
    options_definition.plot_settings.symbol_size=symbsize;
else
    symbsize= options.plot_settings.symbol_size;
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
    [ax2plot]=plot_tern_test([],options);
end

switch options.type.Value
    case 'tern_woenfe'
        if ismember('HighPT_Endmembers',T.Properties.VariableNames) && any(T.HighPT_Endmembers)==1

            XCa=T.apfu_Ca./(T.apfu_Ca+T.apfu_Mg+T.apfu_Fe2+T.apfu_Mn2);
            XFe=T.apfu_Fe2./(T.apfu_Ca+T.apfu_Mg+T.apfu_Fe2+T.apfu_Mn2);

            Xfs=(XFe.*(T.StrctFrm_Xdihd+T.StrctFrm_Xopx))./(T.StrctFrm_Xdihd+T.StrctFrm_Xopx);
            Xwo=(XCa.*(T.StrctFrm_Xdihd+T.StrctFrm_Xopx))./(T.StrctFrm_Xdihd+T.StrctFrm_Xopx);
            X1=0.5.*(Xwo)+(Xfs);
            Y1=(Xwo)*(cos(30*pi()/180));

        elseif ismember('StrctFrm_Xfs',T.Properties.VariableNames) &&ismember('StrctFrm_Xwo',T.Properties.VariableNames)
            %normalize Wo + En + Fs to 1
            Xfs=T.StrctFrm_Xfs./T.StrctFrm_Xquad;
            Xwo=T.StrctFrm_Xwo./T.StrctFrm_Xquad;
            X1=0.5.*(Xwo)+(Xfs);
            Y1=(Xwo)*(cos(30*pi()/180));
        end


    case 'tern_qjdaeg'
         if ismember('HighPT_Endmembers',T.Properties.VariableNames) && any(T.HighPT_Endmembers)==1
            %normalize quad + jd + aeg to 1
            quad=T.StrctFrm_Xdihd+T.StrctFrm_Xopx;

            Xquad=quad./(quad+T.StrctFrm_Xjd+T.StrctFrm_Xaeg);
            Xaeg=T.StrctFrm_Xaeg./(quad+T.StrctFrm_Xjd+T.StrctFrm_Xaeg);

            X1=0.5.*(Xquad)+(Xaeg);
            Y1=(Xquad)*(cos(30*pi()/180));
        elseif ismember('StrctFrm_Xquad',T.Properties.VariableNames) &&ismember('StrctFrm_Xaeg',T.Properties.VariableNames)
            
            %normalize quad + jd + aeg to 1
            Xquad=T.StrctFrm_Xquad./(T.StrctFrm_Xquad+T.StrctFrm_Xjd+T.StrctFrm_Xaeg);
            Xaeg=T.StrctFrm_Xaeg./(T.StrctFrm_Xquad+T.StrctFrm_Xjd+T.StrctFrm_Xaeg);

            X1=0.5.*(Xquad)+(Xaeg);
            Y1=(Xquad)*(cos(30*pi()/180));
        end
    case {'tern_feldspar_simple','tern_feldspar_complex','tern_feldspar_complex2'}
        X1=0.5.*(T.StrctFrm_Xan)+(T.StrctFrm_Xor);
        Y1=(T.StrctFrm_Xan)*(cos(30*pi()/180));
    case 'tern_ol'

        X1=0.5.*(T.StrctFrm_XLrn)+(T.StrctFrm_XFa+T.StrctFrm_XTep);
        Y1=(T.StrctFrm_XLrn)*(cos(30*pi()/180));
    case 'tern_grt'
      
        %Normalization of Xgrs + Xprp + Xalm + Xsps to 1
        Xgrs=T.StrctFrm_Xgrs./(T.StrctFrm_Xgrs+T.StrctFrm_Xprp+T.StrctFrm_Xalm+T.StrctFrm_Xsps);
        Xprp=T.StrctFrm_Xprp./(T.StrctFrm_Xgrs+T.StrctFrm_Xprp+T.StrctFrm_Xalm+T.StrctFrm_Xsps);

        %transforms the data to ternary space
        X1=0.5.*(Xgrs)+(Xprp);
        Y1=(Xgrs)*(cos(30*pi()/180));
 case 'tern_grt_half'
           
        %Normalization of Xgrs + Xprp + Xalm + Xsps to 1
        Xgrs=T.StrctFrm_Xgrs./(T.StrctFrm_Xgrs+T.StrctFrm_Xprp+T.StrctFrm_Xalm+T.StrctFrm_Xsps);
        Xprp=T.StrctFrm_Xprp./(T.StrctFrm_Xgrs+T.StrctFrm_Xprp+T.StrctFrm_Xalm+T.StrctFrm_Xsps);

        %transforms the data to ternary space
        X1=0.5.*(Xgrs)+(Xprp);
        Y1=(Xgrs)*(cos(30*pi()/180));

        X1(Y1>0.34)=nan;
        Y1(Y1>0.34)=nan;

    case 'tern_grtTi'
        %transforms the data to ternary space
        if ismember('StrctFrm_Xmmt',T.Properties.VariableNames) &&ismember('StrctFrm_Xslo',T.Properties.VariableNames)


            %Normalization of Xslo + Xmmt + Xadr to 1
            Xslo=T.StrctFrm_Xslo./(T.StrctFrm_Xslo+T.StrctFrm_Xmmt+T.StrctFrm_Xadr);
            Xmmt=T.StrctFrm_Xmmt./(T.StrctFrm_Xslo+T.StrctFrm_Xmmt+T.StrctFrm_Xadr);

            X1=0.5.*(Xmmt)+(Xslo);
            Y1=(Xmmt)*(cos(30*pi()/180));
        end
    case 'tern_grtoct'
        if ismember('StrctFrm_Fe3_Oct',T.Properties.VariableNames) && ismember('StrctFrm_Al_Oct',T.Properties.VariableNames) && ismember('StrctFrm_Cr_Oct',T.Properties.VariableNames)
            Fe3=T.StrctFrm_Fe3_Oct./(T.StrctFrm_Fe3_Oct+T.StrctFrm_Ti_Oct+T.StrctFrm_Al_Oct+T.StrctFrm_Cr_Oct);
            Cr=(T.StrctFrm_Cr_Oct+T.StrctFrm_Ti_Oct)./(T.StrctFrm_Fe3_Oct+T.StrctFrm_Ti_Oct+T.StrctFrm_Al_Oct+T.StrctFrm_Cr_Oct);
            %Al=T.StrctFrm_Al_Oct./(T.StrctFrm_Fe3_Oct+T.StrctFrm_Ti_Oct+T.StrctFrm_Al_Oct+T.StrctFrm_Cr_Oct);
            %transforms the data to ternary space
            X1=0.5.*(Fe3)+(Cr);
            Y1=(Fe3)*(cos(30*pi()/180));
        end
    case 'tern_mscelprl'

        %only plots Dioctahedral, 50 % rule
        condition=T.StrctFrm_XDiOct >= 0.50;

        %normalize the ternary components to 1
        XAlCelN=T.StrctFrm_XAlcel./(T.StrctFrm_XAlcel+T.StrctFrm_Xprl+T.StrctFrm_Xms);
        XPrlN=T.StrctFrm_Xprl./(T.StrctFrm_XAlcel+T.StrctFrm_Xprl+T.StrctFrm_Xms);

        %transforms the data to ternary space
        X1=0.5.*(XAlCelN(condition))+(XPrlN(condition));
        Y1=(XAlCelN(condition))*(cos(30*pi()/180));

    case 'tern_mica_FOHCl'
             if ismember('StrctFrm_F_A',T.Properties.VariableNames) &&ismember('StrctFrm_OH_A',T.Properties.VariableNames)
   
        F=T.StrctFrm_F_A./(T.StrctFrm_F_A+T.StrctFrm_Cl_A+T.StrctFrm_OH_A);
        OH=T.StrctFrm_OH_A./(T.StrctFrm_F_A+T.StrctFrm_Cl_A+T.StrctFrm_OH_A);
        %transforms the data to ternary space
        X1=0.5.*(F)+(OH);
        Y1=(F)*(cos(30*pi()/180));
             end
    case 'tern_ctd'
        XFe=T.apfu_Fe2./(T.apfu_Fe2+T.apfu_Mn2+T.apfu_Mg);
        XMn=T.apfu_Mn2./(T.apfu_Fe2+T.apfu_Mn2+T.apfu_Mg);
        %transforms the data to ternary space
        X1=0.5.*(XFe)+(XMn);
        Y1=(XFe)*(cos(30*pi()/180));
    case 'tern_sp1'
        Fe3=(T.StrctFrm_Fe3_B)./(T.StrctFrm_Fe3_B+T.StrctFrm_Al_B+T.StrctFrm_Cr_B); %fraction of Fe3+ and Ti
        Al=(T.StrctFrm_Al_B)./(T.StrctFrm_Fe3_B+T.StrctFrm_Al_B+T.StrctFrm_Cr_B); %fraction of Al
        %transforms the data to ternary space
        X1=0.5.*(Fe3)+(Al);
        Y1=(Fe3)*(cos(30*pi()/180));
    case 'tern_sp2'
        Fe2=T.StrctFrm_Fe2_A./(T.StrctFrm_Ni_A+T.StrctFrm_Zn_A+T.StrctFrm_Fe2_A+T.StrctFrm_Mn2_A+T.StrctFrm_Mg_A);
        Mg=T.StrctFrm_Mg_A./(T.StrctFrm_Ni_A+T.StrctFrm_Zn_A+T.StrctFrm_Fe2_A+T.StrctFrm_Mn2_A+T.StrctFrm_Mg_A);
        %transforms the data to ternary space
        X1=0.5.*(Fe2)+(Mg);
        Y1=(Fe2)*(cos(30*pi()/180));
    case 'tern_ap'
        %variables for the ternary plot
        OH = (T.apfu_OH + T.StrctFrm_C4_Z + T.(char('apfu_S1-')) + T.(char('apfu_S2-')))./T.(char('StrctFrm_Anion_Sum')); %OH + S
        F = T.apfu_F./T.(char('StrctFrm_Anion_Sum')); %F
        X1=0.5.*(F)+(OH);
        Y1=(F)*(cos(30*pi()/180));
    case 'tern_scap'
        Cl=T.apfu_Cl./(T.apfu_Cl+T.apfu_S+T.apfu_C);
        Sulfate=T.apfu_S./(T.apfu_Cl+T.apfu_S+T.apfu_C);
        %transforms the data to ternary space
        X1=0.5.*(Cl)+(Sulfate);
        Y1=(Cl)*(cos(30*pi()/180));
    case 'tern_lws'

        %Normalization of Cr + Ti*4 + Fe3 to 1
        XCr=T.apfu_Cr./((T.apfu_Ti.*4)+T.apfu_Cr+T.apfu_Fe3);
        XFe=T.apfu_Fe3./((T.apfu_Ti.*4)+T.apfu_Cr+T.apfu_Fe3);

        %transforms the data to ternary space
        X1=0.5.*(XCr)+(XFe);
        Y1=(XCr)*(cos(30*pi()/180));

    case 'tern_carbonate_Ca_Mg_Fe'

        %Normalization of Ca + Mg+ Fe
        XCa=T.apfu_Ca./((T.apfu_Mg)+T.apfu_Ca+T.apfu_Fe);
        XMg=T.apfu_Fe./((T.apfu_Mg)+T.apfu_Ca+T.apfu_Fe);

        %transforms the data to ternary space
        X1=0.5.*(XCa)+(XMg);
        Y1=(XCa)*(cos(30*pi()/180));

    case 'tern_carbonate_Ca_Mg_Mn'

        %Normalization of Ca + Mg+ Mn
        XCa=T.apfu_Ca./((T.apfu_Mg)+T.apfu_Ca+T.apfu_Mn);
        XMg=T.apfu_Fe./((T.apfu_Mg)+T.apfu_Ca+T.apfu_Mn);

        %transforms the data to ternary space
        X1=0.5.*(XCa)+(XMg);
        Y1=(XCa)*(cos(30*pi()/180));

    case 'custom_tern'
      
        if isfield(options,'X1')&& not(isempty(options.X1))
            X0=T.(char(options.X1));
        else
            X0=(1:size(T,1))';
        end
        if isfield(options,'Y1')&& not(isempty(options.Y1))
            Y0=T.(char(options.Y1));
        else
            Y0=(1:size(T,1))';
        end
        if isfield(options,'Z1')&& not(isempty(options.Z1))
            Z0=T.(char(options.Z1));
        else
            Z0=(1:size(T,1))';
        end

     


        X1=0.5.*(Y0./(X0+Y0+Z0))+(Z0./(X0+Y0+Z0));
        Y1=(Y0./(X0+Y0+Z0))*(cos(30*pi()/180));

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
    scatter(ax2plot,X1(:),Y1(:),symbsize,symb,'filled','MarkerFaceAlpha',mfa,'MarkerEdgeAlpha',mea,'MarkerEdgeColor',mec,'MarkerFaceColor',fil,'LineWidth',mlw)

end


if exist('X1','var') && isfield(options,'legend') && isfield(options.legend,'String')
    legend_str=options.legend.String(any(not(isnan(X1)) & not(isnan(Y1))));
end
end
end

