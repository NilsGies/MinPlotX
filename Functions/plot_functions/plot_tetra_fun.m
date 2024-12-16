function [ax2plot,options_definition]=plot_tetra_fun (ax2plot,options)
if not(exist('ax2plot','var'))
    ax2plot=[];
end
if not(exist('options','var'))
    options=[];
end
%% options definition
options_definition.tetraplot.Value=true;

%change all to linecolor_1 and define as the following to make it failsafe
if  (isfield(options,'default') && options.default==true) || not(isfield(options,'plot_settings')) || isempty(options.plot_settings)  || not(isfield(options.plot_settings,'linecolor_1')) || isempty(options.plot_settings.linecolor_1)
    linecolor_1=[0.3 0.3 0.3];
    options_definition.plot_settings.linecolor_1.Value=linecolor_1;
    options_definition.linecolor_1.description='Linecolor for odd plotting intervals';
else
    linecolor_1=options.plot_settings.linecolor_1.Value;
end

if  (isfield(options,'default') && options.default==true) || not(isfield(options,'plot_settings')) || isempty(options.plot_settings)  || not(isfield(options.plot_settings,'linecolor_2')) || isempty(options.plot_settings.linecolor_2)
    linecolor_2=[0.3 0.3 0.3];
    options_definition.plot_settings.linecolor_2.Value=linecolor_2;
    options_definition.linecolor_2.description='Linecolor for even plotting intervals';
else
    linecolor_2=options.plot_settings.linecolor_2.Value;
end

if  (isfield(options,'default') && options.default==true) || not(isfield(options,'plot_settings')) || isempty(options.plot_settings)  || not(isfield(options.plot_settings,'LineStyle_1')) || isempty(options.plot_settings.LineStyle_1)
    LineStyle_1=':';
    options_definition.plot_settings.LineStyle_1.Value=LineStyle_1;
    options_definition.LineStyle_1.description='LineStyle for odd plotting intervals';
else
    LineStyle_1=options.plot_settings.LineStyle_1.Value;
end

if  (isfield(options,'default') && options.default==true) || not(isfield(options,'plot_settings')) || isempty(options.plot_settings)  || not(isfield(options.plot_settings,'LineStyle_2')) || isempty(options.plot_settings.LineStyle_2)
    LineStyle_2='--';
    options_definition.plot_settings.LineStyle_2.Value=LineStyle_2;
    options_definition.LineStyle_2.description='LineStyle for even plotting intervals';
else
    LineStyle_2=options.plot_settings.LineStyle_2.Value;
end


if  (isfield(options,'default') && options.default==true) || not(isfield(options,'plot_settings')) || isempty(options.plot_settings)  || not(isfield(options.plot_settings,'linewidth_1')) || isempty(options.plot_settings.linewidth_1)
    linewidth_1=0.5;
    options_definition.plot_settings.linewidth_1.Value=linewidth_1;
else
    linewidth_1=options.plot_settings.linewidth_1.Value;
end

if  (isfield(options,'default') && options.default==true) || not(isfield(options,'FontSize')) || isempty(options.FontSize)
    FontSize=12;
    options_definition.FontSize.Value=12;
else
    FontSize=options.FontSize.Value;
end

if  (isfield(options,'default') && options.default==true) || not(isfield(options,'AutoXlim')) || isempty(options.AutoXlim)
    AutoXlim=false;
    options_definition.AutoXlim.Value=AutoXlim;
else
    AutoXlim=options.AutoXlim.Value;
end

options_definition.type.Value='none';
options_definition.type.options={'none';
    'grt_Alm_Py_Grs_Sps';
    'custom_XYZ';
    };

options_definition.type.description={'none';
    'garnet: Almandine-Pyrope-Grossular-Spessartine Tetrahedron';
    'custom: Tetrahedron';
    };


if not(exist('options','var')) || isempty(options)
    if not(exist('ax2plot','var'))
        ax2plot=[];
    end
    return
end

if not(exist('ax2plot','var')) || isempty(ax2plot)
    figure
    ax2plot=nexttile;
end
hold on
       

ax2plot.FontSize=FontSize;

if (isfield(options,'default') && options.default==true)
    type=options.type.Value;
    options=options_definition;
    options.type.Value=type;
end
rotate_state=true;

if options.tetraplot.Value==true
    %% Mineral specific options
    if not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'grt_Alm_Py_Grs_Sps')
   
        str_grt1={'alm';'sps';'grs';'prp'}; %endmember text labels

    elseif not(isempty(options)) && isfield(options,'type') && strcmp(options.type.Value,'custom_XYZ')

    end

    if isfield(options.custom,'Interpreter')
            Interpreter=options.custom.Interpreter;
        else
            Interpreter='none';
        end

x1=[3 2 1 2 3 2 2 1 3];                             %tetrahedron
y1=[1 2.732 1 1.577350 1 1.577350 2.732 1 1 ];       
z1=[1 1 1 2.633 1 2.633  1 1 1 ];              
x2=1;                                      %annotation A
y2=0.85;
z2=1;                                       
x3=3;                                      %annotation B
y3=0.85;
z3=1;
x4=2.07;                                         %annotation C 
y4=2.7;
z4=1.1;
x5=2.03;                                         %annotation D
y5=1.66;
z5=2.72;  

%%
plot3(ax2plot,x1,y1,z1, 'k','handlevisibility','off')                 %plot tetrahedron
hold on 

axis equal 

 set(gcf, 'Color', [1,1,1])
 
hold on
%delete(gca)

% Define the data range
miv=min(0.1);
mav=max(100);
% Get the current colormap

% Plot the points
hold on

view(3)
 

if  not(isempty(options)) && isfield(options,'custom') && isfield(options.custom,'xlabel') && not(isempty(options.custom.xlabel))
    text(ax2plot,x2,y2,z2,options.custom.xlabel,'FontWeight','bold','Interpreter',Interpreter)
else
    text(ax2plot,x2,y2,z2,'A','FontWeight','bold','Interpreter',Interpreter)
end

if  not(isempty(options)) && isfield(options,'custom') && isfield(options.custom,'ylabel') && not(isempty(options.custom.ylabel))
    text(ax2plot,x3,y3,z3,options.custom.ylabel,'FontWeight','bold','Interpreter',Interpreter)
else
    text(ax2plot,x3,y3,z3,'B','FontWeight','bold','Interpreter',Interpreter)
end

if  not(isempty(options)) && isfield(options,'custom') && isfield(options.custom,'zlabel') && not(isempty(options.custom.zlabel))
    text(ax2plot,x4,y4,z4, options.custom.zlabel,'FontWeight','bold','Interpreter',Interpreter)
else
    text(ax2plot,x4,y4,z4, 'C','FontWeight','bold','Interpreter',Interpreter)
end

if  not(isempty(options)) && isfield(options,'custom') && isfield(options.custom,'dlabel') && not(isempty(options.custom.dlabel))
text(ax2plot,x5,y5,z5, options.custom.dlabel,'FontWeight','bold','Interpreter',Interpreter)
else
text(ax2plot,x5,y5,z5, 'D','FontWeight','bold','Interpreter',Interpreter)
end


 
ax2plot.XAxis.Visible=false;
ax2plot.YAxis.Visible=false;
ax2plot.ZAxis.Visible=false;
ax2plot.Color='none';

end
 
 
 
 


