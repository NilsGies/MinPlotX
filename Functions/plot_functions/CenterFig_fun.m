% Custom function to center the figure
function CenterFig_fun(hfig)
try
    screenSize = get(0, 'ScreenSize');
    figPos = hfig.Position;
    figWidth = figPos(3);
    figHeight = figPos(4);
    newX = (screenSize(3) - figWidth) / 2;
    newY = (screenSize(4) - figHeight) / 2;
    hfig.Position=[newX, newY, figWidth, figHeight];
    figure(hfig)
catch ME
    disp('CenterFig_fun failed!')
    disp(ME.message)
end
end
