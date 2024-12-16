



function cmap_struct=ReadColorMaps(file)
if not(exist('file','var'))
    file='XMap_ColorMaps.txt';
end

fid = fopen(file,'r');
            Compt = 0;
            
            cmap_struct(1).Name = 'None';
            cmap_struct(1).Code = 0;
            
            Compt = 0;
            ErrorLoad = 0;
            while 1
                tline = fgetl(fid);
                
                if isequal(tline,-1)
                    break
                end
                
                if length(tline) >= 1
                    if isequal(tline(1),'>')
                        Compt = Compt+1;
                        cmap_struct(Compt).Name = tline(3:end);
                        Row = 0;
                        while 1
                            tline = fgetl(fid);
                            if isequal(tline,-1) || isequal(tline,'')
                                break
                            end
                            NUM = strread(tline,'%f');
                            Row = Row+1;
                            cmap_struct(Compt).Code(Row,1:3) = NUM(1:3);
                        end
                        
                    end
                end
            end
            fclose(fid);

        
                cmap_struct(end+1).Name= 'parula (MATLAB)';
                    cmap_struct(end).Code(:,1:3)=parula(256);
                cmap_struct(end+1).Name= 'turbo (MATLAB)';
                    cmap_struct(end).Code(:,1:3)=turbo(256);
                cmap_struct(end+1).Name= 'hsv (MATLAB)';
                    cmap_struct(end).Code(:,1:3)=hsv(256);
                cmap_struct(end+1).Name= 'hot (MATLAB)';
                    cmap_struct(end).Code(:,1:3)=hot(256);
                cmap_struct(end+1).Name= 'cool (MATLAB)';
                    cmap_struct(end).Code(:,1:3)=cool(256);
                cmap_struct(end+1).Name= 'spring (MATLAB)';
                    cmap_struct(end).Code(:,1:3)=spring(256);
                cmap_struct(end+1).Name= 'summer (MATLAB)';
                    cmap_struct(end).Code(:,1:3)=summer(256);
                cmap_struct(end+1).Name= 'autumn (MATLAB)';
                    cmap_struct(end).Code(:,1:3)=autumn(256);
                cmap_struct(end+1).Name= 'winter (MATLAB)';
                    cmap_struct(end).Code(:,1:3)=winter(256);
                cmap_struct(end+1).Name= 'gray (MATLAB)';
                    cmap_struct(end).Code(:,1:3)=gray(256);
                cmap_struct(end+1).Name= 'bone (MATLAB)';
                    cmap_struct(end).Code(:,1:3)=bone(256);
                cmap_struct(end+1).Name= 'copper (MATLAB)';
                    cmap_struct(end).Code(:,1:3)=copper(256);
                cmap_struct(end+1).Name= 'pink (MATLAB)';
                    cmap_struct(end).Code(:,1:3)=pink(256);
                cmap_struct(end+1).Name= 'jet (MATLAB)';
                    cmap_struct(end).Code(:,1:3)=jet(256);
                cmap_struct(end+1).Name= 'lines (MATLAB)';
                    cmap_struct(end).Code(:,1:3)=lines(256);
                cmap_struct(end+1).Name= 'colorcube (MATLAB)';
                    cmap_struct(end).Code(:,1:3)=colorcube(256);
                cmap_struct(end+1).Name= 'prism (MATLAB)';
                    cmap_struct(end).Code(:,1:3)=prism(256);
                cmap_struct(end+1).Name= 'flag (MATLAB)';
                    cmap_struct(end).Code(:,1:3)=flag(256);
                    cmap_struct(end+1).Name= 'white (MATLAB)';
                    cmap_struct(end).Code(:,1:3)=white(256);
                    cmap_struct(end+1).Name= 'ABC';
                    cmap_struct(end).Code(:,1:3)=123;
                    cmap_struct(end+1).Name= 'custom';
                    cmap_struct(end).Code(:,1:3)=white(256).*0;

end

%     T1Options_ColormapDropDown.Items =  Items; %extractfield(ColorMaps,'Name');


