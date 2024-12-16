function ElOxDataDef=ReadDefFiles
            
            fid = fopen('XMap_Def_Element.txt','r');
            Compt = 0;
            
            ElOxDataDef.ElList = [];
            ElOxDataDef.ElIdx = [];
            ElOxDataDef.ElMass = [];
            
            while 1
                tline = fgetl(fid);
                
                if isequal(tline,-1)
                    break
                end
                
                if length(tline) >= 1
                    if isequal(tline(1),'>')
                        
                        while 1
                            tline = fgetl(fid);
                            
                            if isequal(tline,-1) || isequal(tline,'')
                                break
                            end
                            
                            Compt = Compt + 1;
                            
                            TheStr = strread(tline,'%s');
                            
                            ElOxDataDef.ElList{Compt} = TheStr{1};
                            ElOxDataDef.ElIdx(Compt) = Compt;
                            ElOxDataDef.ElMass(Compt) = str2num(TheStr{2});
                        end
                    end
                end
            end
            
            fclose(fid);
            
            fid = fopen('XMap_Def_Oxide.txt','r');
            Compt = 0;
            
            ElOxDataDef.OxList = [];
            ElOxDataDef.OxElIdx = [];
            ElOxDataDef.OxNbCat = [];
            ElOxDataDef.OxNbOx = [];
            ElOxDataDef.OxMass = [];
            ElOxDataDef.OxFact = [];
            
            while 1
                tline = fgetl(fid);
                
                if isequal(tline,-1)
                    break
                end
                
                if length(tline) >= 1
                    if isequal(tline(1),'>')
                        
                        while 1
                            tline = fgetl(fid);
                            
                            if isequal(tline,-1) || isequal(tline,'')
                                break
                            end
                            
                            Compt = Compt + 1;
                            
                            TheStr = strread(tline,'%s');
                            
                            ElOxDataDef.OxList{Compt} = TheStr{1};
                            
                            Element = TheStr{2};
                            Idx = find(ismember(ElOxDataDef.ElList,Element));
                            
                            ElOxDataDef.OxElIdx(Compt) = Idx;
                            ElOxDataDef.OxNbCat(Compt) = str2num(TheStr{3});
                            ElOxDataDef.OxNbOx(Compt) = str2num(TheStr{4});
                            ElOxDataDef.OxMass(Compt) = str2num(TheStr{5});
                            ElOxDataDef.OxFact(Compt) = str2num(TheStr{6});
                            
                        end
                    end
                end
            end
            fclose(fid);
        end