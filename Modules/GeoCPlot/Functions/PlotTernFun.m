function PlotTernFun(DataIn,maxSP,ax2plot)

for k=1:numel(DataIn.DataMat)


    VarNames=DataIn.VarNames;
    classes=DataIn.classNumber{k};
    ClassNames=DataIn.ClassNames{k};
    cmap=DataIn.cmap{k};
    Symbol=DataIn.Symbol{k};
    SymbolSize=DataIn.SymbolSize{k};
    DataMatrix = DataIn.DataMat{k};
    Active = DataIn.Active{k};



    scale=false;

   


    nclasses_=unique(classes(not(isnan(classes))));


    nDigit = size(DataIn.DataMat{k},2);
    nOnes  = 3;
    index = nchoosek(1:nDigit, nOnes);  % Indices of 2 ones in 4 elements
    counter=100;

    sc=get(0, 'ScreenSize');
    sc(3)=sc(4);
    for n=1:size(index,1)
        disp(n)
        A_idx=index(n,1);
        B_idx=index(n,2);
        C_idx=index(n,3);


        counter=counter+1;
        if counter>maxSP
            counter=1;
            figure
            set(gcf,'position',sc)
        end
%         if exist("ax2plot","var")
%       ax2plt=ax2plot(n);
%         else
        nexttile
        %end



        for j=1:numel(nclasses_)
            if Active(j)==true
                nclasses=nclasses_(j);



            end

            A=DataMatrix(classes==nclasses,A_idx);
            B=DataMatrix(classes==nclasses,B_idx);
            C=DataMatrix(classes==nclasses,C_idx);

            if scale==true
                DataMatrix_2=[normalize(DataMatrix(:,A_idx),'range',[min(DataMatrix(:,[A_idx B_idx]),[],'all') max(DataMatrix(:,[A_idx B_idx]),[],'all')]) normalize(DataMatrix(:,B_idx),'range',[min(DataMatrix(:,[A_idx B_idx]),[],'all') max(DataMatrix(:,[A_idx B_idx]),[],'all')]) DataMatrix(:,C_idx)];
                %  DataMatrix_2=[normalize(DataMatrix(:,A_idx),'range',[min(DataMatrix(:,A_idx)) max(DataMatrix(:,A_idx))]) normalize(DataMatrix(:,B_idx),'range',[min(DataMatrix(:,B_idx)) max(DataMatrix(:,B_idx))]) normalize(DataMatrix(:,C_idx),'range',[min(DataMatrix(:,C_idx)) max(DataMatrix(:,C_idx))])];
                DataMatrix_2=normalize([DataMatrix_2(:,1) DataMatrix_2(:,2) DataMatrix_2(:,3)]','range')';%'.*(1/3);
                % DataMatrix_2=(DataMatrix_3.*DataMatrix_2);
                A=DataMatrix_2(classes==nclasses,1);
                B=DataMatrix_2(classes==nclasses,2);
                C=DataMatrix_2(classes==nclasses,3);
            end

            ternplot(A, B, C,Symbol{j},'MarkerFaceColor',cmap(j,:),'MarkerSize',SymbolSize(j),'MarkerEdgeColor',cmap(j,:));
            ternlabel(replace(VarNames(A_idx),'_',' '),replace(VarNames(B_idx),'_',' '),replace(VarNames(C_idx),'_',' '));
            hold on
        end
    end
end
end