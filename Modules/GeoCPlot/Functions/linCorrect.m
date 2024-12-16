function [y_cor, correction_line]=linCorrect(y)
y_cor=zeros(size(y));
for q=1:size(y,2)
    y_corline=y(not(isnan(y(:,q))),q);
    
        if isempty(y_corline)
            y(:,q)=nan;
            continue

        end

        correction_line=linspace(min(y_corline(1,1)),max(y_corline(end,1)),size(y,1));
   
    y(:,q)=y(:,q)-correction_line';
    y(( y(:,q)<0),q)=0;
    y_cor(:,q)=y(:,q);
end
end
% function y_cor=linCorrect(y)
% y_cor=zeros(size(y));
% for q=1:size(y,2)
%     y_corline=y(not(isnan(y(:,q))));
%     correction_line=linspace(min(y_corline(1,1)),max(y_corline(end,1)),size(y,1));
%     y(:,q)=y(:,q)-correction_line';
%     y(( y(:,q)<0),q)=0;
%     y_cor(:,q)=y(:,q);
% end
% end



