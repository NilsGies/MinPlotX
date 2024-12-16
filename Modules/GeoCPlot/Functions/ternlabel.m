% TERNLABEL label ternary phase diagram
%   TERNLABEL('ALABEL', 'BLABEL', 'CLABEL') labels a ternary phase diagram created using TERNPLOT
%   
%   H = TERNLABEL('ALABEL', 'BLABEL', 'CLABEL') returns handles to the text objects created.
%   with the labels provided.  TeX escape codes are accepted.
%
%   See also TERNPLOT

% Author: Carl Sandrock 20020827

% To Do

% Modifications

% Modifiers

function h = ternlabel(A, B, C)
r(1) = text(0.0, -0.05, A, 'horizontalalignment', 'center');
r(2) = text(1, -0.05, B, 'rotation', 0, 'horizontalalignment', 'center');
r(3) = text(0.5, 0.9, 0.5, C, 'rotation', 0, 'horizontalalignment', 'center');

if nargout > 0
    h = r;
end;