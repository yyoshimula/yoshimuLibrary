%[text] # sample program for reading .obj file
%[text] `.objファイルを読み込んで表示`
%[text] ## note
%[text] NA
%[text] ## references
%[text] NA
%[text] ## revisions
%[text] 20220612  y.yoshimura, y.yoshimula@gmail.com
%[text] See also readSC, showSC.
close all
clear
clc

%% read object
sat = readSC('boxWing4.obj'); %[output:4c0d6f1e]

% check satellite shape
showSC(sat,'Normal','on')
axis equal

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[output:4c0d6f1e]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"関数または変数 'nPolygon' が認識されません。\n\nエラー: <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('calcAreaObj', '\/Users\/yyoshimula\/Dropbox\/MATLAB\/yoshimuLibrary\/object\/calcAreaObj.m', 42)\" style=\"font-weight:bold\">calcAreaObj<\/a> (<a href=\"matlab: opentoline('\/Users\/yyoshimula\/Dropbox\/MATLAB\/yoshimuLibrary\/object\/calcAreaObj.m',42,0)\">行 42<\/a>)\n    if nPolygon == 3\n    ^^^^^^^^^^^^^^^^\nエラー: <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('readSC', '\/Users\/yyoshimula\/Dropbox\/MATLAB\/yoshimuLibrary\/object\/readSC.m', 25)\" style=\"font-weight:bold\">readSC<\/a> (<a href=\"matlab: opentoline('\/Users\/yyoshimula\/Dropbox\/MATLAB\/yoshimuLibrary\/object\/readSC.m',25,0)\">行 25<\/a>)\n[sat.area, sat.pos] = calcAreaObj(sat);  % face areas and position\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"}}
%---
