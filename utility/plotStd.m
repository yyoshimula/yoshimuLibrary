%[text] # plot with standard deviation
%[text] standard deviationも一緒にplot
%[text] ## inputs
%[text] `t`: value along x-axis, Nx1 vector
%[text] `meanVal: mean value, Nx1 vector`
%[text] `stdVal:` standard deviation, Nx1 vector
%[text] ## outputs
%[text] NA (figure is plotted)
%[text] ## note
%[text] NA
%[text] ## references 
%[text] NA
%[text] ## revisions
%[text] 20230822  y.yoshimura, y.yoshimula@gmail.com
%[text] See also fig4Paper.
function plotStd(t_, meanVal, stdVal, lineColor, regionColor)
arguments
    t_ (:,1) {mustBeNumeric}
    meanVal (:,1) {mustBeNumeric}
    stdVal (:,1) {mustBeNumeric}
    lineColor {mustBeText} = 'r' 
    regionColor {mustBeText} = 'k' 
end

plot(t_, meanVal, lineColor);
% semilogx(x, meanVal, lineColor);
hold on
ar = area(t_, [meanVal - stdVal, stdVal + stdVal]);
set(ar(1), 'FaceColor', 'none', 'FaceAlpha', 0.1);

set(ar(2), 'FaceColor', regionColor, 'FaceAlpha', 0.1);

end

%[appendix]{"version":"1.0"}
%---
