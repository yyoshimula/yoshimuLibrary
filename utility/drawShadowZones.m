%[text] # draw shadow (eclipse) zones on the current axes
%[text] sunlitFlag に応じて，xregion で eclipse 中の区間を陰影表示する
%[text] ## inputs
%[text] `t_`: time vector
%[text] `sunlitFlag`: flag per time step, 1: sunlit, 0: shadow (eclipse)
%[text] ## outputs
%[text] NA (shaded regions are drawn on the current axes)
%[text] ## note
%[text] NA
function drawShadowZones(t_, sunlitFlag)
% t_: time vector
% sunlitFlag = 1: sunlit, 0: shadow
% sunlitFlagに応じて，xregionでeclipse中であることを示す

t_ = t_(:);
sunlitFlag = sunlitFlag(:);

% Check if sunlitFlag contains any 0
if all(sunlitFlag)
    return;
end

% Pad with 1 (sunlit) to handle edges where data starts/ends with 0
paddedFlag = [1; sunlitFlag; 1];
d = diff(paddedFlag);

% -1 indicates transition from 1 to 0 (start of eclipse)
% 1 indicates transition from 0 to 1 (end of eclipse)
startIdx = find(d == -1);
endIdx = find(d == 1) - 1;

% Apply xregion
xregion(t_(startIdx), t_(endIdx));

end