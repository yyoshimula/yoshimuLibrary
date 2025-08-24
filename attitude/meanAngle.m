function thetaAve = meanAngle(thetaArray, w)
    thetaArray = thetaArray(:);
    w = w(:);
    w = w ./ sum(w);

    x = sum(w .* cos(thetaArray));
    y = sum(w .* sin(thetaArray));
    
    thetaAve = atan2(y, x);
end