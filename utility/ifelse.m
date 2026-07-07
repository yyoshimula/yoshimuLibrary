%[text] # inline conditional (ternary-like) helper
%[text] condition が真なら trueValue を，偽なら falseValue を返す
%[text] ## inputs
%[text] `condition`: logical scalar
%[text] `trueValue`: value returned when condition is true
%[text] `falseValue`: value returned when condition is false
%[text] ## outputs
%[text] `result`: trueValue if condition is true, otherwise falseValue
%[text] ## note
%[text] NA
function result = ifelse(condition, trueValue, falseValue)
    if condition
        result = trueValue;
    else
        result = falseValue;
    end
end

%[appendix]{"version":"1.0"}
%---
