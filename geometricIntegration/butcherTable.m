%[text] # butcher table for geometric integration
%[text] ## inputs
%[text] `n`: the order of geometric integration method (NOT the number of stage)
%[text] method: RK = Runge–Kutta method
%[text]               CG = Crouch–Grossman method
%[text] ## outputs
%[text] `a, b, c: coefficients`
%[text] ## note
%[text] NA
%[text] ## references 
%[text] to be added
%[text] ## revisions
%[text] 20240109  y.yoshimura y.yoshimula@gmail.com
function [a, b, c] = butcherTable(n, method)


switch method
%[text] ## Crouch–Grossman method
    case 'CG'
%[text] ### three-stage third-order Crouch–Grossman method (CG3)
        if n == 3
            a = [0 0 0
                3/4 0 0
                119/216 17/108 0];
            b = [13/51
                -2/3
                24/17];
            c = [0
                3/4
                17/24];
%[text] ### five-stage fourth-order Crouch–Grossman method (CG4)
        elseif n == 4
            a = [0 0 0 0 0
                0.8177227988124852 0 0 0 0
                0.3199876375476427 0.0659864263556022 0 0 0
                0.9214417194464946 0.4997857776773573 -1.0969984448371582 0 0
                0.3552358559023322 0.2390958372307326 1.3918565724203246 -1.1092979392113565 0];

            b = [0.1370831520630755
                -0.0183698531564020
                0.7397813985370780
                -0.1907142565505889
                0.3322195591068374];

            c = [0.0
                0.8177227988124852
                0.3859740639032449
                0.3242290522866937
                0.8768903263420429];

        else
            error('no data for this order')

        end
%[text] ##  Runge–Kutta method 
    case 'RK'
%[text] ### three-stage third-order Runge–Kutta method (RK3)
        if n == 3
            a = [0 0 0
                1/2 0 0
                -1 2 0];
            b = [1/6
                2/3
                1/6];
            c = [0
                1/2
                1];
%[text] ### four-stage fourth-order Runge–Kutta method (RK4)
        elseif n == 4
            a = [0 0 0 0
                1/2 0 0 0
                0 1/2 0 0
                0 0 1 0];
            b = [1/6
                1/3
                1/3
                1/6];
            c = [0
                1/2
                1/2
                1];

        else
            error('no data for this order')

        end


    otherwise
        error('no data for this method/order')

end

end

%[appendix]{"version":"1.0"}
%---
