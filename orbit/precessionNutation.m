%[text] # precession and nutation based on IAU-2006/2000 theory
%[text] calculating DCM from CIRS to GCRF
%[text] jdTT:  julian centuries of terrestrial time (TT)
%[text] dX, dY: correction term (optional)
%[text] ## note
%[text] `iau06x.dat, iau06y.dat, iau06z.dat`
%[text] are required.
%[text] ## references 
%[text] Vallado, D. A., & McClain, W. D. (2001). Fundamentals of Astrodynamics and Applications. Springer Science & Business Media. 4th edition, p.213
%[text] ## revisions
%[text] 20230608  y.yoshimura, y.yoshimula@gmail.com
%[text] See also readIAU06, leapS.
function dcm = precessionNutation(jdTT, iau06, dX, dY)
% arguments
%     jdTT (:,1) {mustBeNumeric}
%     iau06
%     dX = 0.0
%     dY = 0.0
% end

tTT = jd2jdT(jdTT); % julian centuries of TT (terrestrial time)

%[text] ## delaunay variables
% the delaunay fundamental arguments, deg
l =  134.96340251  + ( 1717915923.2178 * tTT + ...
    31.8792 * tTT^2 + 0.051635 * tTT^3 - 0.00024470 * tTT^4 ) / 3600.0;
l1   =  357.52910918  + (  129596581.0481 *tTT - ...
    0.5532 * tTT^2 - 0.000136 * tTT^3 - 0.00001149 * tTT^4 )  / 3600.0;
f    =   93.27209062  + ( 1739527262.8478 *tTT - ...
    12.7512 *tTT^2 + 0.001037 * tTT^3 + 0.00000417 * tTT^4 )  / 3600.0;
d    =  297.85019547  + ( 1602961601.2090 *tTT - ...
    6.3706 *tTT^2 + 0.006593 * tTT^3 - 0.00003169 * tTT^4 )  / 3600.0;
omega=  125.04455501  + (   -6962890.5431 *tTT + ...
    7.4722 *tTT^2 + 0.007702 * tTT^3 - 0.00005939 * tTT^4 )  / 3600.0;

% the planetary arguments, deg
lonmer  = 252.250905494  + 149472.6746358  * tTT;
lonven  = 181.979800853  +  58517.8156748  * tTT;
lonear  = 100.466448494  +  35999.3728521  * tTT;
lonmar  = 355.433274605  +  19140.299314   * tTT;
lonjup  =  34.351483900  +   3034.90567464 * tTT;
lonsat  =  50.0774713998 +   1222.11379404 * tTT;
lonurn  = 314.055005137  +    428.466998313 * tTT;
lonnep  = 304.348665499  +    218.486200208 * tTT;
precrate= 1.39697137214 * tTT + 0.0003086 * tTT^2;

% deg to rad
l    = deg2rad(rem( l,360.0  )    ); % rad
l1   = deg2rad(rem( l1,360.0  )   );
f    = deg2rad(rem( f,360.0  )    );
d    = deg2rad(rem( d,360.0  )    );
omega= deg2rad(rem( omega,360.0  ));

lonmer= deg2rad(rem( lonmer,360.0 ));  % rad
lonven= deg2rad(rem( lonven,360.0 ));
lonear= deg2rad(rem( lonear,360.0 ));
lonmar= deg2rad(rem( lonmar,360.0 ));
lonjup= deg2rad(rem( lonjup,360.0 ));
lonsat= deg2rad(rem( lonsat,360.0 ));
lonurn= deg2rad(rem( lonurn,360.0 ));
lonnep= deg2rad(rem( lonnep,360.0 ));
precrate= deg2rad(rem( precrate,360.0 ));
%[text] ## X
%[text] fundameの書き方↓をvectorで計算するように変更している
%[text] ```matlabCodeExample
%[text] for i = 1306: -1 : 1
%[text]     tempval = iau06.a0xi(i,1) * l + a0xi(i,2) * l1 + a0xi(i,3) * f + a0xi(i,4) * d + a0xi(i,5) * omega + ...
%[text]         a0xi(i,6) * lonmer  + a0xi(i,7) * lonven  + a0xi(i,8) * lonear  + a0xi(i,9) * lonmar + ...
%[text]         a0xi(i,10) * lonjup + a0xi(i,11) * lonsat + a0xi(i,12) * lonurn + a0xi(i,13) * lonnep + a0xi(i,14) * precrate;
%[text]     xsum0 = xsum0 + axs0(i,1) * sin(tempval) + axs0(i,2) * cos(tempval);
%[text] end
%[text] ```
% ---------------- first find x
% the iers code puts the constants in here, however
% don't sum constants in here because they're larger than the last few terms

ind = 1:1306;
% tmp: (ind x 1) vector
tmp = iau06.a0xi(ind,1) .* l + iau06.a0xi(ind,2) .* l1 + iau06.a0xi(ind,3) .* f + iau06.a0xi(ind,4) .* d + iau06.a0xi(ind,5) .* omega + ...
    iau06.a0xi(ind,6) .* lonmer + iau06.a0xi(ind,7) .* lonven  + iau06.a0xi(ind,8) .* lonear  + iau06.a0xi(ind,9) .* lonmar + ...
    iau06.a0xi(ind,10) .* lonjup + iau06.a0xi(ind,11) .* lonsat + iau06.a0xi(ind,12) .* lonurn + iau06.a0xi(ind,13) .* lonnep + iau06.a0xi(ind,14) .* precrate;
xsum0 = sum(iau06.axs0(ind,1) .* sin(tmp) + iau06.axs0(ind,2) .* cos(tmp));

% note that the index changes here to j. this is because the iau06.a0xi etc
% indicies go from 1 to 1600, but there are 5 groups. the i index counts through each
% calculation, and j takes care of the individual summations. note that
% this same process is used for y and s.

ind = 1307:(1306+253);
% tmp: (ind x 1) vector
tmp = iau06.a0xi(ind,1) .* l + iau06.a0xi(ind,2) .* l1 + iau06.a0xi(ind,3) .* f + iau06.a0xi(ind,4) .* d + iau06.a0xi(ind,5) .* omega + ...
    iau06.a0xi(ind,6) .* lonmer + iau06.a0xi(ind,7) .* lonven  + iau06.a0xi(ind,8) .* lonear  + iau06.a0xi(ind,9) .* lonmar + ...
    iau06.a0xi(ind,10) .* lonjup + iau06.a0xi(ind,11) .* lonsat + iau06.a0xi(ind,12) .* lonurn + iau06.a0xi(ind,13) .* lonnep + iau06.a0xi(ind,14) .* precrate;

xsum1 = sum(iau06.axs0(ind,1) .* sin(tmp) + iau06.axs0(ind,2) .* cos(tmp));


ind = (1306 + 253 + 1):(1306 + 253 + 36);
% tmp: (ind x 1) vector
tmp = iau06.a0xi(ind,1) .* l + iau06.a0xi(ind,2) .* l1 + iau06.a0xi(ind,3) .* f + iau06.a0xi(ind,4) .* d + iau06.a0xi(ind,5) .* omega + ...
    iau06.a0xi(ind,6) .* lonmer + iau06.a0xi(ind,7) .* lonven  + iau06.a0xi(ind,8) .* lonear  + iau06.a0xi(ind,9) .* lonmar + ...
    iau06.a0xi(ind,10) .* lonjup + iau06.a0xi(ind,11) .* lonsat + iau06.a0xi(ind,12) .* lonurn + iau06.a0xi(ind,13) .* lonnep + iau06.a0xi(ind,14) .* precrate;

xsum2 = sum(iau06.axs0(ind,1) .* sin(tmp) + iau06.axs0(ind,2) .* cos(tmp));


ind = (1306 + 253 + 36 + 1):(1306 + 253 + 36 + 4);
% tmp: (ind x 1) vector
tmp = iau06.a0xi(ind,1) .* l + iau06.a0xi(ind,2) .* l1 + iau06.a0xi(ind,3) .* f + iau06.a0xi(ind,4) .* d + iau06.a0xi(ind,5) .* omega + ...
    iau06.a0xi(ind,6) .* lonmer + iau06.a0xi(ind,7) .* lonven  + iau06.a0xi(ind,8) .* lonear  + iau06.a0xi(ind,9) .* lonmar + ...
    iau06.a0xi(ind,10) .* lonjup + iau06.a0xi(ind,11) .* lonsat + iau06.a0xi(ind,12) .* lonurn + iau06.a0xi(ind,13) .* lonnep + iau06.a0xi(ind,14) .* precrate;

xsum3 = sum(iau06.axs0(ind,1) .* sin(tmp) + iau06.axs0(ind,2) .* cos(tmp));

ind = 1306 + 253 + 36 + 4 + 1;
% tmp: scalar
tmp = iau06.a0xi(ind,1) .* l + iau06.a0xi(ind,2) .* l1 + iau06.a0xi(ind,3) .* f + iau06.a0xi(ind,4) .* d + iau06.a0xi(ind,5) .* omega + ...
    iau06.a0xi(ind,6) .* lonmer + iau06.a0xi(ind,7) .* lonven  + iau06.a0xi(ind,8) .* lonear  + iau06.a0xi(ind,9) .* lonmar + ...
    iau06.a0xi(ind,10) .* lonjup + iau06.a0xi(ind,11) .* lonsat + iau06.a0xi(ind,12) .* lonurn + iau06.a0xi(ind,13) .* lonnep + iau06.a0xi(ind,14) .* precrate;

xsum4 = sum(iau06.axs0(ind,1) .* sin(tmp) + iau06.axs0(ind,2) .* cos(tmp));

X = -0.016617 + 2004.191898 * tTT - 0.4297829 * tTT^2 ...
    - 0.19861834 * tTT^3 - 0.000007578 * tTT^4 + 0.0000059285 * tTT^5; % arcsecond
X = arcs2rad(X) + xsum0 + xsum1 * tTT + xsum2 * tTT^2 + xsum3 * tTT^3 + xsum4 * tTT^4;  % rad


%[text] ## Y
ind = 1:962;
% tmp: (ind x 1) vector
tmp = iau06.a0yi(ind,1) .* l + iau06.a0yi(ind,2) .* l1 + iau06.a0yi(ind,3) .* f + iau06.a0yi(ind,4) .* d + iau06.a0yi(ind,5) .* omega + ...
    iau06.a0yi(ind,6) .* lonmer + iau06.a0yi(ind,7) .* lonven  + iau06.a0yi(ind,8) .* lonear  + iau06.a0yi(ind,9) .* lonmar + ...
    iau06.a0yi(ind,10) .* lonjup + iau06.a0yi(ind,11) .* lonsat + iau06.a0yi(ind,12) .* lonurn + iau06.a0yi(ind,13) .* lonnep + iau06.a0yi(ind,14) .* precrate;

ysum0 = sum(iau06.ays0(ind,1) .* sin(tmp) + iau06.ays0(ind,2) .* cos(tmp));


ind = 963:(962+277);
% tmp: (ind x 1) vector
tmp = iau06.a0yi(ind,1) .* l + iau06.a0yi(ind,2) .* l1 + iau06.a0yi(ind,3) .* f + iau06.a0yi(ind,4) .* d + iau06.a0yi(ind,5) .* omega + ...
    iau06.a0yi(ind,6) .* lonmer + iau06.a0yi(ind,7) .* lonven  + iau06.a0yi(ind,8) .* lonear  + iau06.a0yi(ind,9) .* lonmar + ...
    iau06.a0yi(ind,10) .* lonjup + iau06.a0yi(ind,11) .* lonsat + iau06.a0yi(ind,12) .* lonurn + iau06.a0yi(ind,13) .* lonnep + iau06.a0yi(ind,14) .* precrate;

ysum1 = sum(iau06.ays0(ind,1) .* sin(tmp) + iau06.ays0(ind,2) .* cos(tmp));

ind = (962 + 277 + 1):(962+277+30);
% tmp: (ind x 1) vector
tmp = iau06.a0yi(ind,1) .* l + iau06.a0yi(ind,2) .* l1 + iau06.a0yi(ind,3) .* f + iau06.a0yi(ind,4) .* d + iau06.a0yi(ind,5) .* omega + ...
    iau06.a0yi(ind,6) .* lonmer + iau06.a0yi(ind,7) .* lonven  + iau06.a0yi(ind,8) .* lonear  + iau06.a0yi(ind,9) .* lonmar + ...
    iau06.a0yi(ind,10) .* lonjup + iau06.a0yi(ind,11) .* lonsat + iau06.a0yi(ind,12) .* lonurn + iau06.a0yi(ind,13) .* lonnep + iau06.a0yi(ind,14) .* precrate;

ysum2 = sum(iau06.ays0(ind,1) .* sin(tmp) + iau06.ays0(ind,2) .* cos(tmp));

ind = (962+277+30+1):(962+277+30+5);
% tmp: (ind x 1) vector
tmp = iau06.a0yi(ind,1) .* l + iau06.a0yi(ind,2) .* l1 + iau06.a0yi(ind,3) .* f + iau06.a0yi(ind,4) .* d + iau06.a0yi(ind,5) .* omega + ...
    iau06.a0yi(ind,6) .* lonmer + iau06.a0yi(ind,7) .* lonven  + iau06.a0yi(ind,8) .* lonear  + iau06.a0yi(ind,9) .* lonmar + ...
    iau06.a0yi(ind,10) .* lonjup + iau06.a0yi(ind,11) .* lonsat + iau06.a0yi(ind,12) .* lonurn + iau06.a0yi(ind,13) .* lonnep + iau06.a0yi(ind,14) .* precrate;

ysum3 = sum(iau06.ays0(ind,1) .* sin(tmp) + iau06.ays0(ind,2) .* cos(tmp));

ind = 962+277+30+5+1;
% tmp: (ind x 1) vector
tmp = iau06.a0yi(ind,1) .* l + iau06.a0yi(ind,2) .* l1 + iau06.a0yi(ind,3) .* f + iau06.a0yi(ind,4) .* d + iau06.a0yi(ind,5) .* omega + ...
    iau06.a0yi(ind,6) .* lonmer + iau06.a0yi(ind,7) .* lonven  + iau06.a0yi(ind,8) .* lonear  + iau06.a0yi(ind,9) .* lonmar + ...
    iau06.a0yi(ind,10) .* lonjup + iau06.a0yi(ind,11) .* lonsat + iau06.a0yi(ind,12) .* lonurn + iau06.a0yi(ind,13) .* lonnep + iau06.a0yi(ind,14) .* precrate;

ysum4 = sum(iau06.ays0(ind,1) .* sin(tmp) + iau06.ays0(ind,2) .* cos(tmp));

Y = -0.006951 - 0.025896 * tTT - 22.4072747 * tTT^2 ...
    + 0.00190059 * tTT^3 + 0.001112526 * tTT^4 + 0.0000001358 * tTT^5;
Y = arcs2rad(Y) + ysum0 + ysum1 * tTT + ysum2 * tTT^2 + ysum3 * tTT^3 + ysum4 * tTT^4;  % rad
%[text] ## s
ind = 1:33;
% tmp: (ind x 1) vector
tmp = iau06.a0si(ind,1) .* l + iau06.a0si(ind,2) .* l1 + iau06.a0si(ind,3) .* f + iau06.a0si(ind,4) .* d + iau06.a0si(ind,5) .* omega + ...
    iau06.a0si(ind,6) .* lonmer + iau06.a0si(ind,7) .* lonven  + iau06.a0si(ind,8) .* lonear  + iau06.a0si(ind,9) .* lonmar + ...
    iau06.a0si(ind,10) .* lonjup + iau06.a0si(ind,11) .* lonsat + iau06.a0si(ind,12) .* lonurn + iau06.a0si(ind,13) .* lonnep + iau06.a0si(ind,14) .* precrate;

ssum0 = sum(iau06.ass0(ind,1) .* sin(tmp) + iau06.ass0(ind,2) .* cos(tmp));

ind = (33+1):(33+3);
% tmp: (ind x 1) vector
tmp = iau06.a0si(ind,1) .* l + iau06.a0si(ind,2) .* l1 + iau06.a0si(ind,3) .* f + iau06.a0si(ind,4) .* d + iau06.a0si(ind,5) .* omega + ...
    iau06.a0si(ind,6) .* lonmer + iau06.a0si(ind,7) .* lonven  + iau06.a0si(ind,8) .* lonear  + iau06.a0si(ind,9) .* lonmar + ...
    iau06.a0si(ind,10) .* lonjup + iau06.a0si(ind,11) .* lonsat + iau06.a0si(ind,12) .* lonurn + iau06.a0si(ind,13) .* lonnep + iau06.a0si(ind,14) .* precrate;

ssum1 = sum(iau06.ass0(ind,1) .* sin(tmp) + iau06.ass0(ind,2) .* cos(tmp));

ind = (33+3+1):(33+3+25);
% tmp: (ind x 1) vector
tmp = iau06.a0si(ind,1) .* l + iau06.a0si(ind,2) .* l1 + iau06.a0si(ind,3) .* f + iau06.a0si(ind,4) .* d + iau06.a0si(ind,5) .* omega + ...
    iau06.a0si(ind,6) .* lonmer + iau06.a0si(ind,7) .* lonven  + iau06.a0si(ind,8) .* lonear  + iau06.a0si(ind,9) .* lonmar + ...
    iau06.a0si(ind,10) .* lonjup + iau06.a0si(ind,11) .* lonsat + iau06.a0si(ind,12) .* lonurn + iau06.a0si(ind,13) .* lonnep + iau06.a0si(ind,14) .* precrate;

ssum2 = sum(iau06.ass0(ind,1) .* sin(tmp) + iau06.ass0(ind,2) .* cos(tmp));

ind = (33+3+25+1):(33+3+25+4);
% tmp: (ind x 1) vector
tmp = iau06.a0si(ind,1) .* l + iau06.a0si(ind,2) .* l1 + iau06.a0si(ind,3) .* f + iau06.a0si(ind,4) .* d + iau06.a0si(ind,5) .* omega + ...
    iau06.a0si(ind,6) .* lonmer + iau06.a0si(ind,7) .* lonven  + iau06.a0si(ind,8) .* lonear  + iau06.a0si(ind,9) .* lonmar + ...
    iau06.a0si(ind,10) .* lonjup + iau06.a0si(ind,11) .* lonsat + iau06.a0si(ind,12) .* lonurn + iau06.a0si(ind,13) .* lonnep + iau06.a0si(ind,14) .* precrate;

ssum3 = sum(iau06.ass0(ind,1) .* sin(tmp) + iau06.ass0(ind,2) .* cos(tmp));

ind = 33+3+25+4+1;
% tmp: (ind x 1) vector
tmp = iau06.a0si(ind,1) .* l + iau06.a0si(ind,2) .* l1 + iau06.a0si(ind,3) .* f + iau06.a0si(ind,4) .* d + iau06.a0si(ind,5) .* omega + ...
    iau06.a0si(ind,6) .* lonmer + iau06.a0si(ind,7) .* lonven  + iau06.a0si(ind,8) .* lonear  + iau06.a0si(ind,9) .* lonmar + ...
    iau06.a0si(ind,10) .* lonjup + iau06.a0si(ind,11) .* lonsat + iau06.a0si(ind,12) .* lonurn + iau06.a0si(ind,13) .* lonnep + iau06.a0si(ind,14) .* precrate;

ssum4 = sum(iau06.ass0(ind,1) .* sin(tmp) + iau06.ass0(ind,2) .* cos(tmp));

s = 0.000094 + 0.00380865 * tTT - 0.00012268 * tTT^2  ...
    - 0.07257411 * tTT^3 + 0.00002798 * tTT^4 + 0.00001562 * tTT^5;
s = -X * Y * 0.5 + arcs2rad(s) + ssum0 + ssum1 * tTT + ssum2 * tTT^2 + ssum3 * tTT^3 + ssum4 * tTT^4;  % rad


%[text] ### correction
% add corrections if available
X = X + dX; % rad
Y = Y + dY;

%[text] ### a
a = 0.5 + 0.125 * (X * X + Y * Y);

%[text] ## DCM

dcmTmp = [1 - a*X^2, -a*X*Y, X
    -a*X*Y, 1 - a*Y^2, Y
    -X, -Y, 1-a*(X^2 + Y^2)];

dcm = dcmTmp * dcm1axis(3, s);


end



%[appendix]{"version":"1.0"}
%---
