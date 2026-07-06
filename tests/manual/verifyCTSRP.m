%[text] # verifying Cook–Torrance SRP approximation
%[text] Cook–Torrance SRP approximationを数値的なSRP計算と比較して検証
%[text] 太陽方向を\[0, 0, 1\]としてnormal vectorをパラメータとして計算
clc
clear
cls

format long
%[text] ## satelltie configuration
nFacet = 1;
sat.pos = [0 0 0];
thetaN = deg2rad(30);
sat.normal = [0 sin(thetaN) cos(thetaN)];
sat.area = 1;

sat.qlb = [0, 0, 0, 1];
sat.mCT = 0.2 .* ones(nFacet, 1); % roughness
sat.F0 = 0.5 .* ones(nFacet, 1);
sat.rho = 0.5 .* ones(nFacet, 1);
%[text] ## parameters
const = orbitConst;
d = au2km(1.0, const) * 10^3; % m, distance from sat to sun
dAU = km2AU(d ./ 10^3, const); % AU
S0 = const.S0; % Solar constant, W/m^2
c = const.c; % light speed, m/s
coeff = -S0 / c / dAU^2;
%%
%[text] ## diffuse (analytic)
%[text] 太陽定数などの係数を含む値
%[text] ${\\bf f}\_{{\\rm SRP},d} = -\\frac{S\_{0}}{cr^{2}\_{\\rm AU}} A ({\\bf n}^{T}{\\bf s}) \\int \\int ({\\bf n}^{T}{\\bf v})\\frac{\\rho}{\\pi}\\sin{\\theta\_{r}}{\\bf v}\\mathrm{d}\\theta\_{r}\\mathrm{d}\\phi\_{r}$
%[text] sun vector (from satellite to Sun): ${\\bf s} = \[0,0,1\]^T$
%[text] sun-fixed frameで計算
%[text] $ -\\frac{S\_{0}}{cr^{2}\_{\\rm AU}} A ({\\bf n}^{T}{\\bf s}) \\frac{2}{3}\\rho {\\bf n}$
sunB = [0 0 1];

NS = sat.normal * sunB';
sunlitFlag = (sat.normal * sunB' > 0); % nFacetx1 matrix, 1: sunlit, 0: shade

[sat, srpCdOut, ~] = srpCT(sat, sunB, d, const);
cdTrue = srpCdOut;
disp('True:'), disp(cdTrue) %[output:1f7952b6]
disp('analytical solution:'), disp(sunlitFlag .* (2/3*sat.rho) * coeff * sat.area * NS .* sat.normal) %[output:8ccb2f84]
%%
%[text] ## specular
%[text] varying normal vector: ${\\bf n} = \[0, \\sin\\theta\_n,\\cos\\theta\_n\]^T$
thetaN = 0:deg2rad(3):deg2rad(89);
%[text] ### SG 1st-order approximation
srpCsApprox = zeros(length(thetaN),3);
srpCsTrue = zeros(length(thetaN),3);

for i = 1:length(thetaN) %[output:group:67e0223d]
    sat.normal = [0.0, sin(thetaN(i)), cos(thetaN(i))]; % normal vector

    srpCsApprox(i,:) = srpApproxCT(sat, thetaN(i), sunB, d, const);
%[text] ### true (numerically)    
    [~, ~, srpCsOut] = srpCTuni(sat, sunB, d, const, 'Gauss'); %[output:753b28a1] %[output:80a21f0f] %[output:397e1383] %[output:9330974e] %[output:46873fa4] %[output:86127bd9] %[output:15ac5d70] %[output:5d145350] %[output:67cd6f4f] %[output:659ef3fc] %[output:5c58f371] %[output:65665a33] %[output:884956a5] %[output:757c5703] %[output:2ceaaf4a] %[output:23a7e604] %[output:6d874e92] %[output:0b991242] %[output:7874fba6] %[output:428f3d3b] %[output:3ad0a45e] %[output:06ce9531] %[output:328b97a8] %[output:34d9bdbd] %[output:3cee42b4] %[output:2ad5afb5] %[output:80f00f9c] %[output:3c531833] %[output:1d0b207f] %[output:367ffe0a]
    srpCsTrue(i,:) = srpCsOut;
end %[output:group:67e0223d]
%[text] ## show figures
figure %[output:636f2e98]
tiledlayout(3,1),nexttile %[output:636f2e98]
plot(rad2deg(thetaN), srpCsTrue(:,1),'r'), hold on %[output:636f2e98]
plot(rad2deg(thetaN), srpCsApprox(:,1),'g') %[output:636f2e98]
legend('true','approx.') %[output:636f2e98]
nexttile %[output:636f2e98]
plot(rad2deg(thetaN), srpCsTrue(:,2),'r'), hold on %[output:636f2e98]
plot(rad2deg(thetaN), srpCsApprox(:,2),'g') %[output:636f2e98]
nexttile %[output:636f2e98]
plot(rad2deg(thetaN), srpCsTrue(:,3),'r'), hold on %[output:636f2e98]
plot(rad2deg(thetaN), srpCsApprox(:,3),'g') %[output:636f2e98]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":41.5}
%---
%[output:1f7952b6]
%   data: {"dataType":"text","outputData":{"text":"True:\n   1.0e-05 *\n\n                   0  -0.065333891362464  -0.113161619295973\n\n","truncated":false}}
%---
%[output:8ccb2f84]
%   data: {"dataType":"text","outputData":{"text":"analytical solution:\n   1.0e-05 *\n\n                   0  -0.065333891362464  -0.113161619295973\n\n","truncated":false}}
%---
%[output:753b28a1]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:80a21f0f]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:397e1383]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:9330974e]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:46873fa4]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:86127bd9]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:15ac5d70]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:5d145350]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:67cd6f4f]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:659ef3fc]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:5c58f371]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:65665a33]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:884956a5]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:757c5703]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:2ceaaf4a]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:23a7e604]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:6d874e92]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:0b991242]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:7874fba6]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:428f3d3b]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:3ad0a45e]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:06ce9531]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:328b97a8]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:34d9bdbd]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:3cee42b4]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:2ad5afb5]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:80f00f9c]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:3c531833]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:1d0b207f]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:367ffe0a]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["1","1000000"]]}}
%---
%[output:636f2e98]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAaAAAAD7CAYAAAAsCEDZAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnQ2QFEWa\/l\/wC1YEYVFYvEUXAoMLRSU4hPDAGcQ7L2BX3eUQlhPXUE8iUEI4hf8qejPDAXIsCLtzsqwgLCjfgoZ6SoiwM7phMAinIqLgIciHoYffX4wC8o8np7Op7qnuququ6q7MetIwhpnOzqr8ZVY+9Wa++WaLEydOnBAmEiABEiABEigxgRYUoBIT5+VIgARIgAQUAQoQOwIJkAAJkEBZCFCAyoKdFyUBEiABEqAAsQ+QAAmQAAmUhQAFqCzYeVESIAESIAEKEPsACZAACZBAWQhQgMqCnRclARIgARKgALEPkAAJkAAJlIUABags2HlREiABEiABChD7AAmQAAmQQFkIUIDKgp0XJQESIAESoACF3AeOHTsmp5xyivzwww\/qpzMdP35cjhw5Im3atAn5qiyOBEiABMwjQAEKuc1GjBgh7733ntx9990ycuTIdOmrVq2S3\/\/+93LhhRdKhw4dZMKECfLTn\/405KuzOBIgARIwhwAFSERgtZx66qnpn27N9\/XXXyurpm3btp6tO3fuXOnevbv84he\/UHn37t0rd9xxh6xdu1Zat24tf\/nLX+Txxx+XRx991LMsZiABEiABWwkkXoAOHTokAwcOVO07ePBgWbBgQUZbf\/rppzJ27FjZs2ePfPLJJ9K\/f39lyZxzzjk5+wQECJbOkCFDVB6Izccffyzjx49Xv+MEjKuvvlo2btxoa79ivUiABEjAk0DiBQiEsDazYsUKaWhokNra2gxomFLr0aOHTJkyRY4ePSpjxoyRzp07y4wZM3wL0OzZs+VnP\/uZ\/OpXv0p\/55e\/\/KU88sgjeYXMs\/WYgQRIgAQMJkABSjUerJQdO3ZkCIu2jrZu3arWbZCwvgPr5e2335bJkyfLGWecoYTp\/vvvT0\/PZVtAf\/zjH+UnP\/mJXH\/99emuAmtr5cqVFCCDHx7eOgmQQHEErBeg3bt3S8eOHdMCAlybN29WU2nO5CZA69evF1gvGzZsyMjbt29fWbNmjbz\/\/vty2mmnKQuqX79+cvrpp6t82QKEct555530FBzWnDBF97\/\/+7\/SsmXL4lqQ3yYBEiABQwlYL0BY03nyySfVFFu7du1kzpw5au1l3bp1acFA27kJEP72yiuvyLx585pNy40bN04GDBjg2uzZTgiwpG666SZ54oknpH379gKPOEz3PfTQQ4Z2G942CZAACRRPwHoBAqKZM2cqIYGzAcRn2bJlSgi8LKAlS5bIrl27ZPr06Rl54V49evRoGTp0qGsLuK35vPDCC7J8+XKBGP34xz9Wa035HBmKb1qWQAIkQALxJpAIAUITDB8+XLZt2yabNm2SCy64oFmruFlA9fX1yjLK9ozDFNz8+fOlT58+gVv3u+++U+tGTCRAAiSQdAKJECBYMkuXLhUIx759+2Tx4sVqP46XBbR\/\/35l5UC49PrO4cOH1XrPli1b1NoSEwmQAAmQQGEErBcgrLc8\/PDDsnr1aunUqZM88MADcuDAAVm4cKFyINAJArV9+3aZNWtW+m\/YeDps2DAlXBMnTlT7gKZNm6YcB7AXiIkESIAESKBwAtYL0BtvvKHWWrp06aIoQVTglaY3iWp0ixYtkp07d2YIED47ePCgTJo0SXnOIeF72APEeG6Fdzp+kwRIgARAwHoBCquZEUQUwUX1VFxY5bIcEiABEkgqAQpQUlue9SYBEiCBMhOgAJW5AXh5EiABEkgqAQpQUlue9SYBEiCBMhMwVoC++OILadGiRc7jEUaNGpV2HCgz47Jf\/r29e+X37dvL788+u+z3Eocb6N\/YKJtbtYrDrfAeSCAyAgg3hs3vcU7GCRAiE+AwN+zHaWxslHPPPVd5rvXu3TuDc7du3VTg0KSmdP3r6kQGDRKprhapqkoMjpztDxZIlZXW8+AzwDEg7mOgcQKEiNIXXXSRTJ06VUWhRmw3hLnJPluHD1\/q4cOAC+HBzxMnki1ANTWZ9bdckPkMUIAoQCEOeTjU7fLLL1dTa7B8kHDaKI42eOuttzKiGyT94cNG2bvuukukRYsm4dGDr+WDru5u6frrP8ASBIO\/\/KXpLwng0YxBiM+iCUUlvf4mjIFGWUCweHCsQa9evdL9\/5lnnlHx2hDxwJlMgB\/5Q5w96GoxivzCMbtANoeECFDMWoG3U2ICJoyBRglQdvs9\/\/zzcscdd6iwOldddRUFKBuQnn7DekepB10M+vq6JX7wml0um4POkAArqNzoef3yEaAARcQeUQkefPBBgQDhTB0cs5CdAF8nTEWp6aikJTeLJ2orSA\/qECBM95VbhPKJjKkClEtQk9a\/Wd9mBDDt6IxTyTWgkDvJiRMn5Pbbb5fvv\/9eeb\/lOlPHBPUPGU1mcW7TTlFaQRjM4WmH\/52pnGtOuRg47880EdLig\/tOgCdfpM+I5YWbMAYaNwWHoKE4pyc7mrWbBRR39Y+0\/5di2kkP3lp4nGLjZ\/CPFIA0ef75scKitgrDqmc5p1TDqgPLKRkBClDIqPVZPNhcBVdsJES3\/vbbb9XZPM7jFUyAHzKek8V5Df7FDrjZ1k4uK6fY6xQDKIhlEyRvMfdUzHdz3SP+jvbW3n3FXIPftYqACWOgURYQjtUeM2ZMs07yzTffyNNPPy0XX3xx+jMT4EfW270EqJgBN8h3g+QNE0Yh1y2nWHrV3as+ur39WHte1+Ln1hAwYQw0SoCC9AwT4AepT6C8fqaeChlwvQbC7Jv0EsIglUJZSF5ODYVeM2jdgtx7MXn93lcMRIjhr4ppaP\/f9Rtix4QxkALkv93NyOl3APY7sOla+y03m1IhQpddhg6fg7873bv11J9TlPyIb66WDONew+wlhTBH\/cvknGDCgBdm85SrLL+c\/eYrVz1wXWMF6Ouvv1brP23btnXlZwL8SBre76DlN5++yUIH56BC58eK0tZQdmgdLUSFet4Ve69hNmjQ9nFeu0z1SOwzF2a7+yjLL2e\/+XxcMrIsxgnQp59+KmPHjpU9e\/bIJ598IjBH4fee7Y5tAvxIWjWIBeB3oApSph8BCVJxP9fWguRnis7r2oUKrVe5QT\/3U+98ZWpHEW0NeU1det2fj43Fpj9zzz33nJx99tlyxRVXeNEo6+d+OfvNV87KGCdAI0aMkB49esiUKVNUMFI4JXTu3FlmzJiRwdEE+KE3fNC3Zj8C5CePV0WKGdSL+a7Xfbl97re+YF1fH40HWrHi46yX7hP4W1AnBV1HuNlDwLQI4d8VFc3W40x+5rC5vba2Vrp27SrXXXediiuJGZaWLVsKPsO\/zzzzTME+RBwDo5POo3\/\/v\/\/7P2nXrp2cccYZhfQ+X9\/xy9lvPl8XjSiTUQJ06NAhFfVg69at0qFDB4UEe32uvvpq2b17t5x66qlpTCbAD71NgwoQbiDfAF9IecUM6mFbT4UCzsdEuz1nW1uFTvtl36NfAQxaN79ClC062etJTuHNWo\/rdsstxh6B8stf\/lIOHjyotnJAaJ599ln527\/9W0H0fQQ\/hlX0q1\/9SubPny9Lly5V9LEtBCHA3nzzTRUUGWHBsB1k3759MnHiRPnFL34RtJV85fc7tvnN5+uiEWUySoDWr18vs2fPlg0bNmTg6Nu3r6xdu1a9vehkAvzQ27SQN+dcA15Y4oNKFlpWIfUJA2o2k+wNt9lv\/2GJRljl5GOg20KLhxZOWHPa0gnixOBYj+v2\/vvGChCQzZ07V\/7mb\/5G\/vmf\/1kRxBgyadIkue2229TL7euvvy7z5s2TRx55RH2Odeh\/\/Md\/FGwPqampkSuvvFIGDRokX375pfzTP\/2TChUGayjs5Hds85sv7PsLUp5RAoSo12hsdAJnwrTcuHHjZMCAAek\/d\/h9k4WE1K9fP7VWZH2qwaFzWaFw\/FTa7Xt+3Z79lI88hdxbId\/xez9e+fS18VOnfGzDuNcwyvCql\/NzZc2k3Nvx90L6jqO8Jd2WnBQgvQYV5H5KmdflgEY3AXr77bfT02lvvPGGPPzww2kBwv7Df\/iHf5CXX35ZjS8VeDFJJURr+fd\/\/\/dIrKB8wsJYcBF2oiVLlghORJ0+fXrGVUaOHCmjR4+WoUOHZghQ5AFI9SCtH+KKVNTpYhd8C2FYjGBkfzeKgTDo\/UVxD4VwDfKdQu85KJsg91TCvBkCVMLrhnWpbAHCsS+YXtMJAvRf\/\/VfsmDBAvUnTLv9y7\/8i\/z1r39VL78rV67MuBUsE7Rp0yas20uX49ey8Zsv9BsMUKBRFhDeKmAF6Q6g64kpOMzN9unTJ3Aj+WKl572RWQ8Wzvlvp+Doz4Mu+Pq6kTyZipmuck79RDUNFLTcUjsfFMsf3w9aR33NYtoujPsOqQwTBrx8VYUTQvv27QUzKlgLyq6PXvPZsmWLmpK788475X\/+53\/k1Vdfld\/+9rcqEgteho8fPy433nijcozq3r17SHRPFuOXs998od9ggAKNEqD9+\/crK2fbtm1y+umnq2rq+HDoFFgA1Klo+M55fwiMFhltZuezcvwu+AZoqLxZC11jcRaqNzCirKjiivkVlTDqExbboOUEFSFLxAeYin7mgrIOOT\/GlZtvvlmJ0FNPPaWm1F577bUM5yY4ICAQMvYf\/vrXv5Y\/\/elP8tJLL6lxCEfEwGHhggsuUA4Lei0p5Nv0zdmE9jBKgODyOGzYMIHFAy8T7AOaNm2acpV0noFR0MOQ7d2jjxUoxrvJbcE3ium5MAZsPWcP8YniHoNYCKYPyn6FNqhYhT2ShVyeCQOeV5WPHTuWIThu+fPlgWNCFNNuzvvwy9lvPi8mUX5ulAABBFwl4ZmCNw2kIUOGKFM3u9F9wXezclz2N4TSAFFaRaYM2H4G3DDENJQGK6IQP\/X0k6eIWyjHV309c+W4Mcuu6Zez33zlxGOcAGlY2Bx2yimnpKfisiF6wneGcSnGygnaemELkWkDtpd1YIqYerV7PoExrc286pr63POZ81kOs+Un4Jez33zl5G2sAHlBywlfTzW5uGF6lRnq52EJkWmDmdebv5dAhdoIERfmVhfT2isAIhMGvADViW1Wv5z95itnRY0VoC+++EKFxAgUjNRr8CtHSxQrRKZZDEmyDNzqalp7BXgmTBjwAlQntln9cvabr5wVNU6AsA9owoQJyuuksbFRzj33XJk1a5b07t07g2MGfD3IB9nhXepWKUSITH2bzmXl2Dg4O0XIxvo5nhMTBjyvxxrxJRHyCxERdGgvHe\/tu+++E\/zvfOn185kzlhyujw2sX331lYphqVN2TLns35337Zez33xeTKL83DgBQmwmHMc9depUFYx0zpw58sILL8jGjRvdBUjH7ir1vpxCW83pOedVhnYPL+Ualtc9+fnczTIwVUz91Fe7uCOvaW3lp36pPCYMeLmqA8+2\/\/iP\/1BbPLp06SIffPCBLFq0SMWcxMvtDTfcoELrIMHFGnHfsFco32cXXnhhRiw5lI\/9QojmgqjbcNfG3xCu5+c\/\/7kK54MNrdjnuGPHjmaevfre\/XL2my9AE4ee1SgB+vjjj+Xyyy9XHnCwfJCwG3nw4MHy1ltvqQi2GY10\/vllO5wr9JZyK9C5GbYkFwzpIgmbmkpvXo7KvT2kZim2GOeAVyM1Ui0FhIUq9iZ8fh\/3ViVV6dwHDhyQ1atXy7\/927+pqf2ZM2eq\/TbYy4OfiAd3zz33qGjY8MJFAGQEG\/X6zBlL7o9\/\/KO888478rvf\/U45T2FPEeLLPfTQQ4KQP9jYitkcXOfJJ58MtrzgUm8KkM\/O4DcbLB40IEJk6PTMM8+o6AirVq3KLAbTPFHuafF708znTiB7Gs4m54OEtrkJA16+pkGonf\/+7\/9W0awxztx6663ym9\/8RokMorD89Kc\/VV9\/8cUXBWcHQTi8PnPGkkNZEDJE9EdC0NK\/\/\/u\/VyIEj15YXJjZWbFihYpfmSv55ew3Xzm7q1EWUDYomMQwhbEzGWHRnQnwdUJMuMjjwpWzFU28ttMKsnn6zcS2KfCeTRjwclUN4gPrB3EmMcW\/fPlyFYRUCxBEQq\/9YNxBAFLkRZ1zfZYdSw5jEKwmWE9IOFzz7\/7u7+Tdd99VAgTLBxFdIFI33XRTQQLEYKQFdl63r6ExMM+KhAVBbflgUQ9hL9AR8Bai3yiyBQhnBSU1oSPGWnQjXpyPff1L0DFLzcBkAcIMCqyeqqoq+eyzz1SQUcSE0wIEUcDzhLUiCBVC9mDc0VNwuT5zjkGYVsM5QxizfvSjH6n4lVjrQTgfhP7B53\/4wx+USOFvWENCrLlLLrkk44A7v5z95itBV8x5iVhbQIhygDA7SD179lTzo5iDvf322+X7779Xv2cfxa1ragL8KBveiPrraVKIUcjx54yof5QdoAyx2UxmDq9ajCuY5sdaMqbbcCDdv\/7rv6rQX7fccouaGsPLL+JRTp48Wb0c5\/oM60hwUHDGkoN4LV68WIlMq1atlKWFKC4IXoq4cij\/vPPOU9YVRAp5IXJwsvrJT36S7i1+OfvNF3E3zFt8rAXI7c4xT4r5WEy7aevILZ8J8KNseCPqH2H8OSPqH2UHoAAVRNctzptzKi37SO58n+W6AVwDbtY6oHLQG\/Xbt\/3mC3r9MPMbJUA68jXmZ\/H2gISG\/Pbbb1UkbKcgjRo1Kh0vLkxgLCs8And9\/rnc9dln0u1nPwuvUJZUVgK2TXvDOunRo4fs2bNHecc5U77PomwE5\/p2vuvgkDyMlXFORgkQ\/OfHjBnTjCc2dj399NPqPA4mEiCB8hDgS19puJsgLH5JGCVAfivFfCRAAiRAAvEnQAGKfxvxDkmABEjASgJWChAOhcLaUK5Apba1pFdg1qTxcGtf2xkgRhnchzt16tRsrQI8bK8\/6oiNnXASQGgb2\/uAdpbwOhwv3zgYhz5hlQBhY9fYsWPVgiFOS8VcKfZC5HLVNl2IvAKzJo3Hpk2b1CY+7OfQHka2M4DbcHV1tVoDxVHScMRxBue1vf54hrERtLa2VrZv354OUIx4anqx3kYGV1xxhXz44YdqCMt2\/PCqr9fnpRwXrRIgbByDx8qUKVOUPz8cFhBxFr72NiavwKxJ4gEPSewbw4uH84G0nQGC8SIszBNPPKHe\/OfOnSsIT7Vhwwa1h872+uO5vvLKK+Xaa6+V8ePHq+cem0nRDx599FH12NvIAB54u3fvVkFS33zzzYzhzau+Xp+Xcqy0RoAQQh2btrZu3aoi2Oo3A4S9QEPp0OqlhBvltbwCs+ItJyk89Obks846S+0o1wJke5\/AIKStfLwR62koPAODBg1SEZ1t7wPwgMVeHGzexCZOJES0RqQCDMw29wHEmUO8OngH6+RV348++ihWfcIaAVq\/fr3Mnj1bvfk5E3Yqr127Vrp27RqlHpS8bK\/ArEnisWzZMvXmD\/HBgKsFyHYGiASPN39MQe3fv18+\/\/xz6d69e3rt0\/b664cOzziiV0N0kVauXKme+TVr1ojNDNwEyKu+O3fujNU4aY0AISI23gTmzZvXzBwdN26cOmfD5pQdmDUpPLDeh6lIhLZHuCZsUNYCZDsDHEuC6WWcLYM3XzgiIFQMpuXQ322vv36eERkFgTwrKioUA3BBpJRLL73UagZuAuTV5oj0Hadx0hoBWrJkiWBRHhFqnWnkyJEyevRoFb\/JxpQrMGsSeKDuw4cPV\/PgiB6M350CZDsDHAuA+GWICH\/33XcrDzAEu7zvvvuUVYTBKAnPBF46ddw0xIhsaGhQgUNvvPFGsbkPuAmQV30xdR+nPmGNAOEtCA8cvF+yp+AQdbZPnz7W6U++wKxJ4IGXDQRwRGBIhEnB2y88wvA\/nE8QTt\/mPoE1juuuu04NuNrTE265CIIJLliIt7n+eKDxNg+hwboPvACRcLzBNddco87twUFztjJwEyCv5x6u13HiYY0AYQ4cVg46onbB1bHjcKwDYsXZlvIFZk0CD3h+YbpNR0zH2y+mXuD9eP755wsW5m3uE9j3gxcrDLTa5RgvJQjfj9M9zzzzTKvrj+cZU6\/w+sN6jzOh7eENiyMNbO0DbgLk9dwjbmaceFgjQNhwNWzYMBUefeLEiertb9q0aWpwwl4g25JXYFYccJUkHmhfbES87LLL0mtASegTmG5Dve+99161CRVTUY899pgSJXh+2t4HcJ4OHDHgco11L+yD0uuhL730knTp0sVaBhAgTD2\/+uqr6eHNq897fV7qcdIaAQK4gwcPqvPasQiJhH0hWKRt06ZNqblGfj0\/gVmTxAPAIco4GRduyJh+S0KfgLs9jnGGByASLCK8eOHNPwn1Rx3XrVunDqhsbGwUuOIjKDHWBfXJo7Y+BzhrCOt\/TjdsP20eJx5WCZAe9bEYDQug0PM2IlePEl+APEQ5KNjcJ7D+hbdbHKbmlmyvP+qM9Q1YQPrlI5tDEhg46+xVX6\/PSzFMWSlApQDHa5AACZAACRRHgAJUHD9+mwRIgARIoEACFKACwfFrJEACJEACxRGgABXHj98mARIgARIokAAFqEBw\/BoJkAAJkEBxBChAxfHjt0mABEiABAokQAEqEBy\/RgIkQAIkUBwBClBx\/PhtEiABEiCBAglQgAoEx6+RAAmQAAkUR4ACVBw\/fpsESIAESKBAAhSgAsHxayRAAiRAAsURoAAVxy\/j24itlB2L6\/jx4yoGGRMJkAAJkEAmAQpQiD0CR0PjQDCdcP46DkdDyHQmEiABEiABClCzPgDRwNkp+qdbJ0GkXUQbbtu2ra8+tHfvXpkwYYI6mTJXhGJfBTETCZAACVhKIPEW0KFDh2TgwIGqeQcPHtzsSG+ctzJ27FjZs2ePOuSuf\/\/+6oA7fQRyrn4xatQoufPOO9WpnEwkQAIkQALNCSRegIAE6zSwVBoaGqS2tjaD0ogRI6RHjx7qeN+jR4+q4547d+6sDrrLlXBA1COPPCJ\/\/vOf2edIgARIgARyEKAApcA8\/vjjguN9ncKirSOcsNmhQweV87333lMnLeI43MmTJ6vDryBM999\/f3p6DpYP1oP0iYzsfSRAAiRAAgm0gHbv3i0dO3ZMCwgQ4MhuTKU5k5sArV+\/XmbPni0bNmzIyNu3b19Zs2aNvP\/+++oERlhQ\/fr1UyewfvXVV2rabdu2bTyRlU8cCZAACeQhYL0FtGDBAnnyySfVFFu7du1kzpw5snHjRnWOvPPIbjcBwt8wnTZv3rxm03Ljxo2TAQMGNEPL6Tc+byRAAiTgj4D1AgQMM2fOVEICZwOIz7Jly6R9+\/aeFtCSJUtk165dMn369Iy8I0eOlNGjR8vQoUObUV66dKnAcWH8+PH+WoC5SIAESCChBBIhQGjb4cOHq2mxTZs2yQUXXNCsud0soPr6esHfYUU5E6bg5s+fL3369Elot2G1SYAESKB4AokQIFgysEwgHPv27ZPFixc325vjJkD79+9XVo5zPefw4cNqvWfLli1qbYmJBEiABEigMALGCpDfjaGrVq2Shx9+WFavXi2dOnWSBx54QA4cOCALFy5UDgQ6QaC2b98us2bNSv8NG0+HDRumhGvixIlqH9C0adOkZcuWai8QEwmQAAmQQOEEjBOgoBtD33jjDbVptEuXLooSRAXebUOGDMmgtmjRIkHoHKcAIcPBgwdl0qRJynMOCd+Dq3abNm0Kp85vkgAJkAAJiHECVMjG0DDaGYFGEVTU6TkXRrksgwRIgASSSsAoAcq3MRT7fRDPjYkESIAESMAMAkYJUL6NoWvXrpWuXbuaQZ13SQIkQAIkYNYUXJCNoXVSJ5VSySYmARIgARKIKQGjLKAgG0Nbb26tkLf\/Q3tptblVTPHztkiABOJMoH9jo7q9fqmf6t9Hjsgf2reXza3iPa4g3Njy5cvjjNcsCyjIxtBu3brJovcWSY3UqAaokqpEWUSoPwKnJjUlvf5o96Qz8FX\/ujqR+vqTjwl+R9I\/K1OzKPonPquoEHH+HtOHzFf9y3zvRllAQTaGOuFjOi5pQmRC54uy7ye9\/hQghwC7iYxTYCorBWMEUn1FqldmCYz+XPdZE15oTXgGjBKgIBtD3eAnSYiwUfauu+6KcoyPddlJrz8aJ3EMsoRm84wZ0ti\/UeodS8F1DoFxiopeL861blwh+otN3d6E9WUKUARDlN+NofngJ0mIImgCFkkC5SdQ0zS1rmbLpE7qUxaM+r1CpE7PnCmpOKlATiExQUSKAU0BKoaex3e9Nob6gU8hirCBWDQJhEHAYdU4hcZLZGwXFz9o\/YyBfsqJMo9RU3BBQASBnyFEdZVSWVkV5FLMSwIkEBaBlGVTUwHngLqc1gwsGYpMfuhBxsCwmi9oORQgB7G0ENXVSZVUU4iC9ibmJ4EgBFLWDayZ+vpq9c3q1LtftTT9TqEJAjQzLwWocHZFf7MY+BSiovGzABJoTqCuaRJNrdc4rJumVZqmdRp4lzGFQ6CYMTCcO\/AuhRZQHkYZQlRfKZUVVUb4\/3s3O3OQQAkIZAkOrBun2NC6ibYNKEDR8s1bepjwKURlbEhe2hwCDsHBMwNPNFo35Wu+MMfAqGpBCygAWQpRAFjMmggCdXU1TftsaqrV+o0WHFo35W9+ClBEbXDs2DF19IL+6XaZKOFnCFENZuWqjQnPEVGTsNikEKipkRos02QJDtdv4tcBohwDw6qtcRaQPhMIAAYPHiwLFixwZVEK+M2ECAupVVwnCqtzspwYEKirk5rKehUbrbqyKVyN9lCjw0AM2ifPLZRiDCyWgHEChAofP35cVqxYIQ0NDVJbW1s2AdIXphAV2w35\/TgR0NNqdXXVah2HghOn1vF\/LxQg\/6wC58TZQDt27JAZM2aUXYCyhQiCVF0jUoWnlxZR4LblF0pPQAfrhdhwHaf0\/KO4IgUoCqqpMuMoQM7q4oFWD3OdSBXWiTg9F2FvYNFBCegdOdq6oZUTlGD881OAImwjPwKkL4+o0OWKDM3puQg7AYsORMDNyqmcr4AmAAAbd0lEQVRCGAIDzrYJVNEEZ0YEdPyvU9zPBIv1GtCWLVvktNNOUyzh9darV680WD8CFCf4rkKEB9+Qw60S\/EwbW3Wu5RjbdKHcOC2gIjEOGTJEWrZsqUrp2bOnzJo1y1gB0jfebC9R6gBGrhUV2Vn4dUWAaznsCJoABSjCvmCaBZSNQgsRfqbXiSBGsIrouBBhz7Gs6Cw3aa7lWNa+RVSHAlQEPK+vLl26VLZv355hFTm\/YwL8ZlYRgjHiOIialFlEIfLqBon8nFNriWz2wJU2YQyM9RpQPuKLFi2SnTt3WiFAuYQIu83V3j9aRYEfPqu+kCP6AMPdWNXKoVeGAhQ6Uv8FmgA\/V22cB+QhT9oqwvkpFCP\/ncDUnKlzcvShbIixxqk1UxuzfPdtwhhorAXk1awmwPeqAz7PWCuSSk7R+YFmYp6aGtXWOCsn+1A2hrwxsUHLf88mjIEUoPL3E9930OzocKwVOa0ilMQ9Hb55ljWjtnJSJ386I0lTcMraMtZcnAJUxqY0AX6hePJaRbpQOjAUije672Wt5eBCnFqLDnfSSzZhDKQFZHgvdRMjZQjROipvy6YsHAQa0EdQcy2nvE2StKtTgMrY4ibADxtPM+cFqRKpqz\/p1u20jjhdFy7+mpqm8tS+nNSxBVlHUHNqLVzkLC0\/ARPGQGMtoC+++EJatGghbdu2dW0FE+BH+QBlixGCoVbAiU4qRByDpVozYkigYE2Rsm7Ul6qr4S3f9E8c1MZptWAsmTsyAiaMgcYJ0K5du2TChAly+PBhaWxslHPPPVftBerdu3dGQ5oAP7Kel1WwjnyMnyryArzpYB1hqg4v6\/VNB46lHRq0dYQ4dUm3lJxiXVenzsdRR1BTcErVfXmdAgmYMAYaJ0DXX3+9XHTRRTJ16lQ5evSozJkzR1544QXZuHEjBchnR3VaR66ChHK0KOHf2tMO\/4Zzg22ihPrpOuv6pv6mrRus5UB89Fk5CkVKxH1iZzYSKCkBClDIuD\/++GO5\/PLLZfPmzcryQdq7d686mvutt96S1q1bp69oAvyQ8WQUh5Dsfo+gcBMk5wCrzjLSg7TTInAKkZ7Gi4k4Nau\/m8g4xFVZiSnLRouNsnJSbtImRh0I0gei7IvlKjvp9TdhDDTKAoLF884772Qcy\/DMM88IApOuWrWKFpCDQDGdD4Mxko6srK2kZqKkB3BtPehB3mkx6Xty7k\/SU3tarPT3go5UsNLcUl2deknp39h48tPU9WHFqNtNzS7id10\/JbQiYqLYuGEopg8EbYo45mf9u0mcjqRx6yNGCVB2BZ5\/\/nm54447ZOHChXLVVVdRgEISINcx3UWUlH6kBm39bwze6t\/6mAldmFMsnILjJlZ+RzOXTbdKPCsq5T8bZkj\/\/\/dbJS5KJ1M\/9f3aJjYUoOYEKEAUIL9DiWu+XAfSHTlyRB588EGBAD300EMycODAZt8fNWqUegtmipZAY\/+TVkZjv0Y50u+IuqDz7\/i91eZWzW6kdcPJKVO\/d6nL1\/lzXUeX3aqh6bpu1\/d7TeYjARMJ9O\/fX5YvXx7rW4+1BeR2IN2JEyfk9ttvl++\/\/155v51zzjmxBsybayKgLRAnj3qpV9Nd+Ok3aQtL53daYH7LYD4SIIF4EIi1ALkhwjEM9fX1atpNH9cdD5S8CxIgARIggSAEjBIg7P3p16+fMivhio30ww8\/yLfffisdO3akIAVpeeYlARIggTITMEqAXnnlFRkzZkwzZN988408\/fTTcvHFF6vPvv76ayVMuaIklJl56Jf3igqRNB5ugG1n8N1338lnn30mnTp1UhFCspPt9Ud9v\/zyS8EUfbt27VyfMZsYHDt2TE499VTRPwvp83HgYZQAeY3cn376qYwdO1b27Nkjn3zyiWARDnsBbF0n8ooKkTQemzZtkttuu0256p9++umqu9jOAFsTqqur1QtY+\/bt1SyAMzKI7fVHG7\/++utSW1sr27dvT0dHWbBggcALztY+cMUVV8iHH36o6pftau3V5l6fe42zYX5ulQCNGDFCevToIVOmTFFREmAtde7cWWbMmBEms9iU5RUVIkk8MD0LpxW8eDgfSNsZIBLIc889J0888YR68587d65gb9yGDRukZcuWYnv98TBeeeWVcu2118r48ePVc19VVaX6waOPPqqeVRsZHD9+XHbv3i033HCDvPnmmxljkld9vT4v5QBnjQAdOnRIuWNv3bpVOnTokH4zuPrqq1VDwVy1KXlFhcBbTlJ4aM\/Is846S5566qm0ANneJzAIaSsfb8R6GgrPwKBBg+SDDz6wvg9g+r1Xr17y8ssvy3nnnacYbNu2TW6++WY1MNvcB95++2259dZbBUsTOnnV96OPPopVn7BGgNavXy+zZ89Wb37O1LdvX1m7dq107drVJv1Rb3r5okIkiceyZcvUmz\/EB6KrLSDbGSAMFd78MQW1f\/9++fzzz6V79+7ptU\/b668faDzjM2fOVKKLtHLlSvXMr1mzRmxm4CZAXvXduXNnrMZJawQI4XjwJjBv3rxm5ui4ceNkwIABVglQdmWyo0IkhQfW+zAVuXTpUunZs6fyjtQCZDsDbLTG9PLZZ5+t3vThiIBN2piWQ3+3vf76GcC2jHvuuUcqKioUA3DBNo1LL73UagZuAuTV5vv27YvVOGmNAC1ZskSwKD99+vSMsXnkyJEyevRoGTp0qJUClCsqRBJ4oO7Dhw9X8+A33XSTGnydAmQ7gxdffFFtykY4qrvvvlt5gD377LNy3333KasIg1ESngm8dC5evFhZv9ig3tDQoALx3njjjWJzH3ATIK\/6Yuo+Tn3CGgHCWxAeOHi\/ZE\/BzZ8\/X\/r06WOdAOWLCpEEHnjZWLFihUyePFm5HuPtFx5h+B\/OJ2eccYbVfQJrHNddd50acLWnJ9xycTYWuGAh3vZnArMeEBqs+8ALEOndd9+Va665RiDQBw4csJaBmwB5PfdwvY5Tn7BGgDAHDisHHVG74OqNq4gph42qtqV8USGSwAOeX5hug7cXEt5+MfUC78fzzz9fsDBvc5\/Avh+8WGGg1S7HeCm55JJLZPXq1XLmmWdaXX+0OaZe4fWH9R5nQtvDG\/bCCy+0loGbAHk999i0H6dnwhoBwsbTYcOGCRYkJ06cqN7+pk2bpgYn7AWyLXlFhTjllFMSxUN7gF122WXpNaAk9AlMt2ED5r333qs2oWIq6rHHHlOiBM9P25+JHTt2KEcMuFxj3Qv7oPR66EsvvSRdunSxlgEECFPPr776anp48+rzXp+Xepy0RoAA7uDBgzJp0qR0FGzsC8EibZs2bUrNNfLr+YkKkSQeAA5RxrEccEPG9FsS+gTc7XE6MDwAkWAR4cULb\/5JqD\/quG7dOhUdv7GxUeCKj4goWBfEFgybGbz22mtq\/c\/phu2nvnEaF6wSID3qYzEaFoCeiotcDWJ+AfIQ5aBgc5\/A+hfebp2nAju7pe31R12xvgELSL98ZD+WSWAQpM3jwMNKAYq5HvD2SIAESIAERIQCxG5AAiRAAiRQFgIUoLJg50VJgARIgAQoQOwDJEACJEACZSFAASoLdl6UBEiABEiAAsQ+QAIkQAIkUBYCFKCyYOdFSYAESIAEKEDsAyRAAiRAAmUhQAEqC3ZelARIgARIgALEPkACJEACJFAWAhSgsmDnRUmABEiABChA7AMkQAIkQAJlIUABChE7gvtlB4M8fvy4CoLJRAIkQAIkkEmAAhRij7j++usFJ1LqtHPnTnU6J87sYCIBEiABEqAANesDEA0c3qV\/unUShHpHuPu2bdv66kN79+6VCRMmqKORc4XI91UQM5EACZCApQQSbwEdOnRIBg4cqJp38ODBsmDBgoymxoFfY8eOlT179qhTVvv3769OWD3nnHPydolRo0bJnXfeqY6FZiIBEiABEmhOIPECBCRYp4Gl0tDQILW1tRmURowYIT169FDnyx89elTGjBkjnTt3Viet5ko4ofCRRx6RP\/\/5z+xzJEACJEACOQhQgFJgHn\/8ccH58k5h0dYRjnju0KGDyvnee++po35xHvvkyZPV6YsQpvvvvz89PQfLB+tB+khg9j4SIAESIIEEWkC7d++Wjh07pgUECDZv3qym0pzJTYDWr18vs2fPlg0bNmTk7du3r6xZs0bef\/99dQQwLKh+\/fqpI8C\/+uorNe22bds2HgnOJ44ESIAE8hCw3gLCms6TTz6pptjatWsnc+bMkY0bN8q6desyBMJNgPA3TKfNmzev2bTcuHHjZMCAAc3QcvqNzxsJkAAJ+CNgvQABw8yZM5WQwNkA4rNs2TJp3769pwW0ZMkS2bVrl0yfPj0j78iRI2X06NEydOjQZpSXLl0qcFwYP368vxZgLhIgARJIKIFECBDadvjw4WpabNOmTXLBBRc0a243C6i+vl7w92zPOEzBzZ8\/X\/r06ZPQbsNqkwAJkEDxBBIhQLBkYJlAOPbt2yeLFy9utjfHTYD279+vrBznes7hw4fVes+WLVvU2hITCZAACZBAYQSsF6BVq1bJww8\/LKtXr5ZOnTrJAw88IAcOHJCFCxcqBwKdIFDbt2+XWbNmpf+GjafDhg1TwjVx4kS1D2jatGnSsmVLtReIiQRIgARIoHAC1gvQG2+8oTaNdunSRVGCqMC7bciQIRnUFi1aJAid4xQgZDh48KBMmjRJec4h4Xtw1W7Tpk3h1PlNEiABEiABsV6AwmpjBBpFUFG4WjORAAmQAAkUT4ACVDxDlkACJEACJFAAAQpQAdD4FRIgARIggeIJUICKZ8gSSCC+BOrq1L3VVbrcYl29\/\/uurBDJzo+\/iUhl0yWap0q3i\/q\/JHPaT8BYAfI6HgHRqLXjgP3NyBomgUBj\/0ZVzcZ+TT91OtLvSPrf\/RsbXcWm1eZWkSHS95V9gVzC1Koh816Qb3Or1umvt876fHOr6O49MigxKBjhxpYvXx6DO8l9C8YJkN\/jEbp166YChyY1sf7mtH+dNJkQ9XLSItF\/0z\/xeaU0WRT6p+7bFdJkieikPy93H3Deu\/P+nPXE37PzZf\/eTMgqKzOsrooUF0lZZIpRnUi3W27hGBDzMdA4AfJ7PEK5H75yCx\/rHy8BwqCaLTB6oFUDLAZVPZCKiBaVbLEJ0q9s6QNuQuYm1k42zUTMwdaNYVDOVVIVpCnKkteE9jdKgPIdj4Co1zjVVCcT4EfZK7FR9q677oryErEuu1z1dwoN\/p0WGSUvJ9dEIDBaeKICWS4GUdUnaLnp+tfVnZyW1OtY9Y6Fq7o6qc+1XFWR+iBrPUu1n4eoBb3fsPObMAYaJUD5jkdYu3atdO3alQIUdi9mea4E\/AiNCYMUm9dBIOWwof5Sn5oO1X\/TP7UQVVUpqzXOiQIUcusEOR7BBPgh42FxERFwE5sme6ZpAKLQRAQ+jsVmC1Ec7zF1TyaMgUZZQEGORwB8nTAVleTpqBg\/I7G7NafYVEu1uj+KTeyaiTeUgwCmHZ1xKuPuiGWUAAU5HsEE9edTVF4CXpZNVV1F7KdZykuQV48zARPGQKMEKMjxCCbAj3PntfHetOBo5wA4AahptMpKTqPZ2OAJr5MJY6BRAhTkeAQT4Cf8+Yi8+q6CU1ktFUp4aN1E3gC8QFkJmDAGGiVAaE2\/xyOYAL+svdPCizcTnNTqjRKcyvjv27CwSVilMhIwYQw0ToDQnseOHZOjR4\/KiRMn5Ec\/+pFrE5sAv4x904pL5xQcA\/ZoWNEArESsCZgwBhonQHozKlp+8ODBsmDBAgpQrB+DEG+urk5qKutF6uqkurIu7Z1GN+gQGbMoawhQgCJqyuPHj8uKFSukoaFBamtrKUARcS57sVmCg\/vRrtEmhEIpOz\/eQKIJUIAibH5sSt2xY4c6HtstmQA\/QjxmFo2QKYiYVgkjp1qFT6HgmNmUvOvyEzBhDDRuCk43KwWo\/B08lDuoqZEa+AfUVEs1optoxwGu44SCl4UklwAFKMK2pwBFCDfKorXgpNZxOK0WJWyWnWQCFKAiW3\/Lli1y2mmnqVIQ6bpXr17pEv0IkM7MUDxFNkQxX9fTaqnI0JxWKwYmv0sC+QkwFE+IPWTIkCHSsmVLVWLPnj1l1qxZgQQo7nGQQkQVn6IQrLG+ntNq8WkR3klCCdACirDh\/VhAFKAIG0AXnSU4CJ9GK6cE3HkJEvAgQAGKsIssXbpUtm\/fnmEVOS9nAvwI8URbdE2NQGjq66vVT6fgcE9OtOhZOgn4JWDCGGisF9yiRYtk586dFCC\/vbGYfDU16ts1qeMJ6K1WDEx+lwRKQ4ACFCHnL774Qlq0aCFt27Z1vYoJ8CPEU1zRKcHBXhy1J4fTasXx5LdJoAwETBgDjbOAdu3aJRMmTJDDhw9LY2OjnHvuucoK6t27d0YTmwC\/DH3S\/ZJZgoNMtHJi0zq8ERIoiIAJY6BxAnT99dfLRRddJFOnTlUBSefMmSMvvPCCbNy4kQLkt5umBEeqq7H\/UyWn4OB3hrrxC5P5SCCeBChAIbfLxx9\/LJdffrls3rxZWT5Ie\/fuVUFJ33rrLWndunX6iibADxlP7uK04KjYanVpwVHCk1rXoeCUrDV4IRIoCQETxkCjLCBYPO+8807GhtRnnnlG4JK9atUqWkAgkHKLbvpn0xqOtnAoOCV57nkREogFAQpQxM3w\/PPPyx133CELFy6Uq666KnkC5BAbCI8K5ekQHB1XjVNqEXdEFk8CMSRAASqyUXKF4jly5Ig8+OCDAgF66KGHZODAgc2uBPg6WRGKJ0tsIDh6\/UZ7qVFwiuxw\/DoJGE6AoXhCbEC3UDw4BfX222+X77\/\/Xnm\/nXPOOa5XNEH9XW\/cRWiQL1tsnNNp3PwZYqdjUSRgCQETxkCj1oDQL7ABtb6+Xk276UClbv0l9vBzCA2iCuhpNG3ZUGwsGRFYDRIoIYHYj4EiYpQAYe9Pv379ZPny5coVG+mHH36Qb7\/9Vjp27JghSLGA7xQZ3Cx+1z8rK9MeaVJRqcLZYA3HOY1WjGUDUxxTj0lNSa8\/2j3pDJJe\/1iMgR4DkFEC9Morr8iYMWOaVembb76Rp59+Wi6++OL0ZyWBn09goDdwCKislHpJCQ+EBv9pR4GUCzRuuhixMdICjFgZS9L+Edeh2OKTzoD17yZxD8hslAAFeSADdz5tneiL1Nc3vxzypPJpEdHTZSpzlsDAmkHSP6MQmlxMAtc\/CFwD8ia9\/miipDNg\/SlAZRuqWm9uLf0bG9PXr3TRkwocXZOyRnRGrLu4JS04+rNWm1upf7ZuOLn5Fb+3amgl+rOyVZ4XJgESSDyB\/v37q+WKOCdrLSBMdTlTvbgokEvLYCrMLTmtmDg3KO+NBEiABEwhYK0AmdIAvE8SIAESSCoBClBSW571JgESIIEyE7BSgL7++mvlnp3rrKAyMw\/98l5nIyWNhxtg2xl899138tlnn0mnTp3UOVnZyfb6o75ffvmlYKN6u3btXJ8xmxgcO3ZMTj31VNE\/C+nzceBhlQB9+umnMnbsWNmzZ4988skngkU47AXIFS0hdCUocYFeZyMljcemTZvktttuUwFrTz\/9dNUatjNAgN7q6mq1DaF9+\/ZqL5zzfCzb6482fv3116W2tla2b9+ePiNswYIFygvQ1j5wxRVXyIcffqjql+1q7dXmXp+XchizSoBGjBghPXr0kClTpqizgrBnqHPnzjJjxoxSMi3ZtbzORkoSD2xSRugmvHg4H0jbGeA8rOeee06eeOIJ9eY\/d+5cQYT4DRs2SMuWLcX2+uNhu\/LKK+Xaa6+V8ePHq+e+qqpK9YNHH31UPYs2Mjh+\/Ljs3r1bbrjhBnnzzTczxhyv+np9XrIBzLRICPnAHDp0SAUl3bp1q3To0CH9ZnD11VerhoK5alPyOhsJbzlJ4aHjA5511lny1FNPpQXI9j6BQUhb+Xgj1tNQeAYGDRokH3zwgfV9AJvQe\/XqJS+\/\/LKcd955isG2bdvk5ptvVgOzzX3g7bfflltvvVWwQV8nr\/p+9NFHseoT1lhA69evl9mzZ6s3P2fq27evrF27Vrp27WqT\/qg3vXxnIyWJx7Jly9SbP8QHoqstINsZ4DBGvPljCmr\/\/v3y+eefS\/fu3dNrn7bXXz\/QeMZnzpypRBdp5cqV6plfs2aN2MzATYC86rtz585YjZPWCBAOpcObwLx585qZo+PGjZMBAwZYJUDZlck+GykpPLDeh6nIpUuXSs+ePVWMQC1AtjPAycCYXj777LPVmz4cEXBUCabl0N9tr79+BhCc+J577pGKigrFAFwQrPjSSy+1moGbAHm1+b59+2I1TlojQEuWLBEsyk+fPj1jbB45cqSMHj1ahg4daqUA5TobKQk8UPfhw4erefCbbrpJDb5OAbKdwYsvvqiOJsGhjHfffbfyAHv22WflvvvuU1YRBqMkPBN46Vy8eLGyfnFMS0NDgwrEe+ONN4rNfcBNgLzqi6n7OPUJawQIb0F44OD9kj0FN3\/+fOnTp491ApTvbKQk8MDLxooVK2Ty5MnK9Rhvv\/AIw\/9wPjnjjDOs7hNY47juuuvUgKs9PeGW27t3b8UFC\/G2PxOY9YDQYN0HXoBI7777rlxzzTUCgT5w4IC1DNwEyOu5h+t1nPqENQKEOXBYOeiI2gVXH9+Ak1VxXINtKd\/ZSEngAc8vTLfB2wsJb7+YeoH34\/nnny9YmLe5T2DfD16sMNBql2O8lFxyySWyevVqOfPMM62uP9ocU6\/w+sN6jzOh7eENe+GFF1rLwE2AvJ57HF0Tp2fCGgHCxtNhw4YJFiQnTpyo3v6mTZumBifsBbIteZ2NdMoppySKh\/YAu+yyy9JrQEnoE5huwwbMe++9V21CxVTUY489pkQJnp+2PxM7duxQjhhwuca6F\/ZB6fXQl156Sbp06WItAwgQpp5fffXV9PDm1ee9Pi\/1OGmNAAHcwYMHZdKkSWoREgn7QrBI26ZNm1Jzjfx6fs5GShIPAIcoX3XVVcoVH9NvSegTcLefOnWq8gBEgkWEFy+8+Seh\/qjjunXr5MEHH1SbUOGKj3PBsC6ILRg2M3jttdfU+p\/TDdtPfeM0LlglQHrUx2I0LAA9FRe5GsT8AuQhykHB5j6B9S+83bZunXk8SJKeCaxvwALSLx\/Zj6XtfSBofePAw0oBirke8PZIgARIgARsioTA1iQBEiABEjCLAC0gs9qLd0sCJEAC1hCgAFnTlKwICZAACZhFgAJkVnvxbkmABEjAGgIUIGuakhUhARIgAbMIUIDMai\/eLQmQAAlYQ4ACZE1TsiIkQAIkYBYBCpBZ7cW7JQESIAFrCFCArGlKVoQESIAEzCJAATKrvXi3JEACJGANAQqQNU3JipAACZCAWQQoQGa1F++WBEiABKwh8P8BokybSlWjudQAAAAASUVORK5CYII=","height":337,"width":560}}
%---
