% ----------------------------------------------------------------------
%   data handling
%    20190706  y.yoshimura
%    Inputs:
%   Outputs:
%   related function files:
%   note:
%   cf:
%   revisions;
%
%   (c) 2019 yasuhiro yoshimura
%----------------------------------------------------------------------

clc
clear
close all

format long

fileName = {'CowellForm.txt'
    'outputCowell.txt'}

for i = 1:length(fileName)
    output = importdata(fileName{i});
    
    % Time(day), a(km), e, i(deg), n(deg), w(deg), m(deg), hp(km)
    time(:,:,i) = output.data(:,1);
    a(:,:,i) = output.data(:,2);
    e(:,:,i) = output.data(:,3);
    inc(:,:,i) = output.data(:,4);
    n(:,:,i) = output.data(:,5);
    w(:,:,i) = output.data(:,6);
    m(:,:,i) = output.data(:,7);
    hp(:,:,i) = output.data(:,8);
    
end

figure
plot(a(:,1,1)-a(:,1,2), 'r')
hold on
% plot(a(:,1,2), 'g')

[a(:,1,1) a(:,1,2)]