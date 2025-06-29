function DCM = ijk2SEZ(l, lat)
% ----------------------------------------------------------------------
%  Rotational matrix from IJK frame to SEZ frame
%    20170310  y.yoshimura
%    Inputs: l: local sidereal time, rad
%            lat: local latitude, rad            
%    outputs: DCM: %  Rotational matrix from IJK frame to SEZ frame
%   related function files:
%   cf:
%   revisions;
%   function DCM = ijk2SEZ(l, lat)
%   (c) 2021 yasuhiro yoshimura
%----------------------------------------------------------------------

DCM = zyx2DCM(l, pi/2-lat, 0);

end