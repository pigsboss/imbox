function varargout = transmap(varargin)
%TRANSMAP Transform input map with specific coordinates to different coordinate
%system.
%Syntax
% transmap(mapin,cvin,in2out,out2in)
% transmap(mapin,cvin,cvout,in2out,out2in)
% transmap(mapin,cvin,cvout,in2out,out2in,INTERP_SCALE)
% transmap(mapin,cvin,cvout,in2out,out2in,INTERP_SCALE,INTERP_METHOD)
% transmap(mapin,cvin,cvout,in2out,out2in,INTERP_SCALE,INTERP_METHOD,NANVAL)
% {mapout,cvout} = transmap(...)
% [mapout,cvout] = transmap(...)
% [mapout,cvout1,cvout2] = transmap(...)
%
%INPUTS
% mapin    is input map (2-D array).
% cvin     are coordinate variables of mapin (a cell contains multiple 2-D arrays).
% cvout    are coordinate variables of output map (a cell contains multiple 2-D
% arrays).
% in2out   is handle of function convert cvin to cvout.
% out2in   is handle of function convert cvout to cvin.
%
%RETURNS
% mapout   is output map (2-D array).
% cvout    are coordinate variables of output map (a cell contains multiple 2-D
% arrays).
[mapin,cvin,cvout,in2out,~,onmapin,INTERP_METHOD,NANVAL] =...
  parse_inputs(varargin{:});
cvin_out = in2out(cvin);
map = TriScatteredInterp(cvin_out{1}(:),cvin_out{2}(:),mapin(:),INTERP_METHOD);
mapout = NANVAL;
mapout(onmapin) = map(cvout{1}(onmapin),cvout{2}(onmapin));
mapout(isnan(mapout)) = NANVAL(isnan(mapout));
nargoutchk(0,3)
switch nargout
  case 0
    varargout{1} = {mapout,cvout};
  case 1
    varargout{1} = {mapout,cvout};
  case 2
    varargout{1} = mapout;
    varargout{2} = cvout;
  case 3
    varargout{1} = mapout;
    varargout{2} = cvout{1};
    varargout{3} = cvout{2};
end
return

function [mapin,cvin,cvout,in2out,out2in,onmapin,INTERP_METHOD,NANVAL] =...
  parse_inputs(varargin)
cvout = [];
INTERP_SCALE = [];INTERP_SCALE_d = 1;
INTERP_METHOD = [];INTERP_METHOD_d = 'linear';
NANVAL = [];NANVAL_d = 0;
narginchk(4,8)
switch nargin
  case 4
    mapin = varargin{1};
    cvin = varargin{2};
    in2out = varargin{3};
    out2in = varargin{4};
  case 5
    mapin = varargin{1};
    cvin = varargin{2};
    cvout = varargin{3};
    in2out = varargin{4};
    out2in = varargin{5};
  case 6
    mapin = varargin{1};
    cvin = varargin{2};
    cvout = varargin{3};
    in2out = varargin{4};
    out2in = varargin{5};
    INTERP_SCALE = varargin{6};
  case 7
    mapin = varargin{1};
    cvin = varargin{2};
    cvout = varargin{3};
    in2out = varargin{4};
    out2in = varargin{5};
    INTERP_SCALE = varargin{6};
    INTERP_METHOD = varargin{7};
  case 8
    mapin = varargin{1};
    cvin = varargin{2};
    cvout = varargin{3};
    in2out = varargin{4};
    out2in = varargin{5};
    INTERP_SCALE = varargin{6};
    INTERP_METHOD = varargin{7};
    NANVAL = varargin{8};
end
if isempty(INTERP_SCALE),  INTERP_SCALE = INTERP_SCALE_d;end
if isempty(INTERP_METHOD), INTERP_METHOD = INTERP_METHOD_d;end
if isempty(NANVAL),        NANVAL = NANVAL_d;end
if (length(cvin)~=2) || (length(size(mapin))~=2)
  error('Support 2-D map only.')
end
if isvector(cvin{1})
  cvin1v = cvin{1};
  cvin2v = cvin{2};
  [cvin{1},cvin{2}] = meshgrid(cvin1v,cvin2v);
  clear cvin1v cvin2v
end
if isempty(cvout)
  cvout_in = in2out(cvin);
  cvout1_max = max(cvout_in{1}(:));
  cvout2_max = max(cvout_in{2}(:));
  cvout1_min = min(cvout_in{1}(:));
  cvout2_min = min(cvout_in{2}(:));
  [NUMPIX2,NUMPIX1] = size(mapin);
  NUMPIX1 = NUMPIX1*INTERP_SCALE;
  NUMPIX2 = NUMPIX2*INTERP_SCALE;
  [cvout1,cvout2] = meshgrid(1:NUMPIX1,1:NUMPIX2);
  cvout1 = (cvout1-1)/(NUMPIX1-1)*(cvout1_max-cvout1_min)+cvout1_min;
  cvout2 = (cvout2-1)/(NUMPIX2-1)*(cvout2_max-cvout2_min)+cvout2_min;
  cvout = {cvout1,cvout2};
end
cvout_in = out2in(cvout);
cvin1_max = max(cvin{1}(:));
cvin2_max = max(cvin{2}(:));
cvin1_min = min(cvin{1}(:));
cvin2_min = min(cvin{2}(:));
onmapin = (cvout_in{1}<=cvin1_max) & (cvout_in{2}<=cvin2_max) &...
  (cvout_in{1}>=cvin1_min) & (cvout_in{2}>=cvin2_min);
clear cvin1_max cvin2_max cvin1_min cvin2_min
if isscalar(NANVAL),NANVAL = NANVAL*ones(size(cvout{1}));end
return
