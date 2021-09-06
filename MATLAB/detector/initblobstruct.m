function [bmes, bimg] = initblobstruct(NB)
%initialize blob parameters and images
bmes = struct(...
  'xctr',zeros(NB,1),...              % x index of central pixel
  'yctr',zeros(NB,1),...              % y index of central pixel
  'scale',zeros(NB,1),...             % linear-scale of blob (sqrt)
  'rect',zeros(NB,4),...              % rectangle of blob [X Y W H]
  'xctrsubpix',zeros(NB,1),...        % sub-pixel x of blob barycenter
  'yctrsubpix',zeros(NB,1),...        % sub-pixel y of blob barycenter
  'flux',zeros(NB,1),...              % integral flux of blob
  'peak',zeros(NB,1),...              % intensity of brightest pixel of blob
  'eccentricity',zeros(NB,1),...      % eccentricity of border ellipse
  'semiminor',zeros(NB,1),...         % semi-minor axis of border ellipse
  'semimajor',zeros(NB,1),...         % semi-major axis of border ellipse
  'posangle',zeros(NB,1),...          % position angle of major axis
  'elliparam',zeros(NB,5));           % elliptical parameters of blob [xc yc rx ry phi]
bimg = struct(...
  'image',cell(NB,1),...              % image of blob (within the border)
  'mask',cell(NB,1),...               % (boolean) mask of blob (within the border)
  'xbdply',cell(NB,1),...             % x of border polygon vertices
  'ybdply',cell(NB,1),...             % y of border polygon vertices
  'fgrdimg',cell(NB,1),...            % foreground image of blob
  'bgrdimg',cell(NB,1));              % inpainted background image of blob
return
