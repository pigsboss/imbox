function imprint(varargin)
[picname,format,xlabelstr,ylabelstr,zlabelstr,titlestr,clabelstr,labelsize,...
  titlesize,axessize,cmap,fontname,cbar] = parse_inputs(varargin{:});
if ~isempty(cmap)
  colormap(cmap)
end
if ~isempty(clabelstr)
  if isempty(cbar),cbar = colorbar;end
  set(get(cbar,'YLabel'),'string',clabelstr,'fontsize',labelsize,'fontname',fontname);
end
if ~isempty(xlabelstr),xlabel(xlabelstr,'fontsize',labelsize,'fontname',fontname);end
if ~isempty(ylabelstr),ylabel(ylabelstr,'fontsize',labelsize,'fontname',fontname);end
if ~isempty(zlabelstr),zlabel(zlabelstr,'fontsize',labelsize,'fontname',fontname);end
if ~isempty(titlestr), title(titlestr,'fontsize',titlesize,'fontname',fontname);end
set(gca,'fontsize',axessize,'fontname',fontname)
axis xy
axis image
if ~isempty(picname)
  if strcmpi(format,'bitmap')
    print('-dpng',[picname,'.png'])
  else
    if exist('epswrite.m','file')~=2
      print('-depsc2',[picname,'.eps'])
    else
      epswrite([picname,'.eps'])
    end
  end
end
return

function [picname,format,xlabelstr,ylabelstr,zlabelstr,titlestr,clabelstr,labelsize,...
  titlesize,axessize,cmap,fontname,cbar] = parse_inputs(varargin)
picname=[];
format='';format_d='vector';
xlabelstr='';
ylabelstr='';
zlabelstr='';
titlestr=[];
clabelstr='';
labelsize=[];
titlesize=[];
axessize=[];
cmap='';
fontname='';fontname_d='Helvetica';
fontsize=[];
k = 1;
while k < nargin
  if isscalar(varargin{k}) || ischar(varargin{k})
    switch upper(varargin{k})
      case {'PICNAME'}
        k = k+1;
        picname = varargin{k};
      case {'FORMAT'}
        k = k+1;
        format = varargin{k};
      case {'XLABEL','X'}
        k = k+1;
        xlabelstr = varargin{k};
      case {'YLABEL','Y'}
        k = k+1;
        ylabelstr = varargin{k};
      case {'ZLABEL','Z'}
        k = k+1;
        zlabelstr = varargin{k};
      case {'TITLE'}
        k = k+1;
        titlestr = varargin{k};
      case {'CLABEL','C'}
        k = k+1;
        clabelstr = varargin{k};
      case {'LABELSIZE'}
        k = k+1;
        labelsize = varargin{k};
      case {'TITLESIZE'}
        k = k+1;
        titlesize = varargin{k};
      case {'AXESSIZE'}
        k = k+1;
        axessize = varargin{k};
      case {'CMAP','COLORMAP'}
        k = k+1;
        cmap = varargin{k};
      case {'FONTNAME'}
        k = k+1;
        fontname = varargin{k};
      case {'FONTSIZE'}
        k = k+1;
        fontsize = varargin{k};
    end
    k = k+1;
  else
    error('Syntax error.')
  end
end
cbar = findobj(gcf,'Tag','Colorbar');
if isempty(format),   format=format_d;end
if isempty(fontname), fontname=fontname_d;end
if isempty(labelsize)
  if isempty(fontsize)
    labelsize=get(get(gca,'XLabel'),'fontsize');
  else
    labelsize=fontsize;
  end
end
if isempty(axessize)
  if isempty(fontsize)
    axessize=get(gca,'fontsize');
  else
    axessize=fontsize;
  end
end
if isempty(titlesize)
  if isempty(fontsize)
    titlesize=get(get(gca,'Title'),'fontsize');
  else
    titlesize=fontsize;
  end
end
if isempty(xlabelstr),xlabelstr=get(get(gca,'XLabel'),'string');end
if isempty(ylabelstr),ylabelstr=get(get(gca,'YLabel'),'string');end
if isempty(zlabelstr),zlabelstr=get(get(gca,'ZLabel'),'string');end
if isempty(titlestr), titlestr=get(get(gca,'Title'),'string');end
if isempty(clabelstr) && ~isempty(cbar)
  clabelstr=get(get(cbar,'YLabel'),'string');
end
return
