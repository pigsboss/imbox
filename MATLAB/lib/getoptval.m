function varargout = getoptval(optpairs,optnames,optdefaults)
%GETOPTVAL Get values of options from inputs.
%
%Syntax
% optvals = getoptval(optpairs, optnames)
% optvals = getoptval(optpairs, optnames, optdefaults)
%
%Description
%optpairs is user inputs. Like in many command-line tools, user options are
%provided in pairs of option names and values, separated by commas.
%
%optnames is a cell of strings. Each string identifies the name of a option.
%Names are case-insensitive.
%
%optvals is a cell of variabls.
%

optvals = cell(size(optnames));
if nargin == 2
  optdefaults = cell(size(optnames));
end
for k = 1:numel(optnames)
  match = strcmpi(optnames{k}, optpairs);
  kid = find(match, 1, 'last');
  if ~isempty(kid)
    optvals{k} = optpairs{kid+1};
  else
    optvals{k} = optdefaults{k};
  end
end
varargout = optvals;
return