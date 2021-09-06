function data = parse_config_file(filename)
%PARSE_CONFIG_FILE Parse config file and return a matlab struct contains
%all the keys and values defined in the config file.
%
%Definition of config file
%key = value
%key is case insensitive.
%comment line starts with # or %

fid = fopen(filename);
data = struct();
tline = fgetl(fid);
while ischar(tline)
  if numel(tline)>1
    if tline(1)~='#' && tline(1)~='%'
      dlmt = find(tline=='=',1,'first');
      if ~isempty(dlmt)
        key = lower(strrep(tline(1:(dlmt-1)),' ',''));
        val = strrep(tline((dlmt+1):end),' ','');
        try
          data.(key) = val;
        catch err
          if isempty(strfind(err.identifier,'InvalidFieldName'))
            rethrow(err)
          end
        end
      end
    end
  end
  tline = fgetl(fid);
end
fclose(fid);
return
