function [json_struct] = read_json_file(json_fn)
%READ_JSON_FILE
% by Jules Brochard
fid = fopen(json_fn);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);

json_struct=jsondecode(str);

end

