function mmil_struct2csv(info,fname)
%function mmil_struct2csv(info,fname)
%
% Required Input:
%  info: struct array
%  fname: output file name
%
% Created:  10/17/14 by Don Hagler
% Last Mod: 01/30/15 by Don Hagler
%   

if ~mmil_check_nargs(nargin,2), return; end;

tags = fieldnames(info);
data = squeeze(struct2cell(info));
data = cat(2,tags,data)';
mmil_write_csv(fname,data);

