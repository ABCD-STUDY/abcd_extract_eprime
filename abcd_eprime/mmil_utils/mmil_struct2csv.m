function mmil_struct2csv(info,fname,separator)
%function mmil_struct2csv(info,fname,[separator])
%
% Required Input:
%  info: struct array
%  fname: output file name
%
% Optional Input:
%   separator: separator between columns
%     {default: ','}
%
% Created:  10/17/14 by Don Hagler
% Prev Mod: 01/30/15 by Don Hagler
% Last Mod: 05/12/20 by Don Hagler
%   

if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('separator','var') || isempty(separator), separator = ','; end;

tags = fieldnames(info);
data = squeeze(struct2cell(info));
data = cat(2,tags,data)';
mmil_write_csv(fname,data,'separator',separator);

