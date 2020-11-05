function fstem_clean = abcd_clean_fstem(fstem)
%function fstem_clean = abcd_clean_fstem(fstem)
%
% Purpose special characters from output file stem
%
% Created:  11/01/2020 DH
% Last Mod: 11/01/2020 DH
%

if ~mmil_check_nargs(nargin,1), return; end;

fstem_clean = fstem;

% remove spaces
fstem_clean = regexprep(fstem_clean,'\s.+','');
% remove special characters
fstem_clean = regexprep(fstem_clean,'[!@#$%^&*()\\?/]','');
% remove quotes
fstem_clean = regexprep(fstem_clean,'''','');
% remove extra extension
fstem_clean = regexprep(fstem_clean,'\..+','');

return;
