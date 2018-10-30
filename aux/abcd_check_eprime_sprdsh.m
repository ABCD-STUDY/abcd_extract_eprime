function fname_csv = abcd_check_eprime_sprdsh(fname,cols,fields,outdir,forceflag)
%function fname = abcd_check_eprime_sprdsh(fname,[cols],[fields],[outdir],[forceflag])
%
% Purpose: Reads an ASCII file containing an E-Prime spreadsheet. 
%          Extract contents while fixing encoding and format issues. 
%          Check whether eprime file is plain text or unicode UTF-16
%          and if so, convert to ASCII. Outputs _checked.txt file.  
%          Also generates _events.csv file with the information relevant 
%          to the analysis, only.   
%
% Required parameters:
%   fname: full path name of eprime output file
%
% Optional parameters:
%   cols: Columns to keep from the e-prime file for the analysis. 
%         If empty, it does not generate _event.csv files. 
%     {default = []}
%   fields: new names for the columns in _event.csv 
%     {default = []}
%   outdir: output directory
%     output for _checked.csv and _events.csv 
%     {default = pwd}
%   forceflag: [0|1] overwrite existing output
%     {default = 0}
%
% Based on abcd_check_eprime_encoding.m
% Created:  12/21/16 by Don Hagler
% Last Mod: 01/06/17 by Don Hagler
%
% Created:  11/02/17 by Dani Cornejo
% Prev Mod: 06/08/18 by Dani Cornejo
% Last Mod: 07/25/18 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('cols','var') || isempty(cols), cols = []; end;
if ~exist('fields','var') || isempty(fields), fields = []; end;
if ~exist('outdir','var') || isempty(outdir), outdir = pwd; end;
if ~exist('forceflag','var') || isempty(forceflag), forceflag = 0; end;

if ~exist(fname,'file')
  error('file %s not found',fname);
end;

python_flag = check_python;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(fname,'file')
  error('file %s not found',fname);
end;

fname_clean = fname;
fname_clean = regexprep(fname_clean,'(','\\(');
fname_clean = regexprep(fname_clean,')','\\)');
fname_clean = regexprep(fname_clean,' ','\\ ');

[~,fstem,~] = fileparts(fname);
fstem_clean = regexprep(fstem,'\s.+','');

% write event info to file
fname_csv = sprintf('%s/%s_events.csv',outdir,fstem_clean);  
if ~exist(fname_csv,'file') || forceflag

  if python_flag
    fprintf('%s: reading e-prime file with python...\n',mfilename);
    fname_ck = sprintf('%s/%s_checked.txt',outdir,fstem_clean); 
    cmd = sprintf('python3 $MMPS_DIR/python/eprime_sprdsht_get.py %s ExportFile %s',...
      fname_clean,fname_ck);
    [status, cmdout] = unix(cmd);
    if status <= 0 || ~exist(fname_ck,'file')
      fprintf('%s: WARNING: python cmd %s failed:\n%s\n',mfilename,cmd,cmdout);
      python_flag = 0;
    end;
  end;
  
  if ~python_flag
    fprintf('%s: reading e-prime file...\n',mfilename);
    % check for UTF-16, convert to ASCII if necessary
    % read input file, allowing either comma or tab delimited
    fname_ck = abcd_check_eprime_encoding(fname,outdir,forceflag);
  end;

  if ~isempty(cols) || ~isempty(fields)
    [tmp,fstem,fext] = fileparts(fname_ck);
    if strcmp(fext,'.csv')
      A = mmil_readtext(fname_ck,'[,]');
    elseif strcmp(fext,'.txt')
      A = mmil_readtext(fname_ck,'[\t]');
    else
      A = mmil_readtext(fname_ck,'[,\t]');
    end;
    
    % remove usual first row unless it is missing
    if ~strcmp(char(A(1,1)),'ExperimentName')
      A = A(2:end,:);
    end;
    % select columns of interest
    colnames = A(1,:); 
    colnames = colnames(~cellfun('isempty',colnames));
    vals = A(2:end,:);
    [~,i_all,i_sel] = intersect(colnames,cols);
    [i_sel,i_sort] = sort(i_sel);
    i_all = i_all(i_sort);
    % create new matrix with replacement column labels
    vals = vals(:,i_all);
    tags = fields(i_sel);
    B = cat(1,tags,vals);
    % write to csv
    mmil_write_csv(fname_csv,B);
  else
    fprintf('%s: column names and/or fields empty: no event file \n',mfilename);
    fname_csv = []; 
  end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function python_flag = check_python()
  python_flag = 0;
  cmd = 'which python3';
  [s,r] = unix(cmd);
  if ~s && isempty(regexp(r,'not found'))
    python_flag = 1;
  end;
return; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

