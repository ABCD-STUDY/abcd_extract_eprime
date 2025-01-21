function fname_out = abcd_check_eprime_sprdsh(fname,cols,fields,outdir,outstem,forceflag,verbose)
%function fname_out = abcd_check_eprime_sprdsh(fname,[cols],[fields],[outdir],[outstem],[forceflag],[verbose])
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
%   outstem: output file stem
%     if empty, will use filestem of fname
%     {default = []}
%   forceflag: [0|1] overwrite existing output
%     {default = 0}
%   'verbose': [0|1] display messages
%     {default = 1}
%
% Created:  11/02/17 by Dani Cornejo
% Prev Mod: 01/23/19 by Dani Cornejo
% Prev Mod: 01/14/20 by Octavio Ruiz
% Prev Mod: 07/02/20 by Don Hagler
% Prev Mod: 08/28/20 by Don Hagler
% Prev Mod: 11/01/20 by Don Hagler
% Prev Mod: 05/24/24 by Don Hagler
% Last Mod: 08/28/24 by Don Hagler
%

% Based on abcd_check_eprime_encoding.m
% Created:  12/21/16 by Don Hagler
% Last Mod: 01/06/17 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('cols','var') || isempty(cols), cols = []; end;
if ~exist('fields','var') || isempty(fields), fields = []; end;
if ~exist('outdir','var') || isempty(outdir), outdir = pwd; end;
if ~exist('outstem','var'), outstem = []; end;
if ~exist('forceflag','var') || isempty(forceflag), forceflag = 0; end;
if ~exist('verbose','var') || isempty(verbose), verbose = 1; end;

if ~exist(fname,'file')
  error('file %s not found',fname);
end;

python_flag = check_python;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(fname,'file')
  error('file %s not found',fname);
end;

if isempty(outstem)
  [~,outstem] = fileparts(fname);
end
outstem = abcd_clean_fstem(outstem);

% write event info to file
fname_out = sprintf('%s/%s_events.csv',outdir,outstem);

if ~exist(fname_out,'file') || forceflag
  if python_flag
    % find, read, intrepret eprime file, and write contents to specified's name file
    % if the operation is succesful, this file will be read below
    if verbose, fprintf('%s: reading e-prime file using Python...\n',mfilename); end
    diagnos = -256;
    fname_ck = sprintf('%s/%s_checked.txt',outdir,outstem);
    cmd = sprintf('python3 $MMPS_DIR/python/eprime_sprdsht_get.py %s ExportFile %s',...
      clean_fname(fname),clean_fname(fname_ck));

    [status, cmdout] = mmil_unix(cmd);

    if status == 0
      % check if input file was found and interpreted; and get diagnostic variable,
      %   located at the beginning of the returned string (Octavio, 2020jan14)
      try
        [var, moreinfo] = strtok(cmdout, ',');
        diagnos = str2double( strtrim(var) );
      catch
        python_flag = 0;
      end
    else
      if verbose, fprintf('%s: WARNING: python cmd %s failed:\n%s\n',mfilename,cmd,cmdout); end
      python_flag = 0;
    end
    
    if python_flag  &&  diagnos > 0  &&  exist(fname_ck,'file')
      ; % The eprime file was succesfully read, translated, and a copy was rewritten
    else
      if verbose, fprintf('%s: WARNING: eprime_sprdsht_get.py: unable to find, read, intrepret, or write file\n', mfilename); end
      python_flag = 0;
    end
  end

  if ~python_flag
    if verbose, fprintf('%s: reading e-prime file...\n',mfilename); end
    % check for UTF-16, convert to ASCII if necessary
    % read input file, allowing either comma or tab delimited
    fname_ck = abcd_check_eprime_encoding(fname,outdir,outstem,forceflag);
  end;

  if ~isempty(cols) || ~isempty(fields)
    [tmp,fstem,fext] = fileparts(fname_ck);
    if strcmp(fext,'.txt')
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
    if isempty(B)
      if verbose, fprintf('%s: WARNING: failed to read events\n',mfilename); end
      fname_out = [];
      return;
    end;
    mmil_write_csv(fname_out,B);
  else
    if verbose, fprintf('%s: WARNING: column names and/or fields empty: no event file \n',mfilename); end
    fname_out = []; 
  end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function python_flag = check_python()
  python_flag = 0;
  cmd = 'which python3';
  [status,result,errmsg] = jsystem(cmd);
  if ~status && isempty(regexp(result,'not found')) && isempty(regexp(errmsg,'not found'))
    python_flag = 1;
  end;
return; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add escape characters to file name for shell commands
function fname_clean = clean_fname(fname)
  fname_clean = fname;
  % add escape character before parentheses
  fname_clean = regexprep(fname_clean,'(','\\(');
  fname_clean = regexprep(fname_clean,')','\\)');
  % add escape character before spaces
  fname_clean = regexprep(fname_clean,' ','\\ ');
  % add escape character before &
  fname_clean = regexprep(fname_clean,'&','\\&');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
