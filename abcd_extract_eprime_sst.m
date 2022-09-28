function [eprime_nruns,errcode,behav,errmsg] = abcd_extract_eprime_sst(fname,varargin)
%function [eprime_nruns,errcode,behav,errmsg] = abcd_extract_eprime_sst(fname,[options])
%
% Purpose: extract condition time courses
%   from eprime data files for SST task
%
% Required input:
%   fname: name of input eprime file
%
% Optional input:
%   'outdir': output directory
%     {default = pwd}
%   'outstem': output file stem
%     if empty, will use filestem of fname
%     {default = []}
%   'TR': scan repetition time (s)
%     {default = 0.8}
%   'numTRs': number of repetitions
%     {default = 500}
%   'minfrac': minimum fraction of a TR to register an event
%     {default = 0.5}
%   'nskipTRs': number of TRs to remove from beginning of time course
%     {default = 0}
%   'switch_thresh': percent correct below which to swap button assignements
%     {default = 0.3}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%   'timing_files_flag': [0|1] generate timing files
%     {default = 1}
%   'verbose': [0|1] display messages
%     {default = 1}
%
% Output: 
%    eprime_nruns: number of valid runs in the e-prime in the file  
%    errcode: [0|1] whether the file was successfully processed
%    behav: behavioral data
%    errmsg: string describing error if errcode=1
%
%    Since January 2022, behav includes "glitch" variables
%    defined by Hugh Caravan et al

% Created:  10/07/16 by Don Hagler
% Prev Mod: 02/12/19 by Dani Cornejo
% Prev Mod: 03/12/19 by Feng Xue
% Prev Mod: 08/05/19 by Octavio Ruiz
% Prev Mod: 05/03/21 by Don Hagler
% Prev Mod: 08/31/21 by Don Hagler
% Prev Mod: 01/12/22 by Octavio Ruiz
% Prev Mod: 01/14/22 by Octavio Ruiz
% Prev Mod: 02/10/22 by Octavio Ruiz
% Prev Mod: 04/22/22 by Don Hagler
% Prev Mod: 05/02/22 by Don Hagler
% Last Mod: 06/27/22 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: Based on make-SST-refs.py
%         provided by Michael Michael Riedel, Ph.D. from FIU
%       Using style from abcd_extract_eprime_mid.m
%         provided by Donald J Hagler, Ph.D.
%       Includes code for glitch and violator variables by Octavio Ruiz
%         adapted from python code provided by Sage Hahn from UVM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize outputs
eprime_nruns = []; errcode = 0; behav = []; errmsg = [];

% check arguments 
if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
parms = check_input(fname,varargin); 

% create output directory
mmil_mkdir(parms.outdir);

% create structure with info for each event
% write field _events from original; no processing
[event_info,start_time,all_types,errcode,errmsg] = get_event_info(parms);  
if errcode, return; end;

% switch buttons if necessary 
% file _switched.csv will be created 
[event_info,parms.switch_flag,errcode,errmsg] = sst_switch(event_info,parms); 
if errcode, return; end;

% Replace and/or add events names that were not coded 
% correctly by the e-prime software. For example trials 
% identified as late or in the anticipation window. 
% It inputs event_info from previus two steps 
all_types_new = get_sst_types(event_info,parms);

% Write the new events found in previus step 
% to a csv file. No change in the timing. 
write_event_revised(parms,all_types_new);

% Get behav data and write it to a csv file 
% using revised events saved in all_types_new
% timing cells don't change
[behav,eprime_runs,errcode,errmsg] = get_behavioral_data(event_info,parms,all_types_new);
if errcode, return; end;

parms.eprime_runs = eprime_runs; 
eprime_nruns = length(eprime_runs);
parms.eprime_nruns = eprime_nruns; 
if parms.verbose, fprintf('%s: number of e-prime runs = %d \n',mfilename,eprime_nruns); end

if parms.timing_files_flag && eprime_nruns
  % create files for each condition
  for i=1:parms.nconds
    type  = parms.typenames{i};
    cond  = parms.condnames{i}; 
    pat = '\w_(?<Event>\w+)'; 
    get_event = regexp(cond,pat,'names');
    event = get_event.Event;

    % use all_types_new, which is the revised information, to find all
    % types of events 
    ind_type = find(strcmp(type,all_types_new));  
    switch event  
      case 'go'
        % get times of stim onset and offset
        parms.mindur = 0.08;  
        onset  = {event_info(ind_type).(['go_onset_time'])};
        offset = {event_info(ind_type).(['go_offset_time'])};    
        [onset,offset] = check_offsets(onset,offset,event,parms);
      case 'stop'
        % get times of stim onset and offset
        parms.mindur = 0;  
        onset  = [event_info(ind_type).(['stop_start_time'])];
        offset = onset+20;         
    end
    
    % find most recent start time for each event
    [ind_start,event_ref_time] = set_ref(onset,start_time);  
    % create files for each scan
    write_files(cond,onset,offset,ind_start,event_ref_time,parms); 
  end 
end; %timing_files_flag

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname,options)
  parms = mmil_args2parms(options,{...
    'fname',fname,[],...
    ...
    'outdir',pwd,[],...
    'outstem',[],[],...
    'TR',0.8,[],...
    'numTRs',500,[],...
    'minfrac',0.5,[0,1],...
    'nskipTRs',0,[],...
    'switch_thresh',0.3,[0,0.5],...
    'forceflag',false,[false true],...
    'timing_files_flag',true,[false true],...
    'verbose',true,[false true],...
    'eprime_runs',1:2,[],...   
    'eprime_nruns',0,0:2,...
    ...
    'colnames',  {'NARGUID','SessionDate','SessionTime','ExperimentName','ExperimentVersion',...
                  'Subject','Procedure[Block]','Block','Trial','GetReady.RTTime', 'BeginFix.StartTime',...
                  'Go.RT', 'Fix.RT', 'StopSignal.RT','Go.Duration', 'SSDDur',...
                  'StopSignal.Duration','SSD.RT', 'Go.OnsetTime', 'Go.OffsetTime',  'Go.ACC' , 'Go.RESP',...
                  'Fix.RESP', 'SSD.OnsetTime','StopSignal.StartTime','SSD.OffsetTime','StopSignal.ACC', 'Go.CRESP',...
                  'StopSignal.RESP', 'SSD.ACC','TrialCode', 'SSD.RESP' },[],...
    'fieldnames',{'narguid','date','time','experiment','version',...
                  'subject','procedure_block','block','trial','getready_rttime','beginfix_st',...
                  'go_rt','fix_rt','stop_rt','go_dur','ssd_dur',...
                  'stop_dur','ssd_rt','go_onset_time','go_offset_time' , 'go_acc','go_resp',...
                  'fix_resp', 'ssd_onset_time','stop_start_time','ssd_offset_time','stop_acc','go_cresp',...
                  'stop_resp', 'ssd_acc', 'type', 'ssd_resp' },[],...
    'typenames',{'CorrectGo', 'CorrectLateGo', 'IncorrectGo', 'IncorrectLateGo','NoRespGo',... 
                 'CorrectStop', 'IncorrectStop','SSDStop'},[],...
    'condnames',{'correct_go','correctlate_go', 'incorrect_go', 'incorrectlate_go' ...
                  'noresp_go', 'correct_stop', 'incorrect_stop','ssd_stop'},[],...
  });

  % check conditions
  parms.nconds = length(parms.condnames);
  parms.ntypes = length(parms.typenames);
  parms.mindur = parms.minfrac*parms.TR;
  if parms.nconds ~= parms.ntypes
    error('condnames and typenames length mismatch');
  end;
  % check input file
  if ~exist(parms.fname,'file')
    error('file %s not found',parms.fname);
  end;
  [fdir,fstem,fext] = fileparts(parms.fname);
  % remove problematic characters
  if isempty(parms.outstem)
    parms.outstem = abcd_clean_fstem(fstem);
  end;
  % compute time at the end of each TR
  parms.TR_offset = linspace(parms.TR,parms.numTRs * parms.TR,parms.numTRs);
  parms.TR_onset  = parms.TR_offset - parms.TR;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [event_info,start_time,all_types,errcode,errmsg] = get_event_info(parms)
  
  event_info=[]; start_time=[]; all_types=[];
  errcode = 0; errmsg = [];
  
  try 
    % write event info to file
    fname_csv = abcd_check_eprime_sprdsh(parms.fname, parms.colnames,...
                  parms.fieldnames, parms.outdir, parms.forceflag, parms.verbose);
    event_info = mmil_csv2struct(fname_csv);
  catch me
    if parms.verbose, fprintf('%s: ERROR: failed to read or interpret e-prime file %s:\n%s\n',...
      mfilename,parms.fname,me.message); end
    errcode = 1;
    errmsg = 'failed to read or interpret e-prime file';
    return;
  end
  
  % check experiment
  experiment = mmil_getfield(event_info(1),'experiment');
  if isempty(regexpi(experiment,'sst'))
    if parms.verbose, fprintf('%s: ERROR: wrong experiment name in e-prime file %s: %s\n',...
      mfilename,parms.fname,experiment); end
    errcode = 1;
    errmsg = 'wrong experiment name';
    return;
  else
    if parms.verbose, fprintf('%s: experiment name: %s\n',mfilename,experiment); end
  end;
  
  % get start times
  ind_start = find(~cellfun(@isempty,{event_info.beginfix_st})); 
  start_time = [event_info(ind_start).beginfix_st]; 
  % remove non-events
  all_types = {event_info.type}; 
  ind_events = find(~cellfun(@isempty,all_types)); 
  event_info = event_info(ind_events); 
  all_types = {event_info.type};

  % check go_resp
  go_resp = {event_info.go_resp};
  fix_resp = {event_info.fix_resp};
  stop_resp = {event_info.stop_resp};
  go_cresp = {event_info.go_cresp};  
  if any(cellfun(@isstr,go_resp)) || any(cellfun(@iscell,go_resp)) ||...
     any(cellfun(@isstr,fix_resp)) || any(cellfun(@iscell,fix_resp)) ||...
     any(cellfun(@isstr,stop_resp)) || any(cellfun(@iscell,stop_resp))     
    if ~any(cellfun(@isstr,go_cresp))
      if parms.verbose, fprintf('%s: ERROR: string go_resp values without string go_cresp values in e-prime file %s\n',...
        mfilename,parms.fname); end
      errcode = 1;
      errmsg = 'string go_resp values without string go_cresp values';
      return;
    end;
    % get response numbers and names for unique, non-emptpy correct responses
    uniq_go_cresp = unique(go_cresp(find(~cellfun(@isempty,go_cresp))));
    num_cresp = length(uniq_go_cresp);
    uniq_go_cresp_nums = [];
    uniq_go_cresp_names = [];
    for i=1:length(uniq_go_cresp)
      %   allow either '1{LEFTARROW}' or '1,{LEFTARROW}'
      k = regexp(uniq_go_cresp{i},'(?<num>\d+),?{(?<name>\w+)}','names');
      if isempty(k)
        if parms.verbose, fprintf('%s: ERROR: string go_cresp with unexpected pattern (%s) in e-prime file %s\n',...
          mfilename,uniq_go_cresp{i},parms.fname); end
        errcode = 1;
        errmsg = 'string go_cresp with unexpected pattern';
        return;
      end;
      uniq_go_cresp_nums(i) = str2num(k.num);
      uniq_go_cresp_names{i} = k.name;
    end;
    if parms.verbose
      if parms.verbose, fprintf('%s: WARNING: replacing go_resp strings with numeric\n', mfilename); end
    end;
    for i=1:length(event_info)
      % assign numbers to each go_resp
      go_resp = event_info(i).go_resp;
      if iscell(go_resp), go_resp = go_resp{1}; end;
      if isstr(go_resp)
        k = find(strcmp(go_resp,uniq_go_cresp_names));
        go_resp = uniq_go_cresp_nums(k);
      end;
      event_info(i).go_resp = go_resp;
      % assign numbers to each fix_resp
      fix_resp = event_info(i).fix_resp;
      if iscell(fix_resp), fix_resp = fix_resp{1}; end;
      if isstr(fix_resp)
        k = find(strcmp(fix_resp,uniq_go_cresp_names));
        fix_resp = uniq_go_cresp_nums(k);
      end;
      event_info(i).fix_resp = fix_resp;
      % assign numbers to each stop_resp
      stop_resp = event_info(i).stop_resp;
      if iscell(stop_resp), stop_resp = stop_resp{1}; end;
      if isstr(stop_resp)
        k = find(strcmp(stop_resp,uniq_go_cresp_names));
        stop_resp = uniq_go_cresp_nums(k);
      end;
      event_info(i).stop_resp = stop_resp;
      % assign numbers to each ssd_resp
      ssd_resp = event_info(i).ssd_resp;
      if iscell(ssd_resp), ssd_resp = ssd_resp{1}; end;
      if isstr(ssd_resp)
        k = find(strcmp(ssd_resp,uniq_go_cresp_names));
        ssd_resp = uniq_go_cresp_nums(k);
      end;
      event_info(i).ssd_resp = ssd_resp;
    end;
  end;
  
  % check go_cresp
  if any(cellfun(@isstr,{event_info.go_cresp}))
    if parms.verbose, fprintf('%s: WARNING: replacing go_cresp strings with numeric\n', mfilename); end
    for i=1:length(event_info)
      if isstr(event_info(i).go_cresp)
        % remove text description of response
        %   allow either '1{LEFTARROW}' or '1,{LEFTARROW}'
        event_info(i).go_cresp = str2num(regexprep(event_info(i).go_cresp,'[,{].+',''));
      end;
    end;
  end;
  
  % reset variables indicating accurate responses
  event_info = sst_acc(event_info);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function event_info = sst_acc(event_info)
  for i=1:length(event_info)
    if event_info(i).go_resp == event_info(i).go_cresp
      event_info(i).go_acc=1;
    elseif event_info(i).go_resp ~= event_info(i).go_cresp
      event_info(i).go_acc=0;
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [event_info,switch_flag,errcode,errmsg] = sst_switch(event_info,parms)
  [switch_flag,acc1,errcode,errmsg] = sst_switch_flag(event_info,parms); 
  if errcode, return; end;
 
  if switch_flag 
    event_info = sst_switch_event(event_info);
    [~,acc2,~] = sst_switch_flag(event_info,parms);
    if acc2 < acc1
      event_info = sst_switch_event(event_info);
      switch_flag = 0; 
    end 
  end
 
  % write to csv
  if switch_flag 
    fname_csv_out  = sprintf('%s/%s_events_switched.csv',...
        parms.outdir,parms.outstem); 
    if ~exist(fname_csv_out,'file') || parms.forceflag
      mmil_struct2csv(event_info,fname_csv_out)
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function event_info = sst_switch_event(event_info)
  for i=1:length(event_info)
    event_info(i).go_resp = switch_val(event_info(i).go_resp);
    event_info(i).fix_resp = switch_val(event_info(i).fix_resp);
    event_info(i).stop_resp = switch_val(event_info(i).stop_resp);
    event_info(i).ssd_resp = switch_val(event_info(i).ssd_resp);
  end
  event_info = sst_acc(event_info);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function val = switch_val(val)
  if val==1
    val=2;
  elseif val==2
    val=1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [switch_flag,accuracy,errcode,errmsg] = sst_switch_flag(event_info,parms) 
  switch_flag = 0;
  errcode = 0; errmsg = [];
  
  all_types = [event_info.go_acc];
  correct_total = size(all_types,2);
  correct = sum(all_types); 
  accuracy = 100*correct/correct_total;
  if parms.verbose, fprintf('%s: accuracy = %0.1f%%\n',mfilename,accuracy); end
  
  if accuracy == 0
    if parms.verbose, fprintf('%s: ERROR: accuracy equals zero for %s\n',mfilename,parms.fname); end
    errcode = 1;
    errmsg = 'accuracy equals zero';
    return; 
  elseif accuracy < 100*parms.switch_thresh
    if parms.verbose, fprintf('%s: accuracy < %0.1f%%, switching button responses\n',...
      mfilename,100*parms.switch_thresh); end
    switch_flag = 1;
  end;  
return 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function write_event_revised(parms,all_types_new)

  % purpose: Replace and/or add events names that were not coded 
  %          correctly by the e-prime software and write to a csv file. 
  %          all_types_new already has this informantion computed by 
  %          get_sst_types(). 

  % output file name 
  if parms.switch_flag
    fname_csv_in  = sprintf('%s/%s_events_switched.csv',parms.outdir,parms.outstem);
  else
    fname_csv_in  = sprintf('%s/%s_events.csv',parms.outdir,parms.outstem);
  end;
  fname_csv_out = sprintf('%s/%s_events_revised.csv',parms.outdir,parms.outstem); 
  % output file writing 
  if ~exist(fname_csv_out,'file') || parms.forceflag
    % event info 
    event_info = mmil_csv2struct(fname_csv_in);    
    % find events
    all_types = {event_info.type}; 
    ind_events = find(~cellfun(@isempty,all_types));  
    for i=1:length(ind_events)
      event_info(ind_events(i)).type = char(all_types_new(i));     
    end 
    % write to csv
    if ~exist(fname_csv_out,'file') || parms.forceflag
      mmil_struct2csv(event_info,fname_csv_out);
    end; 
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function all_types_new = get_sst_types(event_info,parms)
  all_types_new = {event_info.type};  
  % correct_go
  go_cresp = {event_info.go_cresp};
  go_resp = {event_info.go_resp};  
  correct_go_ind = [];
  for i=1:size(go_cresp,2)
    if go_resp{i} == go_cresp{i}
      correct_go_ind = [correct_go_ind i]; 
      all_types_new{i} = 'CorrectGo';
    end 
  end
  % correct_go_late 
  fix_resp = {event_info.fix_resp}; 
  correct_go_late_ind = [];
  for i=1:size(go_cresp,2)  
    if fix_resp{i} == go_cresp{i} 
      if go_resp{i} == go_cresp{i}
      else 
        correct_go_late_ind = [correct_go_late_ind i]; 
        all_types_new{i} = 'CorrectLateGo';
     end 
    end 
  end
  % no_resp_go 
  go_event_ind = find(~cellfun(@isempty,go_cresp)); 
  noresp_ind = find(cellfun(@isempty,go_resp)); 
  nofix_ind = find(cellfun(@isempty,fix_resp));
  noresp_nofix_ind = [intersect(noresp_ind,nofix_ind)];  
  no_resp_go_ind = [intersect(noresp_nofix_ind,go_event_ind)];   
  for i=no_resp_go_ind
    all_types_new{i} = 'NoRespGo';
  end
  % incorrect_go
  incorrect_go_ind = [];
  incorrect_go_late_ind = [];
  for i=1:size(go_cresp,2)  
    if go_resp{i} ~= go_cresp{i}
      incorrect_go_ind = [incorrect_go_ind i]; 
      all_types_new{i} = 'IncorrectGo'; %legit 1~=2
    elseif go_resp{i} == go_cresp{i}
    else 
      if fix_resp{i} ~= go_cresp{i}
        incorrect_go_late_ind = [incorrect_go_late_ind i]; 
        all_types_new{i} = 'IncorrectLateGo';
      end 
    end 
  end
  % correct & incorrect _stop 
  ssd_rt = {event_info.ssd_rt};
  ssd_dur = {event_info.ssd_dur};
  stop_resp = {event_info.stop_resp};  
  fix_resp = {event_info.fix_resp}; 
  incorrect_stop_ind = [];
  correct_stop_ind = []; 
  ssd_stop_ind = []; 
  for i=1:size(ssd_dur,2) 
    if ssd_dur{i} > 0  % only stop trials  
      if stop_resp{i} > 0
        incorrect_stop_ind = [incorrect_stop_ind i]; 
        all_types_new{i} = 'IncorrectStop';
      elseif fix_resp{i} > 0
        incorrect_stop_ind = [incorrect_stop_ind i]; 
        all_types_new{i} = 'IncorrectStop';     
      else 
        correct_stop_ind = [correct_stop_ind i];
        all_types_new{i} = 'CorrectStop';
      end
      if ssd_rt{i} > 0
        % stop too early 
        ssd_stop_ind = [ssd_stop_ind i];   
        all_types_new{i} = 'SSDStop';
      end
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ind_start,event_ref_time] = set_ref(onset,start_time)
  if ~isempty(onset)
    d = bsxfun(@minus,onset,start_time');
    d(d<0) = Inf;
    [d,ind_start] = min(d,[],1); 
    event_ref_time  = start_time(ind_start);
  else
    ind_start = start_time;
    event_ref_time = 0;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [onset,offset] = check_offsets(onset,offset,eventname,parms)
  ind_empty = find(cellfun(@isempty,onset) | cellfun(@isempty,offset));
  if ~isempty(ind_empty)
    if parms.verbose, fprintf('%s: WARNING: %s event has %d onsets and %d offsets\n',...
      mfilename,eventname,nnz(~cellfun(@isempty,onset)),nnz(~cellfun(@isempty,offset))); end
    ind_keep = setdiff([1:length(onset)],ind_empty);
    onset = onset(ind_keep);
    offset = offset(ind_keep);
  end;
  onset = cell2mat(onset);
  offset = cell2mat(offset);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_files(eventname,onset,offset,ind_start,event_ref_time,parms)
 nscans = parms.eprime_nruns;
 scans = parms.eprime_runs;  
  
  for s=scans
    ind_scan = find(ind_start == s);
    fname_1D = sprintf('%s/%s_scan%d_%s.1D',...
      parms.outdir,parms.outstem,s,eventname);
    fname_txt = sprintf('%s/%s_scan%d_%s.txt',...
      parms.outdir,parms.outstem,s,eventname);
    fname_block = sprintf('%s/%s_scan%d_%s_block.txt',...
      parms.outdir,parms.outstem,s,eventname);  
    fname_block_dur = sprintf('%s/%s_scan%d_%s_block_dur.txt',...
      parms.outdir,parms.outstem,s,eventname);  
    fname_fsl = sprintf('%s/%s_scan%d_%s_fsl.txt',...
      parms.outdir,parms.outstem,s,eventname);
    if ~exist(fname_1D,'file') || ~exist(fname_txt,'file') || ...
       ~exist(fname_block,'file') || ~exist(fname_fsl,'file') || parms.forceflag
      if ~isempty(onset(ind_scan))
        % subtract start time and convert to seconds
        rel_onset  = (onset(ind_scan) - event_ref_time(ind_scan))/1000;
        rel_offset = (offset(ind_scan) - event_ref_time(ind_scan))/1000;
        % loop over events
        TR_mask = zeros(parms.numTRs,1);
        for k=1:length(rel_onset)
          ind_rep = find(rel_onset(k)<(parms.TR_offset-parms.mindur) &...
                         rel_offset(k)>(parms.TR_onset+parms.mindur));
          TR_mask(ind_rep) = 1;
        end;
        % calculate event times from TR_onset
        event_times = parms.TR_onset(find(TR_mask));
        % convert time units to multiples of TRs
        event_TRs = event_times/parms.TR;
        % remove  nskipTRs
        if parms.nskipTRs>0
          TR_mask = TR_mask(parms.nskipTRs+1:end);
          event_times = event_times(parms.nskipTRs+1:end) - parms.nskipTRs*parms.TR;
          event_TRs = event_times/parms.TR; 
        end;
        % write to 1D file
        fid = fopen(fname_1D,'wt');
        if fid<0, error('failed to open %s for writing',fname_1D); end;
        for k=1:length(TR_mask)
          fprintf(fid,'%d\n',TR_mask(k));
        end;
        fclose(fid);
        % write to txt file with stim times in seconds
        fid = fopen(fname_txt,'wt');
        if fid<0, error('failed to open %s for writing',fname_txt); end;
        fprintf(fid,'%s\n',sprintf('%0.2f ',event_TRs*parms.TR));
        fclose(fid);
        % write to block txt file
        fid = fopen(fname_block,'wt');  
        if fid<0, error('failed to open %s for writing',fname_block); end;    
        fprintf(fid,'%s\n',sprintf('%0.2f ',rel_onset));
        fclose(fid);
        fid = fopen(fname_block_dur,'wt');  
        if fid<0, error('failed to open %s for writing',fname_block_dur); end;
        block_len = mean(rel_offset- rel_onset); 
        fprintf(fid,'%s\n',sprintf('%0.2f ',block_len));
        fclose(fid);
        % write to fsl file
        fid = fopen(fname_fsl,'wt');
        if fid<0, error('failed to open %s for writing',fname_fsl); end;
        for i=1:length(rel_onset)
          fprintf(fid, '%0.2f %0.2f 1 \n',rel_onset(i), rel_offset(i)-rel_onset(i));
        end
        fclose(fid); 
      else
        TR_mask = zeros(parms.numTRs,1);
        % write to 0s 1D file
        fid = fopen(fname_1D,'wt');
        if fid<0, error('failed to open %s for writing',fname_1D); end;
        for k=1:length(TR_mask)
          fprintf(fid,'%d\n',TR_mask(k));
        end;
        fclose(fid);
        % write * txt file
        fid = fopen(fname_txt,'wt');
        if fid<0, error('failed to open %s for writing',fname_txt); end;
        fprintf(fid,'*\n');
        fclose(fid);   
        % write * txt file block 
        fid = fopen(fname_block,'wt');
        if fid<0, error('failed to open %s for writing',fname_block); end;
        fprintf(fid,'*\n');
        fclose(fid);  
        fid = fopen(fname_block_dur,'wt');
        if fid<0, error('failed to open %s for writing',fname_block); end;
        fprintf(fid,'0\n');
        fclose(fid); 
        % write empty fsl file
        fid = fopen(fname_fsl,'wt');
        if fid<0, error('failed to open %s for writing',fname_fsl); end;
        fclose(fid);    
      end
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [behav,runs_ok,errcode,errmsg] = get_behavioral_data(event_info,parms,all_types_new) 
  errcode = 0; errmsg = [];

  behav = []; 
  behav.('SubjID') = []; behav.('VisitID') = []; 
  behav.version = mmil_getfield(event_info(1),'version','UNKNOWN');
  if parms.verbose
    fprintf('%s: experiment version = %s\n',mfilename,behav.version);
  end
  behav.switch_flag = parms.switch_flag;
  behav.perform_flag = 1;

  blocks = [event_info.block]; 
  trials = [event_info.trial]; 
  % determine version of task
  n_blocks = size(unique(blocks),2);
  n_trials = size(unique(trials),2);
  if n_trials > 2, 
    blocks = [event_info.block]; 
  else 
    blocks = [event_info.trial];
  end
  
  % check if second runs is truncated 
  runs = blocks; 
  nruns = length(unique(blocks)); 
  runs_ok = []; 
  for i=1:nruns
    run_len = length(find(runs==i));     
    if run_len == 182 % hardcode the run length for now
     runs_ok = [runs_ok i];  
    else 
     if parms.verbose, fprintf('%s: WARNING: run %d is short: %d trials found (<182 trials) \n',mfilename,i,run_len); end
    end 
  end
  nruns = length(runs_ok); 
  
  if nruns == 1
    new_info = find(runs==runs_ok);
    all_types_new = all_types_new(new_info);
    blocks = blocks(new_info); 
    event_info = event_info(new_info);
    runs = runs(new_info);
  elseif nruns == 0
   if parms.verbose, fprintf('%s: ERROR: no valid e-prime runs in %s\n',mfilename,parms.fname); end
   errcode = 1;
   errmsg = 'no valid e-prime runs';
   return;
  end  
  behav.nruns = nruns; 
  
  % event durations
  go_dur = {event_info.go_dur}; 
  time_go_dur = [event_info.go_dur];
  time_go_dur = time_go_dur(1); 
  stop_dur = {event_info.stop_dur}; 
  go_event_ind = find(~cellfun(@isempty,go_dur)); 
  stop_event_ind = find(~cellfun(@isempty,stop_dur));
  % total numbers of trials    
  behav.('Go_total') = size(go_event_ind,2);
  behav.('Stop_total') = size(stop_event_ind,2);
  behav.('Total') = size(go_event_ind,2)+size(stop_event_ind,2);
  
  % number of trials per run
  for j=1:2 %runs 
    % go
    type_counts_run_go = sprintf('Go_total_run_%d',j); 
    behav.(type_counts_run_go) = 'NaN';  
    % stop
    type_counts_run_stop = sprintf('Stop_total_run_%d',j);
    behav.(type_counts_run_stop) = 'NaN';  
  end
  for j=runs_ok %runs 
    % go
    type_counts_run_go = sprintf('Go_total_run_%d',j);
    ind_type_run = [intersect(find(blocks==j),go_event_ind)];   
    behav.(type_counts_run_go) = size(ind_type_run,2); 
    % stop
    type_counts_run_stop = sprintf('Stop_total_run_%d',j);
    ind_type_run = [intersect(find(blocks==j),stop_event_ind)];   
    behav.(type_counts_run_stop) = size(ind_type_run,2); 
  end
  
  % results for each condition
  for i=1:parms.nconds
    % total counts per type
    type  = parms.typenames{i};         
    type_counts_total = sprintf('%s_counts_total',type);  
    ind_type = find(strcmp(type,all_types_new)); 
    behav.(type_counts_total) = size(ind_type,2); 
    % total counts per type per run  
    for j=1:2 % runs        
      type_counts_run = sprintf('%s_counts_run_%d',type,j);   
      behav.(type_counts_run) = 'NaN'; 
    end 
    for j=runs_ok       
      type_counts_run = sprintf('%s_counts_run_%d',type,j); 
      ind_type_run = [intersect(find(blocks==j),ind_type)];  
      behav.(type_counts_run) = size(ind_type_run,2); 
    end 
    
    % selecting go or stop event 
    cond  = parms.condnames{i};  
    pat = '\w_(?<Event>\w+)'; 
    get_event = regexp(cond,pat,'names');
    event = get_event.Event;    
    switch event  
      case 'go'
        event_ind = go_event_ind; 
        event_rt = {event_info.('go_rt')}; 
        if strcmp(cond,'correctlate_go')
           event_rt_fix = {event_info.('fix_rt')};
           event_rt = cellfun(@(x) x+time_go_dur, event_rt_fix,'UniformOutput',false);
        end 
        if strcmp(cond,'incorrectlate_go')
           event_rt_fix = {event_info.('fix_rt')};
           event_rt = cellfun(@(x) x+time_go_dur, event_rt_fix,'UniformOutput',false);
        end 
      case 'stop'
        event_ind = stop_event_ind;
        event_rt = {event_info.('stop_rt')};
    end 
    % total percent per type 
    type_percent_total = sprintf('%s_percent_total',type);  
    behav.(type_percent_total) = size(ind_type,2)./size(event_ind,2);  
    % percent per type per run  
    for j=1:2 % runs 
      type_counts_run = sprintf('%s_percent_run_%d',type,j); 
      behav.(type_counts_run) = 'NaN';  
    end
    for j=runs_ok 
      type_counts_run = sprintf('%s_percent_run_%d',type,j); 
      ind_type_run = [intersect(find(blocks==j),ind_type)];
      ind_type_percent_run = [intersect(find(blocks==j),event_ind)]; 
      behav.(type_counts_run) = size(ind_type_run,2)./size(ind_type_percent_run,2);  
    end
    % total rt per type 
    type_rt_total = sprintf('%s_rt_total',type); 
    all_rt = cell2mat(event_rt(ind_type)); 
    behav.(type_rt_total) = mean(all_rt(all_rt~=0));  
    type_rt_std_total = sprintf('%s_rt_std_total',type);
    behav.(type_rt_std_total) = std(all_rt(all_rt~=0)); 
    
    % total rt per type per run 
    for j=1:2 % runs 
      type_rt_runs = sprintf('%s_rt_run_%d',type,j);
      behav.(type_rt_runs) = 'NaN';  
      type_rt_std_runs = sprintf('%s_rt_std_run_%d',type,j);
      behav.(type_rt_std_runs) = 'NaN';
    end    
    for j=runs_ok % runs 
      type_rt_runs = sprintf('%s_rt_run_%d',type,j);
      ind_type_run = [intersect(find(blocks==j),ind_type)];
      all_rt = cell2mat(event_rt(ind_type_run));  
      behav.(type_rt_runs) = mean(all_rt(all_rt~=0)); 
      type_rt_std_runs = sprintf('%s_rt_std_run_%d',type,j);
      behav.(type_rt_std_runs) = std(all_rt(all_rt~=0));
    end  
    
    
  end
  % ssd and ssrt 
  ssd = {event_info.ssd_dur}; 
  behav.('SSD_mean_total')= mean(cell2mat(ssd));   
  behav.('SSRT_mean_total') = behav.CorrectGo_rt_total-behav.SSD_mean_total; 
  % rm meaningless fields 
  rmfields_names = {'NoRespGo_rt_total','NoRespGo_rt_run_1','NoRespGo_rt_run_2', ... 
                    'NoRespGo_rt_std_total','NoRespGo_rt_std_run_1','NoRespGo_rt_std_run_2', ... 
                    'CorrectStop_rt_total','CorrectStop_rt_run_1', ...
                    'CorrectStop_rt_std_total','CorrectStop_rt_std_run_1', ...
                    'CorrectStop_rt_run_2','SSDStop_rt_total','SSDStop_rt_run_1',...
                    'CorrectStop_rt_std_run_2','SSDStop_rt_std_total','SSDStop_rt_std_run_1',...
                    'SSDStop_rt_run_2','SSDStop_rt_std_run_2'}; 
                  
  behav = rmfield(behav,rmfields_names); 
  
  %behav flag 
  if behav.CorrectGo_percent_total < 0.6   
    if parms.verbose, fprintf('%s: CorrectGo_percent_total < 60%% \n',mfilename); end
    behav.perform_flag = 0; 
  end 
  if behav.IncorrectGo_percent_total > 0.3
    if parms.verbose, fprintf('%s: IncorrectGo_percent_total > 30%% \n',mfilename); end
    behav.perform_flag = 0; 
  end 
  late_go = behav.CorrectLateGo_percent_total+behav.IncorrectLateGo_percent_total;
  if late_go > 0.3 
    if parms.verbose, fprintf('%s: LateGo > 30%% \n',mfilename); end
    behav.perform_flag = 0; 
  end
  if behav.NoRespGo_percent_total > 0.3
    if parms.verbose, fprintf('%s: NoRespGo_percent_total > 30%% \n',mfilename); end
    behav.perform_flag = 0; 
  end 
  if behav.CorrectStop_percent_total < 0.2 || behav.CorrectStop_percent_total > 0.8
    if parms.verbose, fprintf('%s: CorrectStop_percent_total < 20%% or > 80%% \n',mfilename); end
    behav.perform_flag = 0; 
  end 
  if behav.Go_total ~= 300
    if parms.verbose, fprintf('%s: Go_total trials not equal to 300 \n',mfilename); end
    behav.perform_flag = 0; 
  end
  if behav.CorrectGo_rt_total < behav.IncorrectStop_rt_total
    if parms.verbose, fprintf('%s: CorrectGo_rt_total < IncorrectStop_rt_total \n',mfilename);end
    behav.perform_flag = 0;
  end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% NOTE: section below by Octavio Ruiz, 22jan12
  %%  Based on python script eprime_funcs.py by Sage Hahn from Hugh Garavan's lab at UVM
  %%  Reproduces the functions called by the loop inside the function process_event(event, files)
  %%
  %% todo: move section below to new subfunction
  %%       check for redundancy (does it calculate variables that have already been set?)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  behav.SSRT_integrated_total = NaN;
  behav.violatorflag = NaN;
  behav.glitchflag = NaN;
  behav.zeroSSDcount = NaN;
  
  tfmri_sst_all_beh_total_issrt = NaN;
  tfmri_sst_beh_glitchflag = NaN;
  tfmri_sst_beh_glitchcnt = NaN;
  tfmri_sst_beh_0SSDcount = NaN;
  tfmri_sst_beh_0SSD_flag = NaN;  % called "tfmri_sst_beh_0SSD>20flag" in Sage's script
  tfmri_sst_beh_violatorflag = NaN;
  
  % set_correct_go
  manynans = num2cell(NaN(1,length(event_info)));
  [event_info.('correct_go_response')] = manynans{:};

  inds = cellfun(@(x,y)(~isempty(x) & (x==y)),...
                  {event_info.go_resp}, {event_info.go_cresp}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
%   fprintf('num_correct_go_response_1 = %.0f\n', length(inds))
  %% todo: is there a way to do this without looping? (deal?)
  for i = 1:length(inds)
    j = inds(i);
    event_info(j).correct_go_response = 1.0;
  end
  
  inds = cellfun(@(x,y,z)( isempty(x) & (y==z)),...
                  {event_info.go_resp}, {event_info.go_cresp}, {event_info.fix_resp},'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
%   fprintf('num_correct_go_response_1 = %.0f\n', length(inds))
  for i = 1:length(inds)
    j = inds(i);
    event_info(j).correct_go_response = 1.0;
  end

  inds = cellfun(@(x,y,z)(~isempty(x) & ~isequal(x,y) & ~isempty(strfind(z,'Go'))),...
                  {event_info.go_resp}, {event_info.go_cresp}, {event_info.type}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
%   fprintf('num_correct_go_response_0 = %.0f\n', length(inds))
  for i = 1:length(inds)
    j = inds(i);
    event_info(j).correct_go_response = 0.0;
  end
  
  inds = cellfun(@(x,y,z,r)( isempty(x) & ~isequal(y,z) & ~isempty(strfind(r,'Go'))),...
                  {event_info.go_resp}, {event_info.go_cresp}, {event_info.fix_resp}, {event_info.type},...
                  'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
%   fprintf('num_correct_go_response_0 = %.0f\n', length(inds))
  for i = 1:length(inds)
    j = inds(i);
    event_info(j).correct_go_response = 0.0;
  end

  inds = cellfun(@(x,y,z)( isempty(x) & isempty(y) & ~isempty(strfind(z,'Go'))),...
                  {event_info.go_resp}, {event_info.fix_resp}, {event_info.type}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
%   fprintf('num_correct_go_response_omission = %.0f\n', length(inds))
  for i = 1:length(inds)
    j = inds(i);
    event_info(j).correct_go_response = 'omission';
  end

  % set_correct_stop
  manynans = num2cell(NaN(1,length(event_info)));
  [event_info.('correct_stop')] = manynans{:};

  inds = cellfun(@(x,y,z,r)( isempty(x) & isempty(y) & isempty(z) & ~isempty(strfind(r,'Stop')) ),...
            {event_info.stop_resp}, {event_info.fix_resp}, {event_info.ssd_resp}, {event_info.type}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  for i = 1:length(inds)
    j = inds(i);
    event_info(j).correct_stop = 1.0;
  end

  inds = cellfun(@(x,y,z,r)((~isempty(x) | ~isempty(y) | ~isempty(z)) & ~isempty(strfind(r,'Stop')) ),...
            {event_info.stop_resp}, {event_info.fix_resp}, {event_info.ssd_resp}, {event_info.type}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  for i = 1:length(inds)
    j = inds(i);
    event_info(j).correct_stop = 0.0;
  end

  % set_correct_go_rt
  inds = cellfun(@(x,y,z)(~isempty(x) & isempty(y) & ~isempty(strfind(z,'Go')) ),...
                  {event_info.fix_resp}, {event_info.go_resp}, {event_info.type}, 'UniformOutput',0 );
  inds = find(cellfun(@(x)(isequal(x,1)), inds));

  [event_info.('go_rt_adjusted')] = event_info.go_rt;
  for i = 1:length(inds)
    j = inds(i);
    event_info(j).go_rt_adjusted = event_info(j).go_dur + event_info(j).fix_rt;
  end

  % set_correct_stop_rt
  [event_info.('stop_rt_adjusted')] = event_info.stop_rt;

  % set to Stop signal duration + fix rt, as was answered during fix
  inds = cellfun(@(x,y,z,r)(~isempty(x) & isempty(y) & ~isempty(strfind(z,'Stop')) & isequal(r,0) ),...
            {event_info.fix_resp}, {event_info.stop_resp}, {event_info.type}, {event_info.correct_stop}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  for i = 1:length(inds)
    j = inds(i);
    event_info(j).stop_rt_adjusted = event_info(j).stop_dur + event_info(j).fix_rt;
  end

  % adjust for answers during SSD
  inds = cellfun(@(x)(~isempty(x)), {event_info.ssd_resp}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  for i = 1:length(inds)
    j = inds(i);
    event_info(j).stop_rt_adjusted = event_info(j).ssd_rt;
  end

  % set updated for response during stop
  inds = cellfun(@(x,y,z)((~isempty(x) | ~isempty(y)) & ~isempty(z)),...
            {event_info.stop_resp}, {event_info.fix_resp}, {event_info.stop_rt_adjusted}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  for i = 1:length(inds)
    j = inds(i);
    event_info(j).stop_rt_adjusted = event_info(j).stop_rt_adjusted + event_info(j).ssd_dur;
  end

  % check_omissions
  inds = cellfun(@(x,y,z)( isempty(x) & isempty(y) & ~isempty(strfind(z,'Go')) ),...
            {event_info.go_resp}, {event_info.fix_resp}, {event_info.type}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  n_missing = length(inds);

  inds = cellfun(@(z)( ~isempty(strfind(z,'Go')) ), {event_info.type}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  n_go_trials = length(inds);

  if n_missing == n_go_trials
    disp('- All trials in this file are ommissions -')
  end

  %% todo: make below a subfunction?

  % separate just go trials
  inds = cellfun( @(x)( ~isempty(strfind(x,'Go')) ), {event_info.type}, 'UniformOutput',0 );
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  go_trials = event_info(inds);

  % set omissions to max go.rt - if primary rt is_null means omission
  inds = cellfun( @(x)( strcmp(x,'omission') ), {go_trials.correct_go_response}, 'UniformOutput',0 );
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  go_rt_max = max( [go_trials.go_rt] );
  for i = 1:length(inds)
    j = inds(i);
    go_trials(j).go_rt_adjusted = go_rt_max;
  end
  % sort go trials
  sorted_go = sort([go_trials.go_rt_adjusted]);

  % get stop trials
  inds = cellfun( @(x)(~isempty(strfind(x,'Stop'))), {event_info.type}, 'UniformOutput',0 );
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  stop_trials = event_info(inds);

  % calculate prob stop failure
  prob_stop_failure = 1 - ...
    sum( cellfun(@(x)(isequal(x,1)),{stop_trials.correct_stop})) / length( stop_trials );

  % calc go.rt index
  index = prob_stop_failure * length(sorted_go);

  % if prob_stop_failure is 1, use the max go.rt
  if ceil(index) == length(sorted_go)
    index = length(sorted_go);
  else
    index = [floor(index) ceil(index)]+1;
  end

  % calc SSRT
  mean_ssd = mean( [stop_trials.ssd_dur] );
  try
    ssrt = mean(sorted_go(index)) - mean_ssd;
  catch
      fprintf('index = %f\n', index);
      fprintf('len sorted go = %f\n', length(sorted_go));
      fprintf('prob stop fail = %f\n', prob_stop_failure);
      ssrt = NaN;
  end
  tfmri_sst_all_beh_total_issrt = ssrt;
  if parms.verbose, fprintf('%s: tfmri_sst_all_beh_total_issrt = %.4f\n',...
      mfilename,tfmri_sst_all_beh_total_issrt); end

  % number of glitched trials
  inds = arrayfun(@(x,y)( (x < 50) & (x > 0) & (y <= 50)), [event_info.ssd_rt], [event_info.ssd_dur]);
  tfmri_sst_beh_glitchcnt = sum(inds);
  if parms.verbose, fprintf('%s: tfmri_sst_beh_glitchcnt = %.0f\n',...
      mfilename,tfmri_sst_beh_glitchcnt); end
  
  tfmri_sst_beh_glitchflag = (tfmri_sst_beh_glitchcnt > 0);
  if parms.verbose, fprintf('%s: tfmri_sst_beh_glitchflag = %.0f\n',...
      mfilename,tfmri_sst_beh_glitchflag); end

  %  def get_0SSD_cnt(subj)  and   get_0SSD_flag(subj):
  inds = arrayfun(@(x)( x == 0 ), [event_info.ssd_dur]);
  tfmri_sst_beh_0SSDcount = sum(inds);
  if parms.verbose, fprintf('%s: tfmri_sst_beh_0SSDcount = %.0f\n',...
      mfilename,tfmri_sst_beh_0SSDcount); end
  
  tfmri_sst_beh_0SSD_flag = (tfmri_sst_beh_0SSDcount > 20);  % Caled tfmri_sst_beh_0SSD>20flag in Sage's script
  if parms.verbose, fprintf('%s: tfmri_sst_beh_0SSD_flag = %.0f\n',...
      mfilename,tfmri_sst_beh_0SSD_flag); end

  inds = cellfun(@(x,y)( (isequal(x,0) | isequal(x,1))  & ~isempty(strfind(y,'Go')) ),...
            {event_info.correct_go_response}, {event_info.type}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  go_rt = mean( [event_info(inds).go_rt_adjusted] );

  
  inds = cellfun(@(x)( isequal(x,0) ), {event_info.correct_stop}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  stop_rt = mean( [event_info(inds).stop_rt_adjusted] );
  
  if parms.verbose
    fprintf('%s:   go_rt = %.4f\n', mfilename,go_rt);
    fprintf('%s: stop_rt = %.4f\n', mfilename,stop_rt);
  end
  
  tfmri_sst_beh_violatorflag = (stop_rt > go_rt);
  if parms.verbose, fprintf('%s: tfmri_sst_beh_violatorflag = %.0f\n',...
      mfilename,tfmri_sst_beh_violatorflag); end

  behav.SSRT_integrated_total = tfmri_sst_all_beh_total_issrt;
  behav.violatorflag = tfmri_sst_beh_violatorflag;
  behav.glitchflag = tfmri_sst_beh_glitchflag;
  behav.zeroSSDcount  = tfmri_sst_beh_0SSDcount;

  % write to csv
  fname_csv_out = sprintf('%s/%s_behavioral.csv',parms.outdir,parms.outstem); 
  if ~exist(fname_csv_out,'file') || parms.forceflag
    mmil_struct2csv(behav,fname_csv_out)
  end;
return; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function get_behavioral_data_empty(parms) 

  behav = []; 
  behav.('SubjID') = []; behav.('VisitID') = []; 
  behav.version = [];
  behav.switch_flag = 0; 
  behav.perform_flag = 0;
  behav.nruns = 0; 
  
  % total numbers of trials    
  behav.('Go_total') = 'NaN'; 
  behav.('Stop_total') = 'NaN'; 
  behav.('Total') = 'NaN'; 
  % number of trials per run    
  for j=1:2 %runs 
    % go
    type_counts_run_go = sprintf('Go_total_run_%d',j);
    behav.(type_counts_run_go) = 'NaN'; 
    % stop
    type_counts_run_stop = sprintf('Stop_total_run_%d',j);  
    behav.(type_counts_run_stop) = 'NaN';  
  end
  
  % results for each condition
  for i=1:parms.nconds
    % total counts per type
    type  = parms.typenames{i};         
    type_counts_total = sprintf('%s_counts_total',type);  
    behav.(type_counts_total) = 'NaN'; 
    % total counts per type per run  
    for j=1:2 % runs        
      type_counts_run = sprintf('%s_counts_run_%d',type,j);  
      behav.(type_counts_run) = 'NaN'; 
    end 

    % total percent per type 
    type_percent_total = sprintf('%s_percent_total',type);  
    behav.(type_percent_total) = 'NaN'; 
    % percent per type per run  
    for j=1:2 % runs 
      type_counts_run = sprintf('%s_percent_run_%d',type,j); 
      behav.(type_counts_run) = 'NaN'; 
    end
    % total rt per type 
    type_rt_total = sprintf('%s_rt_total',type); 
    behav.(type_rt_total) = 'NaN'; 
    type_rt_std_total = sprintf('%s_rt_std_total',type);
    behav.(type_rt_std_total) = 'NaN'; 
    
    % total rt per type per run 
    for j=1:2 % runs 
      type_rt_runs = sprintf('%s_rt_run_%d',type,j);  
      behav.(type_rt_runs) = 'NaN'; 
      type_rt_std_runs = sprintf('%s_rt_std_run_%d',type,j);
      behav.(type_rt_std_runs) = 'NaN'; 
    end     
  end 
  % ssd and ssrt 
  behav.('SSD_mean_total')= 'NaN'; 
  behav.('SSRT_mean_total') = 'NaN';  
  % rm meaningless fields 
  rmfields_names = {'NoRespGo_rt_total','NoRespGo_rt_run_1','NoRespGo_rt_run_2', ... 
                    'NoRespGo_rt_std_total','NoRespGo_rt_std_run_1','NoRespGo_rt_std_run_2', ... 
                    'CorrectStop_rt_total','CorrectStop_rt_run_1', ...
                    'CorrectStop_rt_std_total','CorrectStop_rt_std_run_1', ...
                    'CorrectStop_rt_run_2','SSDStop_rt_total','SSDStop_rt_run_1',...
                    'CorrectStop_rt_std_run_2','SSDStop_rt_std_total','SSDStop_rt_std_run_1',...
                    'SSDStop_rt_run_2','SSDStop_rt_std_run_2'}; 
                  
  behav = rmfield(behav,rmfields_names); 
  % write to csv
  fname_csv_out = sprintf('%s/%s_behavioral.csv',parms.outdir,parms.outstem); 
  if ~exist(fname_csv_out,'file') || parms.forceflag
    mmil_struct2csv(behav,fname_csv_out)
  end;
  
return; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
