function [eprime_nruns,eprime_runs,errcode,behav,errmsg] = abcd_extract_eprime_nback(fname,varargin)
%function [eprime_nruns,eprime_runs,errcode,behav,errmsg] = abcd_extract_eprime_nback(fname,[options])
%
% Purpose: extract condition time courses
%   from eprime data files for NBACK task
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
%   'error_flag': [0|1] encode error trials as separate events
%     {default = 0}
%   'extra_flag': [0|1] create timing files for extra, orthogonal event encoding
%     (i.e. 'target','lure','nonlure')
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%   'timing_files_flag': [0|1] generate timing files
%     {default = 1}
%   'verbose': [0|1] display messages
%     {default = 1}
%
%  Output: 
%    eprime_nruns: number of valid runs in the e-prime in the file  
%    errcode: [0|1] whether the file was successfully processed
%    behav: behavioral data
%    errmsg: string describing error if errcode=1
%
% Created : 01/06/17 by Jose Teruel 
% Prev Mod: 01/23/19 by Dani Cornejo
% Prev Mod: 01/29/20 by Don Hagler
% Prev Mod: 05/13/20 by Octavio Ruiz
% Prev Mod: 05/03/21 by Don Hagler
% Prev Mod: 08/18/21 by Don Hagler
% Prev Mod: 02/10/22 by Octavio Ruiz
% Prev Mod: 04/22/22 by Don Hagler
% Prev Mod: 06/14/24 by Don Hagler
% Last Mod: 06/24/24 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: based on Nback_scanner.py provided by Eric Feczko from OHSU
%       Using style from abcd_extract_eprime_mid.m provided by Donald Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize outputs
eprime_nruns = []; eprime_runs = []; errcode = 0; behav = []; errmsg = [];

% check arguments 
if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
parms = check_input(fname,varargin); 

% create output directory
mmil_mkdir(parms.outdir);

% create struct with info for each event
[event_info,start_time,all_types,all_stims,all_targets,all_procs,errcode,errmsg] = get_event_info(parms);
if errcode, return; end;

% check for no responses 
[nresp, errcode, errmsg] = check_no_responses(event_info,parms);
if errcode, return; end;
        
% switch buttons if necesary 
[event_info,parms.switch_flag,errcode,errmsg] = nback_switch(event_info,parms); 
if errcode, return; end;

% remove non-events
[event_info,event_info_proc] = remove_non_events(event_info);

% get behav data and write it to a csv 
[behav,eprime_runs,errcode,errmsg] = get_behavioral_data_nback(event_info,parms,start_time);
if errcode, return; end;

parms.eprime_runs = eprime_runs; 
eprime_nruns = length(eprime_runs);
parms.eprime_nruns = eprime_nruns; 
if parms.verbose, fprintf('%s: number of e-prime runs = %d \n',mfilename,eprime_nruns); end

if parms.timing_files_flag && eprime_nruns 
  % write files for each condition
  for i=1:parms.nconds
    type = parms.typenames{i}; 
    cond = parms.condnames{i}; 
    ind_type = find(strcmp(type,all_types));
    for j=1:parms.nstims
      for corr_flag = parms.corr_flags
        stim = parms.stimnames{j};
        eventname = sprintf('%s_%s',cond,stim);  
        if (corr_flag==0), eventname = [eventname '_err']; end;
        ind_stim = find(strcmp(stim,lower(all_stims)));
        if isempty(ind_stim)
          ind_stim = find(strcmp(stim,lower(all_targets)));
        end
        ind_inter = intersect(ind_type,ind_stim);
        if parms.error_flag
          acc = [event_info.stim_acc];
          ind_acc = find(acc==corr_flag);
          ind_inter = intersect(ind_inter,ind_acc);
        end;
        % get response_times
        resp_times = {event_info(ind_inter).stim_rt};
        ind_noresp = find(cellfun(@isempty,{event_info(ind_inter).stim_resp}));
        resp_times(ind_noresp) = {NaN};
        % get times of stim onset and offset
        onset = {event_info(ind_inter).stim_onset}; 
        offset = {event_info(ind_inter).stim_offset};
        % check offsets match onsets
        [onset,offset,resp_times] = check_offsets(onset,offset,resp_times,eventname,parms);
        % find most recent start time for each event
        [ind_start,event_ref_time] = set_ref(onset,start_time); 
        % create files for each scan
        write_files(eventname,onset,offset,resp_times,ind_start,event_ref_time,parms);
      end;
    end;
  end

  % create timing files for cue trials (pre-block instructions)
  onset = []; offset = []; resp_times = [];
  eventname = 'cue';
  for i=1:parms.ncues
    ind_cues = find(strcmp(parms.cues{i},all_procs));
    % onset of fixation block cue
    tmp_onset = {event_info_proc(ind_cues).cuefix_onset};
    % offset of stimulus block cue
    tmp_offset = {event_info_proc(ind_cues).([parms.cuenames{i} '_offset'])};
    % set tmp_resp_times to all NaNs (because no response to cue)
    tmp_resp_times = repmat({NaN},size(tmp_onset));
    % check offsets match onsets
    [tmp_onset,tmp_offset,tmp_resp_times] = check_offsets(tmp_onset,tmp_offset,tmp_resp_times,eventname,parms);
    % concatenate with other cue trials    
    onset = [onset tmp_onset];
    offset = [offset tmp_offset];
    resp_times = [resp_times tmp_resp_times];
  end

  % sort onset and offset times
  [onset,ind_sort] = sort(onset); 
  offset = offset(ind_sort);
  resp_times = resp_times(ind_sort);

  % find most recent start time for each event
  [ind_start,event_ref_time] = set_ref(onset,start_time);
 
  % create files for each scan
  write_files(eventname,onset,offset,resp_times,ind_start,event_ref_time,parms);

end %timing_files_flag

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
    'error_flag',false,[false true],...
    'extra_flag',false,[false true],...
    'forceflag',false,[false true],...
    'timing_files_flag',true,[false true],...
    'verbose',true,[false true],...
    ...
    'eprime_runs',1:2,[],...
    'eprime_nruns',0,0:2,...
    'expected_numtrials',80,[],...
    'expected_start_delay_GE',12,[],...
    'expected_start_delay_Siemens',6.4,[],...
    'expected_init_delay_GE_run1',0.5,[],...
    'expected_init_delay_GE_run2',0,[],...
    'expected_init_delay_Siemens',0,[],...
    'start_delay_Siemens',6.4,[],...
    'max_start_diff',0.5,[],...
    'max_init_diff',5.0,[],...
    ...
    'colnames',  {'NARGUID','SessionDate','SessionTime','ExperimentName','ExperimentVersion',...
                  'Procedure[Block]','BlockType','StimType',...
                  'GetReady.RTTime','GetReady2.RTTime','Stim.OnsetTime','Stim.OffsetTime','Stim.ACC',...
                  'Cue2Back.OnsetTime','Cue2Back.OffsetTime','CueTarget.OnsetTime','CueTarget.OffsetTime',...
                  'CueFix.OnsetTime','CueFix.OffsetTime','CueFix.StartTime',...
                  'Fix15sec.OnsetTime', 'Fix15sec.OffsetTime','TargetType','Stim.RT',...
                  'CorrectResponse','Stim.RESP', 'Stimulus'},[],...
    'fieldnames',{'narguid','date','time','experiment','version',...
                  'procedure_block','block_type','stim_type',...
                  'getready_rttime','getready2_rttime','stim_onset','stim_offset','stim_acc',...
                  'cue2back_onset','cue2back_offset','cue0back_onset','cue0back_offset',...
                  'cuefix_onset','cuefix_offset','cuefix_start',...
                  'fixation_onset', 'fixation_offset','target_type','stim_rt',...
                  'correct_response','stim_resp', 'stim'},[],...
    'typenames',{'2-Back','0-Back'},[],...
    'condnames',{'2_back','0_back'},[],...
    'stimnames',{'posface','neutface','negface','place'},[],...
    'extra_stimnames',{'target','lure','nonlure'},[],...
    'cues',{'Cue0BackPROC','Cue2BackPROC'},[]...
    'cuenames',{'cue0back','cue2back'},[],...
    'procedures',{'TRSyncPROC', 'TRSyncPROCR2'},[]...
  });

  %% todo: option to specify which type(s) of file to write?
  
  if parms.error_flag
    parms.corr_flags = [0,1];
  else
    parms.corr_flags = 2;
  end;
  if parms.extra_flag
    parms.stimnames = cat(2,mmil_rowvec(parms.stimnames),...
                            mmil_rowvec(parms.extra_stimnames));
  end;
  parms.ncues = length(parms.cues);
  parms.nprocedures = length(parms.procedures);  
  parms.nconds = length(parms.condnames);
  parms.ntypes = length(parms.typenames);
  parms.nstims = length(parms.stimnames);
  parms.nextra = length(parms.extra_stimnames); 
  parms.mindur = parms.minfrac*parms.TR;
  if parms.nconds ~= parms.ntypes
    error('condnames and typenames length mismatch');
  end;
  if ~exist(parms.fname,'file')
    error('file %s not found',parms.fname);
  end;
  [fdir,fstem,fext] = fileparts(parms.fname);
  % remove problematic characters
  if isempty(parms.outstem)
    parms.outstem = abcd_clean_fstem(fstem);
  end;
  % calculate time at the end of each TR
  parms.TR_offset = linspace(parms.TR,parms.numTRs * parms.TR,parms.numTRs);
  parms.TR_onset = parms.TR_offset - parms.TR;
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

function [onset,offset,resp_times] = check_offsets(onset,offset,resp_times,eventname,parms)
  ind_empty = find(cellfun(@isempty,onset) | cellfun(@isempty,offset));
  if ~isempty(ind_empty)
    if parms.verbose
      fprintf('%s: WARNING: %s event has %d onsets and %d offsets\n',...
        mfilename,eventname,nnz(~cellfun(@isempty,onset)),nnz(~cellfun(@isempty,offset)));
      % NOTE: this may indicate file truncation
    end
    ind_keep = setdiff([1:length(onset)],ind_empty);
    onset = onset(ind_keep);
    offset = offset(ind_keep);
    resp_times = resp_times(ind_keep);
  end;
  onset = cell2mat(onset);
  offset = cell2mat(offset);
  if iscell(resp_times), resp_times = cell2mat(resp_times); end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_files(eventname,onset,offset,resp_times,ind_start,event_ref_time,parms)
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
    fname_rt = sprintf('%s/%s_scan%d_%s_rt.txt',...
      parms.outdir,parms.outstem,s,eventname);

    if ~exist(fname_1D,'file') || ~exist(fname_txt,'file') || ...
       ~exist(fname_block,'file') || ~exist(fname_fsl,'file') ||...
       ~exist(fname_rt,'file') || parms.forceflag
      if ~isempty(onset(ind_scan))
        % subtract start time and convert to seconds
        rel_onset  = (onset(ind_scan) - event_ref_time(ind_scan))/1000;
        rel_offset = (offset(ind_scan) - event_ref_time(ind_scan))/1000;
        % get resp_times for this scan
        rel_resp_times = resp_times(ind_scan)/1000;
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
        if ~strcmp(eventname, 'cue')
          fid = fopen(fname_block,'wt');  
          if fid<0, error('failed to open %s for writing',fname_block); end;    
          fprintf(fid,'%s\n',sprintf('%0.2f ',rel_onset(1)));
          fclose(fid);
          fid = fopen(fname_block_dur,'wt');  
          if fid<0, error('failed to open %s for writing',fname_block_dur); end;
          block_len = rel_offset(length(rel_offset))-rel_onset(1);
          fprintf(fid,'%s\n',sprintf('%0.2f ',block_len));
          fclose(fid);
        else
          fid = fopen(fname_block,'wt');
          if fid<0, error('failed to open %s for writing',fname_txt); end;
          fprintf(fid,'%s\n',sprintf('%0.2f ',rel_onset));
          fclose(fid);
          fid = fopen(fname_block_dur,'wt');  
          if fid<0, error('failed to open %s for writing',fname_block_dur); end;
          block_len = mean(rel_offset-rel_onset);
          fprintf(fid,'%s\n',sprintf('%0.2f ',block_len));
          fclose(fid);
        end 
        % write to fsl file
        fid = fopen(fname_fsl,'wt');
        if fid<0, error('failed to open %s for writing',fname_fsl); end;
        for i=1:length(rel_onset)
          fprintf(fid, '%0.2f %0.2f 1 \n',rel_onset(i), rel_offset(i)-rel_onset(i));
        end
        fclose(fid);
        % write to rt file
        fid = fopen(fname_rt,'wt');
        if fid<0, error('failed to open %s for writing',fname_rt); end;
        for i=1:length(rel_resp_times)
          fprintf(fid, '%0.3f\n',rel_resp_times(i));
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
        % write empty rt file
        fid = fopen(fname_rt,'wt');
        if fid<0, error('failed to open %s for writing',fname_rt); end;
        fclose(fid);    
      end
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [event_info,start_time,all_types,all_stims,all_targets,all_procs,errcode,errmsg] = get_event_info(parms)
  
  event_info=[]; start_time=[];
  all_types=[]; all_stims = []; all_targets = []; all_procs = [];
  errcode = 0; errmsg = [];
  
  try
    % write event info to file
    fname_csv = abcd_check_eprime_sprdsh(parms.fname, parms.colnames,...
                  parms.fieldnames, parms.outdir, parms.outstem, parms.forceflag, parms.verbose);
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
  if isempty(regexpi(experiment,'nback'))
    if parms.verbose, fprintf('%s: ERROR: wrong experiment name in e-prime file %s: %s\n',...
      mfilename,parms.fname,experiment); end
    errcode = 1;
    errmsg = 'wrong experiment name';
    return;
  else
    if parms.verbose, fprintf('%s: experiment name: %s\n',mfilename,experiment); end
  end
  
  % check for behavioral or GE-specific experiment
  behav_flag = 0; ge_flag = 0;
  if ~isempty(regexpi(experiment,'behavioral'))
    behav_flag = 1;
    expected_start_delay = [0,0];
    expected_init_delay = [0,0];
  elseif ~isempty(regexp(experiment,'_GE'))
    ge_flag = 1;
    expected_start_delay = parms.expected_start_delay_GE*[1,1];
    expected_init_delay = [parms.expected_init_delay_GE_run1,parms.expected_init_delay_GE_run2];
  else
    expected_start_delay = parms.expected_start_delay_Siemens*[1,1];
    expected_init_delay = parms.expected_init_delay_Siemens*[1,1];
  end

  % get trigger times
  trig_time = [];
  if isfield(event_info,'getready_rttime') || isfield(event_info,'getready2_rttime')
    % get initial trigger times
    trig_time = [];
    if isfield(event_info,'getready_rttime')
      ind_trig = find(~cellfun(@isempty,{event_info.getready_rttime}));
      if ~isempty(ind_trig)
        trig_time = event_info(ind_trig(1)).getready_rttime;
      end
    end
    if isfield(event_info,'getready2_rttime')
      ind_trig = find(~cellfun(@isempty,{event_info.getready2_rttime}));
      if ~isempty(ind_trig)
        trig_time = [trig_time,event_info(ind_trig(1)).getready2_rttime];
      end
    end
    % exclude zeros
    trig_time = nonzeros(trig_time)';
  end
  
  % get start index
  all_procs = {event_info.procedure_block};
  ind_start=[];
  start_proc=[];
  for i=1:parms.nprocedures
    proc = parms.procedures{i};
    ind_proc = find(strcmp(proc,all_procs));
    if ~isempty(ind_proc) && ind_proc(end)<length(event_info)
      ind_start(i) = ind_proc(end);
      start_proc{i} = proc;
    end
  end

  % set start time based on trigger
  start_time = [];
  ind_proc = [];
  for i=1:length(ind_start)
    ind_proc(i) = find(strcmp(start_proc{i},parms.procedures));
    if ind_proc(i)==1
      start_time = [start_time,event_info(ind_start(i)).getready_rttime];
    else
      start_time = [start_time,event_info(ind_start(i)).getready2_rttime];
    end
  end
  % exclude zeros
  start_time = nonzeros(start_time)';

  % reset expected delays if only one run
  if length(ind_proc)==1
    expected_start_delay = expected_start_delay(ind_proc);
    expected_init_delay = expected_init_delay(ind_proc);
  end

  if ~behav_flag && ~ge_flag
    % NOTE: for behavioral only, only one trigger is sent/recorded
    %       set start_time to trigger
    % NOTE: for Siemens/Philips, only one trigger is sent/recorded
    %       set start_time relative to trigger with the standard start delay for Siemens/Philips
    start_time = start_time + 1000*parms.start_delay_Siemens;
  end

  % get init_time
  ind_init = ind_start + 1;
  init_time = [event_info(ind_init).cuefix_onset];
  % exclude zeros
  init_time = nonzeros(init_time)';

  % remove extra triggers
  if length(trig_time) > length(init_time)
    trig_time = trig_time(1:length(init_time));
  end  
  % remove extra starts
  if length(start_time) > length(init_time)
    start_time = start_time(1:length(init_time));
  end

  if ge_flag
    % check whether start delay diff is an integer multiple of TR
    %   indicating missed triggers
    if length(start_time)==1
      ex_start_delay = expected_start_delay(1);
    else
      ex_start_delay = expected_start_delay;
    end
    start_delay_diff = (start_time - trig_time) - 1000*ex_start_delay;
    start_delay_diff_TR = start_delay_diff / (1000 * parms.TR);
    start_delay_diff_TR_int = round(start_delay_diff_TR);
    idx_miss = find(ismember(start_delay_diff_TR_int,[1:3]) & abs(start_delay_diff_TR - start_delay_diff_TR_int) <= 0.01);
    if ~isempty(idx_miss)
      if parms.verbose
        for i=1:length(idx_miss)
          fprintf('%s: WARNING: adjusting start_time for run %d by %0.2f seconds due to indications of missed triggers\n',...
            mfilename,idx_miss(i),start_delay_diff(idx_miss(i))/1000);
        end
      end
      % adjust start_time by subtracting start_delay_diff
      start_time(idx_miss) = start_time(idx_miss) - start_delay_diff(idx_miss);
      %% NOTE: start_time represents the start of the non-dummy imaging volume
      %%       init_time is when the stimulus starts
    end
  end

  % check for disparity between trig_time and init_time
  abcd_check_eprime_timing(trig_time,start_time,init_time,...
    expected_start_delay,expected_init_delay,...
    parms.outdir,parms.outstem,...
    parms.max_start_diff,parms.max_init_diff,...
    parms.verbose);

  % get block, stim, and target types for events
  event_info_events = remove_non_events(event_info);
  all_types = {event_info_events.block_type};
  all_stims = {event_info_events.stim_type};
  all_targets = {event_info_events.target_type};

  % check stim_resp
  stim_resp = {event_info_events.stim_resp};
  correct_response = {event_info_events.correct_response};  
  if any(cellfun(@isstr,stim_resp)) || any(cellfun(@iscell,stim_resp))    
    if ~any(cellfun(@isstr,correct_response))
      if parms.verbose, fprintf('%s: ERROR: string stim_resp values without string correct_response values in e-prime file %s\n',...
        mfilename,parms.fname); end
      errcode = 1;
      errmsg = 'string stim_resp values without string correct_response values';
      return;
    end;
    % get response numbers and names for unique, non-empty correct responses
    uniq_correct_response = unique(correct_response(find(~cellfun(@isempty,correct_response))));
    num_cresp = length(uniq_correct_response);
    uniq_correct_response_nums = [];
    uniq_correct_response_names = [];
    for i=1:length(uniq_correct_response)
      k = regexp(uniq_correct_response{i},'(?<num>\d+),{(?<name>\w+)}','names');
      if ~isempty(k)
        uniq_correct_response_nums(i) = str2num(k.num);
        uniq_correct_response_names{i} = k.name;
      else
        % check for single character r, g, b, or y
        k = regexp(uniq_correct_response{i},'^(?<char>[rgby])$','names');
        if ~isempty(k)
          uniq_correct_response_nums(i) = i;
          uniq_correct_response_names{i} = k.char;
        else
          if parms.verbose, fprintf('%s: ERROR: string correct_response with unexpected pattern (%s) in e-prime file %s\n',...
            mfilename,uniq_correct_response{i},parms.fname); end
          errcode = 1;
          errmsg = sprintf('string go_cresp with unexpected pattern (%s)',uniq_correct_response{i});
          return;
        end
      end
    end;
    if parms.verbose, fprintf('%s: WARNING: replacing stim_resp strings with numeric\n', mfilename); end
    for i=1:length(event_info)
      % assign numbers to each stim_resp
      stim_resp = event_info(i).stim_resp;
      if iscell(stim_resp), stim_resp = stim_resp{1}; end;
      if isstr(stim_resp)
        k = find(strcmp(stim_resp,uniq_correct_response_names));
        stim_resp = uniq_correct_response_nums(k);
      end;
      event_info(i).stim_resp = stim_resp;
    end;
  end;

  % check corect_response
  if any(cellfun(@isstr,{event_info.correct_response}))
    if parms.verbose
      fprintf('%s: WARNING: replacing correct_response strings with numeric\n',...
        mfilename);
    end;
    for i=1:length(event_info)
      if isstr(event_info(i).correct_response)
        % remove text description of response
        event_info(i).correct_response = str2num(regexprep(event_info(i).correct_response,',.+',''));
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check for no responses
function [nresp, errcode, errmsg] = check_no_responses(event_info,parms)
  nresp = 0; errcode = 0; errmsg = [];
  % remove non-events
  event_info_events = remove_non_events(event_info);
  % check stim_resp
  stim_resp = {event_info_events.stim_resp};
  % count number of responses
  nresp = nnz(~cellfun(@isempty,stim_resp));
  % if no responses, set errcode and errmsg
  if nresp == 0
    if parms.verbose, fprintf('%s: ERROR: no responses in %s\n',mfilename,parms.fname); end
    errcode = 1;
    errmsg = 'no responses';
    return;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [event_info,switch_flag,errcode,errmsg] = nback_switch(event_info,parms)

  [switch_flag,acc1,errcode,errmsg] = nback_switch_flag(event_info,parms); 
  if errcode, return; end;

  if switch_flag 
    event_info = nback_switch_event(event_info);
    [~,acc2,~] = nback_switch_flag(event_info,parms);
    if acc2 < acc1
      event_info = nback_switch_event(event_info);
      switch_flag = 0; 
    end 
  end
 
  % write to csv
  if switch_flag 
    fname_csv_out  = sprintf('%s/%s_events_switched.csv',...
          parms.outdir,parms.outstem); 
    if ~exist(fname_csv_out,'file') || parms.forceflag
      mmil_struct2csv(event_info,fname_csv_out);
    end; 
  end  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [switch_flag,accuracy,errcode,errmsg] = nback_switch_flag(event_info,parms) 
  errcode = 0; errmsg = [];
  switch_flag = 0;

  % remove non-events
  event_info = remove_non_events(event_info);
  % get responses
  resp = {event_info.stim_resp}; 
  % get correct responses
  correct_resp = {event_info.correct_response};
  correct_total = numel(correct_resp);
  % count number of responses that match correct response
  correct = 0;
  for i=1:correct_total
    if correct_resp{i}==resp{i}
      correct = correct+1;
    end;
  end;
  % calculate accuracy
  accuracy = 100*correct/correct_total; 
  if parms.verbose, fprintf('%s: accuracy = %0.1f%%\n',mfilename,accuracy); end
  if accuracy < 100*parms.switch_thresh
    if parms.verbose, fprintf('%s: accuracy < %0.1f%%, switching button responses\n',...
      mfilename,100*parms.switch_thresh); end
    switch_flag = 1;
  end;
return 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function event_info = nback_switch_event(event_info)
  for i=1:length({event_info.stim_resp})
    if event_info(i).stim_resp == 1
      event_info(i).stim_resp = 2; 
    elseif event_info(i).stim_resp == 2
        event_info(i).stim_resp = 1;
    end
    if (event_info(i).stim_acc >= 0)
      if event_info(i).stim_resp==event_info(i).correct_response
        event_info(i).stim_acc = 1;
      else 
        event_info(i).stim_acc = 0;
      end
    end
  end 
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [event_info,event_info_proc] = remove_non_events(event_info)
  event_info_proc = event_info;
  all_types = {event_info.block_type};
  all_stims = {event_info.stim_type};
  all_targets = {event_info.target_type};
  % NOTE: if block_type is not empty, but stim_type or target_type is empty
  %       that may reflect a truncated file due to premature file transfer
  %       i.e., uploading file before it is fully written
  ind_events = find(~cellfun(@isempty,all_types) &...
                    ~cellfun(@isempty,all_stims) &...
                    ~cellfun(@isempty,all_targets));
  event_info = event_info(ind_events);
return;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [behav,runs_ok,errcode,errmsg] = get_behavioral_data_nback(event_info,parms,start_time) 
  errcode = 0; errmsg = [];
  
  behav = []; 
  behav.('SubjID') = []; behav.('VisitID') = [];
  behav.version = mmil_getfield(event_info(1),'version','UNKNOWN');
  if parms.verbose
    fprintf('%s: experiment version = %s\n',mfilename,behav.version);
  end
  behav.switch_flag = parms.switch_flag;
  behav.perform_flag = 1;

  % check if runs are truncated 
  onset_all = [event_info.stim_onset];
  [ind_start,~] = set_ref(onset_all,start_time);   
  nruns = length(unique(ind_start)); 
  runs = ind_start; 
  runs_ok = []; 
  for i=1:nruns
    run_len = length(find(runs==i));   
    if run_len == parms.expected_numtrials
      runs_ok = [runs_ok i];  
    else
      if parms.verbose
        fprintf('%s: WARNING: run %d is short: %d trials found (expected %d)\n',...
          mfilename,i,run_len,parms.expected_numtrials);
      end
    end 
  end
  nruns = length(runs_ok);  
  if nruns == 1
    ind_ok = find(runs==runs_ok);
    event_info = event_info(ind_ok);
    ind_start = ind_start(ind_ok); 
  elseif nruns == 0
    if parms.verbose, fprintf('%s: ERROR: no valid e-prime runs in %s\n',mfilename,parms.fname); end
    errcode = 1;
    errmsg = 'no valid e-prime runs';
    return;
  end  
  behav.nruns = nruns; 
  
  all_types = {event_info.block_type}; 
  all_targets = {event_info.target_type};
  all_stims = {event_info.stim_type};
  acc = [event_info.stim_acc]; 
  ind_acc = find(acc==1);
  rt = [event_info.stim_rt]; 

  %total
  eventname = sprintf('all_total');
  ind_type = 1:length(all_types); 
  behav.(eventname) = size(ind_type,2);
  for run=1:2
    eventname = sprintf('all_run_%d',run);
    behav.(eventname) = 'NaN';     
  end 
  for run=runs_ok
    eventname = sprintf('all_run_%d',run);
    ind_type_run = [intersect(find(ind_start==run),ind_type)]; 
    behav.(eventname) = size(ind_type_run,2);    
  end 
  
  eventname = sprintf('all_total_correct');
  ind_type_correct = intersect(ind_type,ind_acc);   
  behav.(eventname) = size(ind_type_correct,2);    
  eventname_acc = sprintf('all_total_acc');
  behav.(eventname_acc) = size(ind_type_correct,2)/size(ind_type,2); 
  eventname_rt = sprintf('all_total_correct_rt');
  behav.(eventname_rt) = mean(rt(ind_type_correct)); 
  eventname_std = sprintf('all_total_correct_std');
  behav.(eventname_std) = std(rt(ind_type_correct)); 
  
  
  for i=1:parms.nconds
    type = parms.typenames{i}; 
    cond = parms.condnames{i};  
    
    % 0back and 2back total 
    eventname = sprintf('block_%s_total',cond);
    ind_type = find(strcmp(type,all_types)); 
    behav.(eventname) = size(ind_type,2);
    for run=1:2
      eventname = sprintf('block_%s_run_%d',cond,run);
      behav.(eventname) = 'NaN';    
    end 
    for run=runs_ok
      eventname = sprintf('block_%s_run_%d',cond,run);
      ind_type_run = [intersect(find(ind_start==run),ind_type)]; 
      behav.(eventname) = size(ind_type_run,2);    
    end 
    
    % 0back and 2back total correct ACC and RT 
    eventname = sprintf('block_%s_total_correct',cond);
    ind_type_correct = intersect(ind_type,ind_acc);   
    behav.(eventname) = size(ind_type_correct,2);    
    eventname_acc = sprintf('block_%s_total_acc',cond);
    behav.(eventname_acc) = size(ind_type_correct,2)/size(ind_type,2); 
    eventname_rt = sprintf('block_%s_total_correct_rt',cond);
    behav.(eventname_rt) = mean(rt(ind_type_correct)); 
    eventname_std = sprintf('block_%s_total_correct_std',cond);
    behav.(eventname_std) = std(rt(ind_type_correct)); 
    for run=1:2
      eventname = sprintf('block_%s_run_%d_correct',cond,run);
      behav.(eventname) = 'NaN';     
      eventname_acc = sprintf('block_%s_run_%d_acc',cond,run);
      behav.(eventname_acc) = 'NaN';  
      eventname_rt = sprintf('block_%s_run_%d_correct_rt',cond,run);
      behav.(eventname_rt) = 'NaN';  
      eventname_std = sprintf('block_%s_run_%d_correct_std',cond,run);
      behav.(eventname_std) = 'NaN';  
    end
    for run=runs_ok
      eventname = sprintf('block_%s_run_%d_correct',cond,run);
      ind_type = find(strcmp(type,all_types)); 
      ind_type_correct = intersect(intersect(find(ind_start==run),ind_type),ind_acc);   
      behav.(eventname) = size(ind_type_correct,2);    
      eventname_acc = sprintf('block_%s_run_%d_acc',cond,run);
      behav.(eventname_acc) = size(ind_type_correct,2)/size(intersect(find(ind_start==run),ind_type),2); 
      eventname_rt = sprintf('block_%s_run_%d_correct_rt',cond,run);
      behav.(eventname_rt) = mean(rt(ind_type_correct)); 
      eventname_std = sprintf('block_%s_run_%d_correct_std',cond,run);
      behav.(eventname_std) = std(rt(ind_type_correct)); 
    end
    
    for j=1:parms.nstims
      stim = parms.stimnames{j}; 
      % total each stim 
      eventname = sprintf('%s_total',stim); 
      ind_stim = find(strcmp(stim,lower(all_stims))); %find posface
      if isempty(ind_stim)
        ind_stim = find(strcmp(stim,lower(all_targets)));
      end
      behav.(eventname) = size(ind_stim,2); 
      for run=1:2
        eventname = sprintf('%s_run_%d',stim,run); 
        behav.(eventname) = 'NaN'; 
      end 
      for run=runs_ok
        eventname = sprintf('%s_run_%d',stim,run); 
        behav.(eventname) = size(intersect(find(ind_start==run),ind_stim),2);
      end 
      
      % total each stim and correct, ACC and RT
      eventname = sprintf('%s_correct',stim);
      ind_stim_correct = intersect(ind_stim,ind_acc); 
      behav.(eventname) = size(ind_stim_correct,2); 
      eventname = sprintf('%s_acc',stim); 
      behav.(eventname) = size(ind_stim_correct,2)/size(ind_stim,2); 
      eventname = sprintf('%s_correct_rt',stim); 
      behav.(eventname) = mean(rt(ind_stim_correct)); 
      eventname = sprintf('%s_correct_std',stim); 
      behav.(eventname) = std(rt(ind_stim_correct)); 
      for run=1:2
        eventname = sprintf('%s_run_%d_correct',stim,run);
        behav.(eventname) = 'NaN'; 
        eventname_acc = sprintf('%s_run_%d_acc',stim,run); 
        behav.(eventname_acc) = 'NaN';  
        eventname_rt = sprintf('%s_run_%d_correct_rt',stim,run); 
        behav.(eventname_rt) = 'NaN';  
        eventname_std = sprintf('%s_run_%d_correct_std',stim,run); 
        behav.(eventname_std) = 'NaN';  
      end 
      for run=runs_ok
        eventname = sprintf('%s_run_%d_correct',stim,run);
        ind_stim_correct = intersect(intersect(find(ind_start==run),ind_stim),ind_acc); 
        behav.(eventname) = size(ind_stim_correct,2); 
        eventname_acc = sprintf('%s_run_%d_acc',stim,run); 
        behav.(eventname_acc) = size(ind_stim_correct,2)/size(intersect(find(ind_start==run),ind_stim),2); 
        eventname_rt = sprintf('%s_run_%d_correct_rt',stim,run); 
        behav.(eventname_rt) = mean(rt(ind_stim_correct)); 
        eventname_std = sprintf('%s_run_%d_correct_std',stim,run); 
        behav.(eventname_std) = std(rt(ind_stim_correct)); 
      end 
          
      % total each stim && 0back and 2back
      eventname = sprintf('block_%s_%s_total',cond,stim);
      ind_inter_total = intersect(ind_type,ind_stim); % 0back pos
      behav.(eventname) = size(ind_inter_total,2);
      for run=1:2
        eventname = sprintf('block_%s_%s_run_%d',cond,stim,run); 
        behav.(eventname) = 'NaN';  
      end 
      for run=runs_ok
        eventname = sprintf('block_%s_%s_run_%d',cond,stim,run);
        ind_inter_total_run = intersect(intersect(find(ind_start==run),ind_type),ind_stim);   
        behav.(eventname) = size(ind_inter_total_run,2);   
      end 
      
      % total each stim && 0back and 2back and correct, ACC and RT   
      eventname = sprintf('block_%s_%s_correct',cond,stim);
      ind_inter_correct = intersect(ind_inter_total, ind_acc); 
      behav.(eventname) = size(ind_inter_correct,2);
      eventname_acc = sprintf('block_%s_%s_acc',cond,stim);
      behav.(eventname_acc) = size(ind_inter_correct,2)/size(ind_inter_total,2); 
      eventname_rt = sprintf('block_%s_%s_correct_rt',cond,stim);
      behav.(eventname_rt) = mean(rt(ind_inter_correct)); 
      eventname_std = sprintf('block_%s_%s_correct_std',cond,stim);
      behav.(eventname_std) = std(rt(ind_inter_correct)); 
      for run=1:2
        eventname = sprintf('block_%s_%s_run_%d_correct',cond,stim,run);
        behav.(eventname) = 'NaN'; 
        eventname_acc = sprintf('block_%s_%s_run_%d_acc',cond,stim,run);
        behav.(eventname_acc) = 'NaN';
        eventname_rt = sprintf('block_%s_%s_run_%d_correct_rt',cond,stim,run);
        behav.(eventname_rt) = 'NaN';
        eventname_std = sprintf('block_%s_%s_run_%d_correct_std',cond,stim,run);
        behav.(eventname_std) = 'NaN';
      end
      for run=runs_ok
        eventname = sprintf('block_%s_%s_run_%d_correct',cond,stim,run);
        ind_inter_correct_run = intersect(intersect(find(ind_start==run),ind_inter_total),ind_acc); 
        behav.(eventname) = size(ind_inter_correct_run,2);
        eventname_acc = sprintf('block_%s_%s_run_%d_acc',cond,stim,run);
        behav.(eventname_acc) = size(ind_inter_correct_run,2)/size(intersect(find(ind_start==run),ind_inter_total),2); 
        eventname_rt = sprintf('block_%s_%s_run_%d_correct_rt',cond,stim,run);
        behav.(eventname_rt) = mean(rt(ind_inter_correct_run)); 
        eventname_std = sprintf('block_%s_%s_run_%d_correct_std',cond,stim,run);
        behav.(eventname_std) = std(rt(ind_inter_correct_run));
      end
      
      % target lure nonlure total 
      for k=1:parms.nextra 
        extra = parms.extra_stimnames{k};
        ind_extra = find(strcmp(extra,all_targets));
        ind_extra_types = intersect(ind_extra,ind_inter_total);
        eventname = sprintf('block_%s_%s_%s_total',cond,stim,extra);
        behav.(eventname) = size(ind_extra_types,2);
        eventname = sprintf('block_%s_%s_%s_correct',cond,stim,extra);
        ind_extra_correct = intersect(ind_extra_types,ind_acc);  
        behav.(eventname) = size(ind_extra_correct,2);
        eventname = sprintf('block_%s_%s_%s_acc',cond,stim,extra);
        behav.(eventname) = size(ind_extra_correct,2)/size(ind_extra_types,2);
        eventname = sprintf('block_%s_%s_%s_rt',cond,stim,extra);
        behav.(eventname) = mean(rt(ind_extra_correct));
        eventname = sprintf('block_%s_%s_%s_std',cond,stim,extra);
        behav.(eventname) = std(rt(ind_extra_correct));
        for run=1:2
          eventname = sprintf('block_%s_%s_%s_run_%d_total',cond,stim,extra,run);
          behav.(eventname) = 'NaN'; 
          eventname = sprintf('block_%s_%s_%s_run_%d_correct',cond,stim,extra,run);
          behav.(eventname) = 'NaN'; 
          eventname = sprintf('block_%s_%s_%s_run_%d_acc',cond,stim,extra,run);
          behav.(eventname) = 'NaN';
          eventname = sprintf('block_%s_%s_%s_run_%d_correct_rt',cond,stim,extra,run);
          behav.(eventname) = 'NaN';
          eventname = sprintf('block_%s_%s_%s_run_%d_correct_std',cond,stim,extra,run);
          behav.(eventname) = 'NaN';
        end
        for run=runs_ok
          eventname = sprintf('block_%s_%s_%s_run_%d_total',cond,stim,extra,run);
          ind_extra_run = intersect(find(ind_start==run),ind_extra_types);
          behav.(eventname) = size(ind_extra_run,2);
          ind_extra_correct_run = intersect(intersect(find(ind_start==run),ind_extra_types),ind_acc);   
          eventname = sprintf('block_%s_%s_%s_run_%d_correct',cond,stim,extra,run);
          behav.(eventname) = size(ind_extra_correct_run,2);
          eventname = sprintf('block_%s_%s_%s_run_%d_acc',cond,stim,extra,run);
          behav.(eventname) = size(ind_extra_correct_run,2)/size(ind_extra_run,2);
          eventname = sprintf('block_%s_%s_%s_run_%d_correct_rt',cond,stim,extra,run);
          behav.(eventname) = mean(rt(ind_extra_correct_run));
          eventname = sprintf('block_%s_%s_%s_run_%d_correct_std',cond,stim,extra,run);
          behav.(eventname) = std(rt(ind_extra_correct_run));
        end
      end
    end
  end
  
  % performace
  if (behav.block_2_back_total_acc <= 0.6) || ... 
    (behav.block_0_back_total_acc <= 0.6)
    if parms.verbose, fprintf('%s: block_2_back_total_acc or block_0_back_total_acc <= 0.6 \n',mfilename); end
    behav.perform_flag = 0; 
  end 
  
  % write to csv
  fname_csv_out = sprintf('%s/%s_behavioral.csv',parms.outdir,parms.outstem); 
  if ~exist(fname_csv_out,'file') || parms.forceflag
    mmil_struct2csv(behav,fname_csv_out)
  end;
  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function get_behavioral_data_nback_empty(parms) 
  
  behav = []; 
  behav.('SubjID') = []; behav.('VisitID') = [];
  behav.version = [];
  behav.switch_flag = 0; 
  behav.perform_flag = 0; 
  behav.nruns = 0; 
  
  %total
  eventname = sprintf('all_total');
  behav.(eventname) = 'NaN'; 
  for run=1:2
    eventname = sprintf('all_run_%d',run);
    behav.(eventname) = 'NaN';     
  end 
  eventname = sprintf('all_total_correct');
  behav.(eventname) = 'NaN'; 
  eventname_acc = sprintf('all_total_acc');
  behav.(eventname_acc) = 'NaN';  
  eventname_rt = sprintf('all_total_correct_rt');
  behav.(eventname_rt) = 'NaN'; 
  eventname_std = sprintf('all_total_correct_std');
  behav.(eventname_std) = 'NaN'; 
  
  for i=1:parms.nconds
    type = parms.typenames{i}; 
    cond = parms.condnames{i};  
    
    % 0back and 2back total 
    eventname = sprintf('block_%s_total',cond);
    behav.(eventname) = 'NaN'; 
    for run=1:2
      eventname = sprintf('block_%s_run_%d',cond,run);
      behav.(eventname) = 'NaN'; 
    end 
    
    % 0back and 2back total correct ACC and RT 
    eventname = sprintf('block_%s_total_correct',cond);  
    behav.(eventname) = 'NaN'; 
    eventname_acc = sprintf('block_%s_total_acc',cond);
    behav.(eventname_acc) = 'NaN';  
    eventname_rt = sprintf('block_%s_total_correct_rt',cond);
    behav.(eventname_rt) = 'NaN'; 
    eventname_std = sprintf('block_%s_total_correct_std',cond);
    behav.(eventname_std) = 'NaN'; 
    for run=1:2
      eventname = sprintf('block_%s_run_%d_correct',cond,run); 
      behav.(eventname) = 'NaN';  
      eventname_acc = sprintf('block_%s_run_%d_acc',cond,run);
      behav.(eventname_acc) = 'NaN';  
      eventname_rt = sprintf('block_%s_run_%d_correct_rt',cond,run);
      behav.(eventname_rt) = 'NaN';  
      eventname_std = sprintf('block_%s_run_%d_correct_std',cond,run);
      behav.(eventname_std) = 'NaN'; 
    end
    
    for j=1:parms.nstims
      stim = parms.stimnames{j}; 
      % total each stim 
      eventname = sprintf('%s_total',stim); 
      behav.(eventname) = 'NaN'; 
      for run=1:2
        eventname = sprintf('%s_run_%d',stim,run); 
        behav.(eventname) = 'NaN'; 
      end 
      
      % total each stim and correct, ACC and RT
      eventname = sprintf('%s_correct',stim);
      behav.(eventname) = 'NaN'; 
      eventname = sprintf('%s_acc',stim); 
      behav.(eventname) = 'NaN';  
      eventname = sprintf('%s_correct_rt',stim); 
      behav.(eventname) = 'NaN'; 
      eventname = sprintf('%s_correct_std',stim); 
      behav.(eventname) = 'NaN'; 
      for run=1:2
        eventname = sprintf('%s_run_%d_correct',stim,run);
        behav.(eventname) = 'NaN'; 
        eventname_acc = sprintf('%s_run_%d_acc',stim,run); 
        behav.(eventname_acc) = 'NaN';  
        eventname_rt = sprintf('%s_run_%d_correct_rt',stim,run); 
        behav.(eventname_rt) = 'NaN'; 
        eventname_std = sprintf('%s_run_%d_correct_std',stim,run); 
        behav.(eventname_std) = 'NaN'; 
      end 
      
      % total each stim && 0back and 2back
      eventname = sprintf('block_%s_%s_total',cond,stim);
      behav.(eventname) = 'NaN'; 
      for run=1:2
        eventname = sprintf('block_%s_%s_run_%d',cond,stim,run);
        behav.(eventname) = 'NaN';    
      end 

      % total each stim && 0back and 2back and correct, ACC and RT   
      eventname = sprintf('block_%s_%s_correct',cond,stim);
      behav.(eventname) = 'NaN'; 
      eventname_acc = sprintf('block_%s_%s_acc',cond,stim);
      behav.(eventname_acc) = 'NaN'; 
      eventname_rt = sprintf('block_%s_%s_correct_rt',cond,stim);
      behav.(eventname_rt) = 'NaN';  
      eventname_std = sprintf('block_%s_%s_correct_std',cond,stim);
      behav.(eventname_std) = 'NaN'; 
      for run=1:2
        eventname = sprintf('block_%s_%s_run_%d_correct',cond,stim,run); 
        behav.(eventname) = 'NaN'; 
        eventname_acc = sprintf('block_%s_%s_run_%d_acc',cond,stim,run);
        behav.(eventname_acc) = 'NaN'; 
        eventname_rt = sprintf('block_%s_%s_run_%d_correct_rt',cond,stim,run);
        behav.(eventname_rt) = 'NaN'; 
        eventname_std = sprintf('block_%s_%s_run_%d_correct_std',cond,stim,run);
        behav.(eventname_std) = 'NaN'; 
      end
      
      % target lure nonlure total 
      for k=1:parms.nextra 
      extra = parms.extra_stimnames{k};
      eventname = sprintf('block_%s_%s_%s_total',cond,stim,extra);
      behav.(eventname) = 'NaN'; 
      eventname = sprintf('block_%s_%s_%s_correct',cond,stim,extra);
      behav.(eventname) = 'NaN';
      eventname = sprintf('block_%s_%s_%s_acc',cond,stim,extra);
      behav.(eventname) = 'NaN';
      eventname = sprintf('block_%s_%s_%s_rt',cond,stim,extra);
      behav.(eventname) = 'NaN';
      eventname = sprintf('block_%s_%s_%s_std',cond,stim,extra);
      behav.(eventname) = 'NaN';
      for run=1:2
        eventname = sprintf('block_%s_%s_%s_run_%d_total',cond,stim,extra,run);
        behav.(eventname) = 'NaN'; 
        eventname = sprintf('block_%s_%s_%s_run_%d_correct',cond,stim,extra,run);
        behav.(eventname) = 'NaN'; 
        eventname = sprintf('block_%s_%s_%s_run_%d_acc',cond,stim,extra,run);
        behav.(eventname) = 'NaN';
        eventname = sprintf('block_%s_%s_%s_run_%d_correct_rt',cond,stim,extra,run);
        behav.(eventname) = 'NaN';
        eventname = sprintf('block_%s_%s_%s_run_%d_correct_std',cond,stim,extra,run);
        behav.(eventname) = 'NaN';
      end
     end 
    end
  end
  
  
  % write to csv
  fname_csv_out = sprintf('%s/%s_behavioral.csv',parms.outdir,parms.outstem); 
  if ~exist(fname_csv_out,'file') || parms.forceflag
    mmil_struct2csv(behav,fname_csv_out)
  end;
  
return;
