function [eprime_nruns,errcode] = abcd_extract_eprime_nback(fname,varargin)
%function abcd_extract_eprime_nback(fname,[options])
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
%
%  Output: 
%    eprime_nruns: number of valid runs in the e-prime in the file  
%    errcode: [0|1] whether the file was successfully processed  
%
% Created : 01/06/17 by Jose Teruel 
% Prev Mod: 04/03/18 by Dani Cornejo
% Last Mod: 07/25/18 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: based on Nback_scanner.py provided by Eric Feczko from OHSU
%       Using style from abcd_extract_eprime_mid.m provided by Donald Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check arguments 
if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
parms = check_input(fname,varargin); 
errcode = parms.errcode; 

% create output directory
mmil_mkdir(parms.outdir);

try 
% create struct with info for each event
[event_info,event_info_proc,start_time,all_types,all_stims,all_targets,all_procs] = get_event_info(parms); 

% switch buttons if necesary 
[event_info,event_info_proc,parms.switch_flag] = nback_switch(event_info,event_info_proc,parms); 

% get behav data and write it to a csv 
[~,eprime_runs] = get_behavoral_data_nback(event_info,parms,start_time);  
parms.eprime_runs = eprime_runs; 
eprime_nruns = length(eprime_runs);
parms.eprime_nruns = eprime_nruns; 
fprintf('%s: number of e-prime runs = %d \n',mfilename,eprime_nruns);

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
        % get times of stim onset and offset
        onset = [event_info(ind_inter).stim_onset]; 
        offset = [event_info(ind_inter).stim_offset];
        % find most recent start time for each event
        [ind_start,event_ref_time] = set_ref(onset,start_time); 
        % create files for each scan
        write_files(eventname,onset,offset,ind_start,event_ref_time,parms);
      end;
    end;   
  end

  % create timing files for cue trials (pre-block instructions)
  onset = []; offset = [];
  eventname = 'cue';
  for i=1:parms.ncues
    ind_cues = find(strcmp(parms.cues{i},all_procs));
    % onset of fixation block cue
    onset = [onset event_info_proc(ind_cues).cuefix_onset];
    % offset of fixation block cue
    % offset = [offset event_info_proc(ind_cues).cuefix_offset];
    % onset of stimulus block cue
    % onset = [onset event_info_proc(ind_cues).([parms.cuenames{i} '_onset'])];
    % offset of stimulus block cue
    offset = [offset event_info_proc(ind_cues).([parms.cuenames{i} '_offset'])];
  end
  % sort onset and offset times
  onset = sort(onset); offset = sort(offset);  
  % find most recent start time for each event
  [ind_start,event_ref_time] = set_ref(onset,start_time);
  % create files for each scan
  write_files(eventname,onset,offset,ind_start,event_ref_time,parms);

end %timing_files_flag

catch 
  
  % get behav data and write it to a csv
  get_behavioral_data_nback_empty(parms); 
  parms.eprime_nruns=0; eprime_nruns=parms.eprime_nruns; 
  parms.errcode=1; errcode=parms.errcode; 
  fprintf('%s: Empty behavioral file created \n',mfilename)
  
end %try 
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
    'eprime_runs',1:2,[],...
    'eprime_nruns',0,0:2,...
    'errcode',0,0:1,...
    ...
    'colnames',  {'NARGUID','SessionDate','SessionTime','ExperimentName','ExperimentVersion',...
                  'Procedure[Block]','BlockType','StimType','Stim.OnsetTime','Stim.OffsetTime','Stim.ACC',...
                  'Cue2Back.OnsetTime','Cue2Back.OffsetTime','CueTarget.OnsetTime','CueTarget.OffsetTime',...
                  'CueFix.OnsetTime','CueFix.OffsetTime','CueFix.StartTime',...
                  'Fix15sec.OnsetTime', 'Fix15sec.OffsetTime','TargetType','Stim.RT',...
                  'CorrectResponse','Stim.RESP'},[],...
    'fieldnames',{'narguid','date','time','experiment','version',...
                  'procedure_block','block_type','stim_type','stim_onset','stim_offset','stim_acc',...
                  'cue2back_onset','cue2back_offset','cue0back_onset','cue0back_offset',...
                  'cuefix_onset','cuefix_offset','cuefix_start',...
                  'fixation_onset', 'fixation_offset','target_type','stim_rt',...
                  'correct_response','stim_resp'},[],...
    'typenames',{'2-Back','0-Back'},[],...
    'condnames',{'2_back','0_back'},[],...
    'stimnames',{'posface','neutface','negface','place'},[],...
    'extra_stimnames',{'target','lure','nonlure'},[],...
    'cues',{'Cue0BackPROC','Cue2BackPROC'},[]...
    'cuenames',{'cue0back','cue2back'},[],...
    'procedures',{'TRSyncPROC', 'TRSyncPROCR2'},[]...
  });
  
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
  if isempty(parms.outstem)
    %parms.outstem = fstem;
    parms.outstem = regexprep(fstem,'\s.+','');
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

function [event_info,event_info_proc,start_time,all_types,all_stims,all_targets,all_procs] = get_event_info(parms)
  
  event_info=[]; event_info_proc = []; start_time=[];
  all_types=[]; all_stims = []; all_targets = []; all_procs = [];
  
  try 
    % write event info to file
    fname_csv = abcd_check_eprime_sprdsh(parms.fname,parms.colnames,...
      parms.fieldnames,parms.outdir,parms.forceflag);  
    event_info = mmil_csv2struct(fname_csv);
  catch 
    fprintf('%s: ERROR: Failed reading e-prime file (format issues) \n',mfilename);  
    return; 
  end
  
  % get start times
  all_procs = {event_info.procedure_block};
  ind_start=[];
  for i=1:parms.nprocedures
    proc = parms.procedures{i};
    ind_proc = find(strcmp(proc,all_procs));
    if ~isempty(ind_proc)
      ind_start(i) = ind_proc(end) + 1;
    end;
  end;
  start_time = [event_info(ind_start).cuefix_start];  
  
  % remove non-events
  all_types = {event_info.block_type};
  ind_events = find(~cellfun(@isempty,all_types));
  event_info_proc = event_info;
  event_info = event_info(ind_events);
  all_types = {event_info.block_type};
  all_stims = {event_info.stim_type};
  all_targets = {event_info.target_type};

  return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [event_info,event_info_proc,switch_flag] = nback_switch(event_info,event_info_proc,parms)
  
 [switch_flag,acc1,~] = nback_switch_flag(event_info,parms); 

 if switch_flag 
   [event_info, event_info_proc] = nback_switch_event(event_info,event_info_proc);
   [~,acc2,~] = nback_switch_flag(event_info,parms);
   if acc2 < acc1
     [event_info, event_info_proc] = nback_switch_event(event_info,event_info_proc);
      switch_flag = 0; 
   end 
 end
 
 % write to csv
 if switch_flag 
   fname_csv_out  = sprintf('%s/%s_events_switched.csv',...
        parms.outdir,parms.outstem); 
   if ~exist(fname_csv_out,'file') || parms.forceflag
     mmil_struct2csv(event_info_proc,fname_csv_out);
   end; 
 end  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [switch_flag,accuracy,errcode] = nback_switch_flag(event_info,parms) 

  errcode = 0;
  switch_flag = 0;

  correct_resp = {event_info.stim_resp}; 
  correct_total = size(cell2mat({event_info.correct_response}),2); 
  resp = {event_info.correct_response}; 
  correct = 0;  
  for i=1:correct_total
    if correct_resp{i}==resp{i}
      correct = correct+1;
    end;
  end;
  accuracy = 100*correct/correct_total; 
  fprintf('%s: accuracy = %0.1f%%\n',mfilename,accuracy);  

  if accuracy == 0
    error('%s: accuracy equals zero, probably format error \n',mfilename);
    errcode = 1; 
    return; 
  elseif accuracy < 100*parms.switch_thresh
    fprintf('%s: accuracy < %0.1f%%, switching button responses\n',...
      mfilename,100*parms.switch_thresh);
    switch_flag = 1;
  end;
 
return 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [event_info, event_info_proc] = nback_switch_event(event_info,event_info_proc)

 for i=1:length({event_info.stim_resp})
   if event_info(i).stim_resp == 1
    event_info(i).stim_resp = 2;
   elseif event_info(i).stim_resp == 2
     event_info(i).stim_resp = 1;
   end
   if event_info(i).stim_resp==event_info(i).correct_response
     event_info(i).stim_acc = 1;
   else 
     event_info(i).stim_acc = 0;
   end  
 end
 for i=1:length({event_info_proc.stim_resp})
   if event_info_proc(i).stim_resp == 1
     event_info_proc(i).stim_resp = 2; 
   elseif event_info_proc(i).stim_resp == 2
      event_info_proc(i).stim_resp = 1;
   end
   if (event_info_proc(i).stim_acc >= 0)
     if event_info_proc(i).stim_resp==event_info_proc(i).correct_response
       event_info_proc(i).stim_acc = 1;
     else 
       event_info_proc(i).stim_acc = 0;
     end 
   end
 end 
 
return; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [behav, runs_ok] = get_behavoral_data_nback(event_info,parms,start_time) 
  
  behav = []; 
  behav.('SubjID') = []; behav.('VisitID') = []; 
  behav.switch_flag = parms.switch_flag;
  behav.perform_flag = 1;
  
  % check if second runs is truncated 
  onset_all = [event_info.stim_onset];
  [ind_start,~] = set_ref(onset_all,start_time);   
  nruns = length(unique(ind_start)); 
  runs = ind_start; 
  runs_ok = []; 
  for i=1:nruns
    run_len = length(find(runs==i));   
    if run_len == 80 %hardcode the run lenght for now
     runs_ok = [runs_ok i];  
    else 
     fprintf('%s: run %d is short: %d trials found (<80 trials) \n',mfilename,i,run_len);
    end 
  end
  nruns = length(runs_ok);  
  if nruns == 1
   new_info = find(runs==runs_ok);
   event_info = event_info(new_info);
   ind_start = ind_start(new_info); 
  elseif nruns == 0
   fprintf('%s: ERROR: no valid e-prime runs \n',mfilename);
   error('%s: ERROR: no valid e-prime runs \n',mfilename); 
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
    fprintf('%s: block_2_back_total_acc or block_0_back_total_acc <= 0.6 \n',mfilename); 
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
