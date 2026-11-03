function abcd_check_eprime_timing(trig_time,start_time,init_time,expected_start_delay,expected_init_delay,outdir,outstem,max_start_diff,max_init_diff,verbose)
%function abcd_check_eprime_timing(trig_time,start_time,init_time,expected_start_delay,expected_init_delay,outdir,outstem,max_start_diff,max_init_diff,verbose)
%
% Created:  05/16/2024 by Don Hagler
% Prev Mod: 05/17/2024 by Don Hagler
% Prev Mod: 06/06/2024 by Don Hagler
% Last Mod: 06/25/2024 by Don Hagler
%

if ~mmil_check_nargs(nargin,5), reurn; end;
if ~exist('outdir','var') || isempty(outdir), outdir = pwd; end
if ~exist('outstem','var') || isempty(outstem), outstem = 'eprime'; end
if ~exist('max_start_diff','var') || isempty(max_start_diff), max_start_diff = 0.5; end
if ~exist('max_init_diff','var') || isempty(max_init_diff), max_init_diff = 5.0; end
if ~exist('verbose','var') || isempty(verbose), verbose = true; end

% check for disparity between trig_time, start_time, and init_time
if length(trig_time) ~= length(start_time) || length(trig_time) ~= length(init_time)
  % report mismatch between number of trigger, start, and init events (should not happen)
  if verbose
    fprintf('%s: WARNING: found %d trigger events, %d start events, and %d init events\n',...
      mfilename,length(trig_time),length(start_time),length(init_time));
  end
  fname_timing_err = sprintf('%s/%s_timing.err',outdir,outstem);
  fid = fopen(fname_timing_err,'wt');
  if fid<0, error('failed to open %s for writing',fname_timing_err); end;
  fprintf(fid,'WARNING: found %d trigger events, %d start events, and %d init events\n',...
    mfilename,length(trig_time),length(start_time),length(init_time));
  fclose(fid);
else
  timing_check_info = [];
  start_delay = (start_time - trig_time)/1000;
  init_delay = (init_time - start_time)/1000;
  for r=1:length(trig_time)
    if verbose
      % report delay between initial trigger time and start time
      %   i.e., duration of calibration/"dummy" volumes
      if abs(start_delay(r) - expected_start_delay(r)) > max_start_diff
        fprintf('%s: WARNING: start delay relative to trigger for run %d = %0.3f seconds\n',...
          mfilename,r,start_delay(r));
      else
        fprintf('%s: start delay relative to trigger for run %d = %0.3f seconds\n',...
          mfilename,r,start_delay(r));
      end
      % report delay between start time and init time
      %   i.e., time between start of first non-dummy volume and onset of initial event
      if abs(init_delay(r) - expected_init_delay(r)) > max_init_diff
        fprintf('%s: WARNING: init delay relative to start for run %d = %0.3f seconds\n',...
          mfilename,r,init_delay(r));
      else
        fprintf('%s: init delay relative to start for run %d = %0.3f seconds\n',...
          mfilename,r,init_delay(r));
      end
    end
    % record timing info
    timing_check_info(r).run = r;
    timing_check_info(r).trig_time = trig_time(r)/1000;
    timing_check_info(r).start_time = start_time(r)/1000;
    timing_check_info(r).start_delay = start_delay(r);
    timing_check_info(r).expected_start_delay = expected_start_delay(r);
    timing_check_info(r).start_delay_diff = start_delay(r) - expected_start_delay(r);
    timing_check_info(r).init_time = init_time(r)/1000;
    timing_check_info(r).init_delay = init_delay(r);
    timing_check_info(r).expected_init_delay = expected_init_delay(r);
    timing_check_info(r).init_delay_diff = init_delay(r) - expected_init_delay(r);
  end
  % check if these delays match expectations (e.g., difference < 10 msec)?
  max_start_delay_diff = max(abs(start_delay - expected_start_delay));
  max_init_delay_diff = max(abs(init_delay - expected_init_delay));
  if max_start_delay_diff>max_start_diff || max_init_delay_diff>max_init_diff
    fname_timing_err = sprintf('%s/%s_timing.err',outdir,outstem);
    fid = fopen(fname_timing_err,'wt');
    if fid<0, error('failed to open %s for writing',fname_timing_err); end;
    for r=1:length(init_delay)
      % report start delay
      if abs(start_delay(r) - expected_start_delay(r)) > max_start_diff
        fprintf(fid,'WARNING: start delay relative to trigger for run %d = %0.3f seconds -- expected %0.3f seconds\n',...
          r,start_delay(r),expected_start_delay(r));
      else
        fprintf(fid,'NOTE:    start delay relative to trigger for run %d = %0.3f seconds -- expected %0.3f seconds\n',...
          r,start_delay(r),expected_start_delay(r));
      end
      % report init delay    
      if abs(init_delay(r) - expected_init_delay(r)) > max_init_diff
        fprintf(fid,'WARNING: init delay relative to start for run %d = %0.3f seconds -- expected %0.3f seconds\n',...
          r,init_delay(r),expected_init_delay(r));
      else
        fprintf(fid,'NOTE:    init delay relative to start for run %d = %0.3f seconds -- expected %0.3f seconds\n',...
          r,init_delay(r),expected_init_delay(r));
      end
    end
    if max_start_delay_diff>max_start_diff
      if verbose
        fprintf('%s: WARNING: max start delay difference (%0.3f) greater than allowed (%0.3f)\n',...
          mfilename,max_start_delay_diff,max_start_diff);
      end
      fprintf(fid,'WARNING: max start delay difference (%0.3f) greater than allowed (%0.3f)\n',...
        max_start_delay_diff,max_start_diff);
    end
    if max_init_delay_diff>max_init_diff
      if verbose
        fprintf('%s: WARNING: max init delay difference (%0.3f) greater than allowed (%0.3f)\n',...
          mfilename,max_init_delay_diff,max_init_diff);
      end
      fprintf(fid,'WARNING: max init delay difference (%0.3f) greater than allowed (%0.3f)\n',...
        max_init_delay_diff,max_init_diff);
    end
    fclose(fid);
  end
  
  % write timing check info to csv file
  fname_timing_csv = sprintf('%s/%s_timing.csv',outdir,outstem);
  mmil_struct2csv(timing_check_info,fname_timing_csv);
end

