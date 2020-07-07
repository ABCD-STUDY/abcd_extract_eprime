#!/usr/bin/env python3

import sys, io, os
import datetime, dateutil.parser
import numpy as np

# Import pandas, avoiding warning
import warnings
warnings.simplefilter(action='ignore', category=UserWarning)
import pandas as pd
pd.set_option('display.width', 1024)
pd.set_option('max_colwidth',   256)
pd.set_option('display.max_columns', 30)

import glob, json

from difflib import SequenceMatcher
import copy

# Task keys and corresponding file-name identifiers
tsk_fid_dict = {'MID':   ['MID'],
                'nBack': ['NBACK','WM','REC'],
                'SST':   ['SST'] }

encd_optns_list = ['utf-8','utf-16']

# Acceptable extreme values
rows_n_min      = 10
rows_head_n_max =  5
cols_n_min      = 10
nruns_max       =  2   # Maximumn number of runs accepted within a file

head_vars = ['diagnos', 'dir_found', 'file_found', 'pGUIDmatch', 'contents_ok', 'exper', 'exper_ok', 'fname_exp_match', 'datime_ok', 'exp_t0', 'run', 'run_t0', 'tdiff', 'naming_ok', 'fname', 'msg']

Verbose = False

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# ---------------------------------------------------------------------------------------------------------------
def program_description():
    print()

    print('Read E-Prime files saved as csv or tsv spreadsheets, detecting encoding and format.')
    print('Find runs in file, extract their starting time, and calculate their difference relative to a specified date & time.')
    print('Reads E-Prime files in a directory, and pick the file-run closest to a specified date & time.')
    print('Spreadsheets with practice experiments are considered invalid.')
    print('                        Octavio Ruiz.  2017jun-nov01, 2018feb-jun25, 2019apr-jul23, nov-2020feb14, 2020mar06-apr14')
    print()
    print('Usage:')
    print('  ./eprime_sprdsht_get.py                           Print this help')
    print('  ./eprime_sprdsht_get.py file Summary              Read file, print summary of file encoding and contents')
    print('  ./eprime_sprdsht_get.py file Info                 Read file, print diagnostics and file information')
    print('  ./eprime_sprdsht_get.py file InfoCheckName        Read file, print diagnostics and file information,')
    print('                                                    check that file name matches experiment in spreadsheet')
    print('  ./eprime_sprdsht_get.py file ExportFile outfile   Read file, and write a new file containing the interpreted EPrime data in our standard format')
    print('                                                    (tab-separated, no comments before header line, extension = .txt)')
    print('                                                    outfile should be entered with no extension')
    print()
    print('  ./eprime_sprdsht_get.py dir                                                            Read files in directory and print a table with diagnostics, experiment date, time, and filename')
    print('  ./eprime_sprdsht_get.py dir PickFile                                                   Idem, suppress listing, pick file with best format and diagnostics, print diagnostics and file information')
    print('  ./eprime_sprdsht_get.py dir PickFile "YYYYmmdd HHMMSS"                                 Idem, and calculate difference (in minutes) relative to given time')
    print('  ./eprime_sprdsht_get.py dir PickFile "YYYYmmdd HHMMSS" pGUID                           Idem, order file names by specified subject')
    print('  ./eprime_sprdsht_get.py dir PickFile "YYYYmmdd HHMMSS" pGUID Task Info                 Idem, filter for files containing specified task, pick file with smallest pGUID and time differences')
    print('                                                                                         Suppress table listing, print diagnostics and file information')
    print('  ./eprime_sprdsht_get.py dir PickFile "YYYYmmdd HHMMSS" pGUID Task ExportFile outfile   Idem, write picked file using our standard format')
    print('                                                                                         (tab-separated, no comments before header line, extension = .txt)')
    print('                                                                                         outfile should be entered with no extension')
    print('Arguments:')
    print('  pGUID    Participant short keyname or ""')
    # print('  Task     Either "" or one of', tsk_fid_dict.keys() )
    print('  Task     Either "" or one of:  MID, SST, nBack-WM, nBack-Rec.')
    print('           If given nBack only, the program will return the nBack experiment better matching the other arguments')
    print('  -h       Print variable names above output line')
    print('  -v       Verbose operation')
    print('  Info     Print an output line containing:')
    print('     ',  '  '.join(head_vars) )
    print()
    print('Output:')
    print('  diagnos    Sum of values from the following code:')
    print('    = 0  File not found')
    print('    > 0  File was found and read succesfully. Encoding and format are:')
    print('           1   =>   Encoding = utf-8')
    print('           2   =>   Encoding = utf-16')
    print('           4   =>   Separator = Tab, instead of comma')
    print('           8   =>   Quoted rows')
    print('          16   =>   One or more rows before column-names header')
    print('          64   =>   Experiment in spreadsheet matches file name (returned only if option = "FileNameCheck")')
    print('    < 0  File was found and read, but it is not acceptable:')
    print('     -1..-16   =>   Format as described for diag > 0')
    print('         -32   =>   Unable to recognize experiment in spreadsheet')
    print('         -64   =>   Experiment in spreadsheet does not match file name, or practice experiment (returned only if option = "FileNameCheck")')
    print('        -128   =>   Start_time_info not found or unable to extract')
    print('        -256   =>   Non specified error')
    print()
    print('  dir_found        Directory found')
    print('  file_found       Found at least one file in directory')
    print('  pGUIDmatch       Subject in file name matches the specified pGUID')
    print('  contents_ok      File contains E-Prime data')
    print()
    print('  exper            Experiment found inside the file (short code), one of:  MID, SST, nBack-WM, nBack-Rec')
    print('  exper_ok         File contains the specified task')
    print('  fname_exp_match  Experiment in file matches file name')
    print()
    print('  datime_ok        Date and time can be extracted')
    print('  exp_t0           Experiment starting time, according to file contents')
    print('  tdiff            Time difference relative to provided date & time (minutes)')
    print()
    print('  naming_ok        File name starts with "NDAR_INV"')
    print('  fname            Full file fname')
    print()
    print('Use: "echo $?" to check exit status code')
    print()
    print('Examples:')
    print('  ./eprime_sprdsht_get.py   partic1  PickFile  "20170520 164900"  ""  ""  Info')
    print('  ./eprime_sprdsht_get.py   partic1  PickFile  "20170520 164900"  NDAR_INVJPLWZ1Z0  ""   Info')
    print('  ./eprime_sprdsht_get.py   partic1  PickFile  "20170520 164900"  NDAR_INVJPLWZ1Z0  MID  Info')
    print('  ./eprime_sprdsht_get.py   /space/syn07/1/data/ABCD/DAL_ABCD_QC/aux_incoming/CUB/NDAR_INVJPLZW1Z0/baseline_year_1_arm_1/sst-exported_NDAR_INVJPLZW1Z0_baseline_year_1_arm_1_SessionC')
    print('  ./eprime_sprdsht_get.py  "/space/syn07/1/data/ABCD/DAL_ABCD_QC/aux_incoming/*/*INVB9CDPZUA*/*/*exported*"')
    print('  ./eprime_sprdsht_get.py  "/space/syn05/1/data/ABCD/DAL_ABCD_QC/aux_incoming/*/*INVB9CDPZUA*/*/*exported*"')

    print('  ./eprime_sprdsht_get.py  "/space/syn05/1/data/MMILDB/DAL_ABCD_QC/aux_incoming/*/*INVE51WF5J0*/baseline_year_1_arm_1/*INVE51WF5J0*/"')

    print()
# ---------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------
def command_line_get_variables():
    cmnd_syntx_ok = True
    fname = ''
    optn  = ''
    fname_out = ''
    subj = ''
    ref_time  = ''
    task  = '' 
    optn2 = ''
    verbose = False
    header_show = False

    if 2 <= len(sys.argv) and len(sys.argv) <= 9:
        if '-v' in sys.argv:
            verbose = True
            sys.argv.remove('-v')
        if '-h' in sys.argv:
            header_show = True
            sys.argv.remove('-h')

        fname = sys.argv[1]

        if len(sys.argv) == 2:
            optn = 'ListFiles'

        elif len(sys.argv) >= 3:
            optn = sys.argv[2]

            if optn in ['Summary', 'Info', 'InfoCheckName']:
                pass

            elif optn == 'ExportFile' and len(sys.argv) == 4:
                fname_out = sys.argv[3]

            elif optn == 'PickFile':
                if len(sys.argv) == 3:
                    pass

                elif len(sys.argv) == 4:
                    ref_time = sys.argv[3]

                elif len(sys.argv) == 5:
                    ref_time = sys.argv[3]
                    subj = sys.argv[4]
                    pass

                elif len(sys.argv) == 6:
                    cmnd_syntx_ok = False
                    program_description()

                elif len(sys.argv) >= 7:
                    ref_time = sys.argv[3]
                    subj = sys.argv[4]
                    task = sys.argv[5]

                    if len(sys.argv) == 7 and sys.argv[6] in ['Info', 'InfoCheckName']:
                        optn2 = sys.argv[6]

                    elif len(sys.argv) == 8 and sys.argv[6] == 'ExportFile':
                        optn2     = sys.argv[6]
                        fname_out = sys.argv[7]
                    else:
                        cmnd_syntx_ok = False
                        program_description()
                else:
                    cmnd_syntx_ok = False
                    program_description()
            else:
                cmnd_syntx_ok = False
                program_description()
    else:
        cmnd_syntx_ok = False
        program_description()

    return cmnd_syntx_ok, fname, optn, fname_out, ref_time, subj, task, optn2, verbose, header_show
# ---------------------------------------------------------------------------------------------------------------
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# ---------------------------------------------------------------------------------------------------------------
def File_Diagnostics( EPrime_Info ):
     diagnos  = 0
     diagnos += encd_optns_list.index( EPrime_Info['encoding'] ) + 1  if  EPrime_Info['encoding'] in encd_optns_list  else  0
     diagnos += EPrime_Info['sep_tab']     *  4
     diagnos += EPrime_Info['quoted_rows'] *  8
     diagnos += 16 if (EPrime_Info['hoffs'] > 0) else 0

     EPrime_Info['diagnos'] = diagnos

     return EPrime_Info
# ---------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------
def file_read( fname ):
    # Read file into a text buffer.  Determine encoding (UTF-16 or -8), number of rows,
    # if rows are encased by quotation marks, and what are the field separators
    buf = None
    encoding  = None
    txt_lngth = None
    num_lines = None
    sep_tab   = None
    itab      = None
    quoted_rows = None
    iquot       = None
    line_typic_seps_num = None

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    txt = ''
    ienc = 0
    txt_lines = []

    for i, encd in enumerate(encd_optns_list):
        try:
            ff = open(fname, encoding=encd)
            txt = ff.read()
            ff.close()
            
            txt_lngth = len(txt)
            txt_lines = txt.split('\n')
            num_lines = len(txt_lines)
            if num_lines >= rows_n_min:   # and got to this point:
                encoding = encd
                ienc = i+1
                break
        except:
            # File reading failed => not the right encoding, try another one
            pass

    if not encoding:
        return buf, encoding, txt_lngth, num_lines, sep_tab, itab, quoted_rows, iquot, line_typic_seps_num
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Tab- or comma-separated values?
    if txt.count('\t') > cols_n_min:
        sep = '\t'
        sep_tab = True
        itab = 1
    else:
        sep = ','
        sep_tab = False
        itab = 2
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #            Are rows encased by quotation marks? - If so, remove quotes
    # check a "typical" line, after min. number of initial rows
    i_typ = int( (rows_n_min + num_lines) / 2 )
    line = txt_lines[i_typ]
    if line.startswith('"') and line.endswith('"'):
        quoted_rows = True
        iquot = 1
    else:
        quoted_rows = False
        iquot = 0

    if quoted_rows:
        txt = ''
        for line in txt_lines:
            txt = txt + line.strip('"') + '\n'
        txt_lines = txt.split('\n')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # line_typic_seps_num = txt_lines[i_typ].count( sep )

    # Many times, when an EPrime file is uploaded at the site, its encoding is Microsoft, not Unix
    # When I read those files, there are blank lines between each spreadsheet line.
    # Here I check for this case
    # print()
    # print('i_typ=', i_typ )
    # print( [i  for i in range(i_typ-2, i_typ+3)] )
    # print( [txt_lines[i].count( sep )  for i in range(i_typ-2, i_typ+3)] )
    # print( pd.Series( [txt_lines[i].count( sep )  for i in range(i_typ-2, i_typ+3)] ) )
    # print( pd.Series( [txt_lines[i].count( sep )  for i in range(i_typ-2, i_typ+3)] ).min() )
    # print()
    line_typic_seps_num  = pd.Series( [txt_lines[i].count( sep )  for i in range(i_typ-2, i_typ+3)] ).min()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Store "de-quoted" text as a readable memory object, that can be read into a pd.DataFrame
    buf = io.StringIO( txt )
    buf.flush()
    buf.seek(0)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    return buf, encoding, txt_lngth, num_lines, sep_tab, itab, quoted_rows, iquot, line_typic_seps_num
# ---------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------
def sprdsht_read( sprdsh, sep_tab, num_lines, line_typic_seps_num ):
    # Some files have comments in the first line or lines, instead of the column names.
    # Read buffer using increasing row offsets until starting a set of column names consistent with an EPrime file:
    # first row defining  number of columns > cols_n_min  and  number of columns ~ line_typic_seps_num
    data  = pd.DataFrame()
    hoffs = None
    ok  = True
    msg = ''

    for hoffs in range(0, rows_head_n_max):
        data = pd.DataFrame()
        sprdsh.seek(0)
        try:
            if sep_tab:
                data = pd.read_csv(sprdsh, header=hoffs, sep='\t')
            else:
                data = pd.read_csv(sprdsh, header=hoffs )
        except:
            pass

        if (len(data.columns) > cols_n_min) and (len(data.columns) >= (0.9 * line_typic_seps_num)):
            # If row contains column names, there must be none or very few empty (unnamed) cells
            if ['Unnamed' in s for s in data.columns].count(True) < 5:
                break

    # Many times, when an EPrime file is uploaded at the site, its encoding is Microsoft, not Unix
    # When I read those files, there are blank lines between each spreadsheet line.
    if line_typic_seps_num > 0:

        if (len(data) > 0.8 * num_lines) and (len(data.columns) >= (0.8 * line_typic_seps_num)):
            ok = True
        else:
            ok = False
            msg = "Error: unable to interpret file (1). "

    else:
        if (len(data) > 0.4 * num_lines) and len(data.columns) > cols_n_min:
            ok = True
        else:
            ok = False
            msg = "Error: unable to interpret file (2). "

    return data, hoffs, ok, msg
# ---------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------
def ExperimentCheck( data, fname ):
    exp_in_file = None
    exper = 'Unknown'
    fname_exp_match = 0
    ok  = True
    msg = ''
    exp_diagnos = 0
    
    if data.shape[0] >= rows_n_min and data.shape[1] >= cols_n_min:
        pass
    else:
        ok = False
        msg = 'Unable to interpret file as an EPrime spreadsheet. '
        exp_diagnos = -32
        return exp_in_file, exper, fname_exp_match, ok, msg, exp_diagnos


    # Find experiment reported in spreadsheet
    if 'ExperimentName' in data.columns[0]:   # It can be 'ExperimentName' or '0-ExperimentName' or ...?...
        ExpNameCol = data.columns[0]
        exp_in_file = data[ExpNameCol].unique()[0]

        if Verbose:
            print( 'exp_in_file:', exp_in_file )

        if isinstance(exp_in_file, str):   # Some files have a number instead of an experim.description

            if 'practice' in exp_in_file.lower():
                exper = 'Practice'
                ok = False
                msg = 'Practice experiment. '
                exp_diagnos = -64

            else:
                exp_fnameID_list = []
                for tsk in tsk_fid_dict.keys():
                    if tsk.lower() in exp_in_file.lower():
                        exper = tsk
                        exp_fnameID_list = tsk_fid_dict[tsk]

                        if 'nBack'.lower() in exp_in_file.lower():
                            if 'Nback_Rec'.lower() in exp_in_file.lower():
                                exper = 'nBack_Rec'
                            else:
                                exper = 'nBack_WM'
                        break

                if Verbose:
                    print('exp_fnameID_list =', exp_fnameID_list )

                if len(fname):
                    # Check if experiment matches file name identifier
                    for exp_fn in exp_fnameID_list:
                        if fname.lower().rfind(exp_fn.lower()) > (len(fname)-13):   # (look for task id near the end of file name)
                            fname_exp_match = 1
                            break
                else:
                    pass

                if len(fname):
                    if fname_exp_match:
                        exp_diagnos =  64
                # else:
                #   Not an error; but we can check later fname_exp_match if requested.

        else:
            ok = False
            msg = 'Invalid experiment information in spreadsheet. '
            exp_diagnos = -32

    else:
        ok = False
        msg = 'Unable to recognize experiment in spreadsheet. '
        exp_diagnos = -32
    
    return exp_in_file, exper, fname_exp_match, ok, msg, exp_diagnos
# ---------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------
def extract_ini_time( data, Info, exp_diagnos ):
    # This is the way time_0 was calculated before 2019 revision
    ok = False
    msg = ''

    if 'GetReady.RTTime' in data.columns:
        delay = data.iloc[0]['GetReady.RTTime']
        if Verbose:
            print( data.iloc[[0,1,-2,-1]][['SessionDate', 'SessionTime', 'GetReady.RTTime']] )

    elif 'Wait4Scanner.RTTime' in data.columns:
        delay = data.iloc[0]['Wait4Scanner.RTTime']
        if Verbose:
            print( data.iloc[[0,1,-2,-1]][['SessionDate', 'SessionTime', 'Wait4Scanner.RTTime']] )
    else:
        delay = float('nan')
        if Verbose:
            print( data.iloc[[0,1,-2,-1]][['SessionDate', 'SessionTime']] )
            print('Unable to find GetReady.RTTime or Wait4Scanner.RTTime in data.columns')

    if delay == delay:   # not nan
        try:
            if isinstance( delay, str ):
                delay = int( delay )
            delay = int( delay / 1000 )   # seconds
            Info['delay']  = delay
            Info['exp_t0'] = Info['exp_datime'] + datetime.timedelta(seconds=delay)
            Info['datime_ok'] = 1
            ok = True
        except Exception as err:
            Info['delay']  = float('nan')
            Info['exp_t0'] = '0000-00-00 00:00:00'
            ok = False
            exp_diagnos += -128
            msg = 'Unable to assess starting time: %s. ' % str(err) 
    else:
        Info['delay']  = float('nan')
        Info['exp_t0'] = '0000-00-00 00:00:00'
        ok = False
        exp_diagnos += -128
        msg = 'Unable to extract internal delay. '
        
    return Info, ok, msg, exp_diagnos
# ---------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------
def extract_nruns_delays_and_times( data, exper, exp_datime ):
    # Find number of runs in a file, Locate starting time of each run in file
    nruns = 0
    exp_delays = []
    exp_t0s    = []
    ok = False
    msg = ''

    if exper in ['MID','SST']:
        if exper == 'MID':
            col = 'PrepTime.OnsetTime'
        else:  # exper == 'SST':
            col = 'BeginFix.StartTime'

        if col in data.columns:
            starting_rows = data.dropna( axis=0, subset=[col] )
            exp_delays = starting_rows[col]
            nruns = len(exp_delays)
            # print('nruns =', nruns, '.  exp_delays:')
            # print( exp_delays )
            exp_delays = exp_delays.tolist()
        else:
            msg = 'Unable to find column with starting times per run. '
            return nruns, exp_delays, exp_t0s, ok, msg

    # elif exper == 'nBack':
    elif 'nBack' in exper:
        # Look for variables, signals, and values found in nBack experiment files
        col = 'Procedure[Block]'
        if col in data.columns:
            signal_list = ['TRSyncPROC','TRSyncPROCR2']
            for signal in signal_list:
                if Verbose:
                    print('Locate starting time of this run: value in CueFix.StartTime, in the row that follows the last occurrence of the specified signal')
                ind = data.index[ data[col] == signal ].tolist()
                if Verbose:
                    print('nruns =', nruns, '.  ', signal, 'occurs in rows:', ind )
                if len(ind) > 0:
                    if Verbose:
                        print('Last row: CueFix.StartTime =', data.iloc[ind[-1]  ]['CueFix.StartTime'] )
                        print('Next row: CueFix.StartTime =', data.iloc[ind[-1]+1]['CueFix.StartTime'] )
                    try:
                        delay = data.iloc[ind[-1]+1]['CueFix.StartTime']
                        nruns += 1
                        if Verbose:
                            print('here: nruns =', nruns, '.  delay =', delay, '\n')
                        exp_delays.append( delay )
                    except:
                        msg += 'Experiment terminated before last run completed. '
                        if Verbose:
                            print( msg, '\n')
        else:
            msg = 'Unable to find column with starting times per run. '
            return nruns, exp_delays, exp_t0s, ok, msg
    else:
        msg = 'Unrecognized experiment type. '
        return nruns, exp_delays, exp_t0s, ok, msg

    if nruns:
        exp_t0s = []
        try:
            for j in range(0,len(exp_delays)):
                delay = int( exp_delays[j] )
                delay = int( delay / 1000 )   # ms to seconds
                exp_delays[j] = delay
                exp_t0 = exp_datime + datetime.timedelta( seconds=delay )
                exp_t0s.append( exp_t0 )
            ok = True
        except Exception as err:
            msg = 'Unable to assess starting time: %s. ' % str(err) 
    else:
        msg = 'Unable to find starting times per run. '

    return nruns, exp_delays, exp_t0s, ok, msg
# ---------------------------------------------------------------------------------------------------------------
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# ---------------------------------------------------------------------------------------------------------------
#                                                     Main process

# Some of the variables in this structure are internal, used to interpret the spreadsheet;
# the others (output) are explained in the program description.
Info_Prot = {
    'dir_found':  0,
    'file_found': 0,
    'pGUIDmatch': 0.0,
    'naming_ok':  0,
    'fname':      None,

    'encoding':    None,
    'sep_tab':     None,
    'quoted_rows': None,
    'hoffs':       None,
    'n_rows':      None,
    'n_cols':      None,

    'contents_ok': 0,
    'exp_in_file': None,
    'exper':       None,
    'exper_ok':    0,
    'fname_exp_match': 0,
    'diagnos':       -256,

    'datime_ok':  0,
    'exp_datime': None,
    'exp_t0':     '0000-00-00 00:00:00',
    'delay' :     float('nan'),
    'tdiff' :     float('nan'),

    'nruns':     0,
    'run':       0,
    'exp_times': [],
    'run_t0':    '0000-00-00 00:00:00',

    'fout_ok' : None,
    'fout_msg': '',

    'ok' : None,
    'msg': ''
}

def EPrime_Info_and_Data_get( fname ):
    Info = copy.deepcopy( Info_Prot )
    data = pd.DataFrame()
    exp_diagnos = 0

    # # Does file name contain official pGUID format?
    # Info['naming_ok'] = 0
    # f_nopath_name = os.path.basename( fname )
    # if f_nopath_name.startswith('NDAR_INV') and f_nopath_name.count('_') == 2:
    #     Info['naming_ok'] = 1

    # Does file name complies with ABCD format: official pGUID format + _taskname ?
    Info['naming_ok'] = 0
    f_nopath_name = os.path.basename( fname )
    if f_nopath_name.startswith('NDAR_INV') and f_nopath_name.count('_') == 2:
        pGUID_short = f_nopath_name.replace('NDAR_INV','')
        if len(pGUID_short.split('_')[0]) == 8:
            Info['naming_ok'] = 1

    # ---------------------------------------------------------------------------------------------------------------
    #                                 Read file and determine encoding and format

    sprdsh, encoding, txt_lngth, lines_num, sep_tab, itab, quoted_rows, iquot, line_typic_seps_num  =  file_read( fname )
    
    Info['fname']       = fname
    Info['encoding']    = encoding
    Info['sep_tab']     = sep_tab
    Info['quoted_rows'] = quoted_rows

    if Verbose:
        print('encoding =', encoding, ',   txt_lngth =', txt_lngth, ',   lines_num =', lines_num,
           ',   sep_tab =', sep_tab, ',   quoted_rows =', quoted_rows, ',   line_typic_seps_num =', line_typic_seps_num )

    if encoding:
        Info['dir_found'] = 1
        Info['file_found'] = 1
        Info['fname'] = fname
    else:
        Info['ok']  = False
        Info['msg'] = 'Error: unable to read file'
        Info['diagnos'] = 0
        if Verbose:
            print( Info['msg'], '\n')
        return Info, data
    # ---------------------------------------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------------------------------------
    #                                  Read spreadsheet into a pandas DataFrame

    data, hoffs, ok, msg  =  sprdsht_read( sprdsh, sep_tab, lines_num, line_typic_seps_num )

    Info['hoffs']  = hoffs
    Info['n_rows'] = data.shape[0]
    Info['n_cols'] = data.shape[1]
    Info['ok']  = ok
    Info['msg'] = msg

    if ok:
        # Further checks: discard plain-text files that are not really spreadsheets
        if Info['n_rows'] < rows_n_min or Info['n_cols'] < cols_n_min:
            Info['ok'] = False
            Info['msg'] += 'Invalid num. of rows or cols. '

    if not Info['ok']:
        if Verbose:
            print('sprdsh_ok =', ok )
            print('Info =', Info, '\n' )

        Info = File_Diagnostics( Info )
        Info['diagnos'] = -Info['diagnos']
        return Info, data

    # If arrived here, file contanins an EPrime spreadsheet
    Info['contents_ok'] = 1

    if Verbose:
        print('Spreadsheet:  shape =', data.shape, ',  hoffs =', hoffs )
        print('First two and last two rows:')
        print( data.iloc[[0,1,-2,-1]] )
        print()
    # ---------------------------------------------------------------------------------------------------------------


    # ---------------------------------------------------------------------------------------------------------------
    #             Check that file contains a known task

    exp_in_file, exper, fname_exp_match, ok, msg, exp_diagnos  =  ExperimentCheck( data, fname )

    Info['exp_in_file']     = exp_in_file
    Info['exper']           = exper
    Info['fname_exp_match'] = fname_exp_match
    Info['ok']  = ok
    Info['msg'] = Info['msg'] + msg

    if Verbose:
        print('experim_check_ok =', Info['ok'] )

    if not Info['ok']:
        if Verbose:
            print('Info:')
            print( Info, '\n')

        Info = File_Diagnostics( Info )

        # Add exp_diagnos to File_Diagnostics.
        # If exp_diagnos < 0, make the file-encoding number negative
        if exp_diagnos > 0:
            Info['diagnos'] = exp_diagnos + Info['diagnos']
        else:
            Info['diagnos'] = exp_diagnos - Info['diagnos']
        return Info, data

    if Verbose:
        print('Task reported in file:', exper, '. ', msg)
        print()
    # ---------------------------------------------------------------------------------------------------------------


    # ---------------------------------------------------------------------------------------------------------------
    #       Extract experiment time information and starting time of each run in file

    exp_date_str = data.iloc[0]['SessionDate']
    exp_time_str = data.iloc[0]['SessionTime']
    delay = float('nan')
    try:
        exp_datime = dateutil.parser.parse( exp_date_str + ' ' + exp_time_str )
    except:
        exp_diagnos += -128
        Info['ok']  = False
        Info['msg'] += 'Unable to extract or interpret date or time. '

    if Info['ok']:
        Info['exp_datime'] = exp_datime

        if  Info['exper'] == 'nBack_Rec':
            # File contains an nBack recall experiment, performed outside the scanner,
            # therefore there is no waiting for a sync.pulse from the scanner to begin the experiment
            Info['delay']  = float('nan')
            Info['exp_t0'] = Info['exp_datime']
            Info['datime_ok'] = 1
            ok = True
        else:
            Info, ok, msg, exp_diagnos  =  extract_ini_time( data, Info, exp_diagnos )

            if Info['datime_ok']:

                nruns, exp_delays, exp_t0s, ok, msg  =  extract_nruns_delays_and_times( data, Info['exper'], exp_datime )

                Info['nruns']     = nruns
                Info['exp_times'] = exp_t0s
                # Info['run']    = exp_t0s

        Info['ok']  = ok
        Info['msg'] = Info['msg'] + msg

    else:
        if Verbose:
            print( Info['msg'], '\n')
            print('fname:', fname )
            print('exp_date_str:', exp_date_str )
            print('exp_time_str:', exp_time_str )
            print()
    # ---------------------------------------------------------------------------------------------------------------


    # ---------------------------------------------------------------------------------------------------------------
    #                                          Encode diagnostics

    # Summarize file-encoding in variable Info['diagnos']
    Info = File_Diagnostics( Info )

    if Verbose:
        print('Info:')
        print( Info, '\n')

    # If exp_diagnos < 0, combine the negative summary of the file-encoding number with exp_diagnos
    if exp_diagnos < 0:
        Info['diagnos'] = exp_diagnos - Info['diagnos']
    # ---------------------------------------------------------------------------------------------------------------

    return Info, data    # End of EPrime_Info_and_Data_get
# ---------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------
def ExportFile( data, fname_out ):
    ok = False
    msg = ''

    fxpos = fname_out.rfind('.')
    if fxpos > 0:
        if fxpos > (len(fname_out) - 5):   # The rightmost '.' is near the end of filename, and thus consistent with an extension
            fname_bas = fname_out[0:fxpos]
        else:
            fname_bas = fname_out
    else:
        if fxpos < 0:
            fname_bas = fname_out
        else:
            msg = 'Error: invalid output file name'
            return ok, msg

    if len(fname_bas) > 0:
        fname_out = fname_bas + '.txt'

    if Verbose:
        print( '\nWriting tab-separated file:', fname_out, '...' )

    try:
        data.to_csv( fname_out, index=False, sep='\t' )
        ok = True
    except Exception as err:
        msg = 'Unable to write output file: %s: %s' % (fname_out + str(err))

    if Verbose:
        print( msg )

    return ok, msg
# ---------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------
def List_and_Pick_Files( fname, optn, fname_out, subj, ref_time_arg, task_arg, optn2, header_show, quiet=False ):

    Files = pd.DataFrame()
    EPrime_Info = copy.deepcopy( Info_Prot )
    EPrime_Data = pd.DataFrame()

    f_path = fname
    path_to_search = f_path + '/*'
    f_list = glob.glob( path_to_search )

    # Extract "kernel" of subject ID for comparison with file name
    subj = subj.replace('NDAR_INV','').replace('INV','')
    if Verbose:
        print('subj=', subj)

    if not f_list:
        EPrime_Info['fname'] = fname
        EPrime_Info['ok'] = False
        EPrime_Info['msg'] = 'Error: argument is not a directory. '
        # if Verbose or not optn2:
        if Verbose:
            print('f_list:', f_list )
            print( EPrime_Info['msg'] )

        # return ok, msg, Files, EPrime_Info, EPrime_Data
        return Files, EPrime_Info, EPrime_Data


    for j, fname_full in enumerate( f_list ):

        path  = os.path.dirname( fname_full)
        fname = os.path.basename(fname_full)

        if Verbose:
            print('Reading:', fname_full )

        EPrime_Info, EPrime_Data  =  EPrime_Info_and_Data_get( fname_full )

        if EPrime_Info['fname']:

            # If a subject ID was provided: search for a perfect-match substring in the file name then evaluate general string similarity
            pGUIDmatch = 0.0
            if len(subj)>7:
                subloc = fname.find(subj)
                if subloc >= 0:
                    pGUIDmatch = 10.0
            if len(subj)>3:
                pGUIDmatch += round( SequenceMatcher( a=subj, b=fname ).ratio(), 2 )

            if EPrime_Info['exper'] == 'nBack_Rec' and  EPrime_Info['datime_ok']  and  EPrime_Info['nruns'] == 0:
                if Verbose:
                    print('Experiment is a nBack - recall')

                record = pd.DataFrame( dict({'file': fname,
                                            'encoding':    EPrime_Info['encoding'],
                                            'sep_tab':     EPrime_Info['sep_tab'],
                                            'quoted_rows': EPrime_Info['quoted_rows'],
                                            'hoffs':       EPrime_Info['hoffs'],
                                            'n_rows':      EPrime_Info['n_rows'],
                                            'n_cols':      EPrime_Info['n_cols'],
                                            'contents_ok': EPrime_Info['contents_ok'],
                                            'exp_in_file': EPrime_Info['exp_in_file'],
                                            'exper':       EPrime_Info['exper'],
                                            'fname_exp_match': EPrime_Info['fname_exp_match'],
                                            'exp_datime':  EPrime_Info['exp_datime'],
                                            'delay':       EPrime_Info['delay'],
                                            'exp_t0':      EPrime_Info['exp_t0'],
                                            'nruns':       EPrime_Info['nruns'],
                                            'run':         0,
                                            'run_t0':      EPrime_Info['exp_t0'],
                                            'ok':          EPrime_Info['ok'],
                                            'naming_ok':   EPrime_Info['naming_ok'],
                                            'pGUIDmatch':  pGUIDmatch,
                                            'diagnos':     EPrime_Info['diagnos'],
                                            'msg':         EPrime_Info['msg'],
                                            'path': path }),  index=[j] )
                Files = Files.append( record, ignore_index=True )

            else:
                nruns = min( [EPrime_Info['nruns'], nruns_max] )
                for ri in range( 0, nruns ):
                    record = pd.DataFrame( dict({'file': fname,
                                                'encoding':    EPrime_Info['encoding'],
                                                'sep_tab':     EPrime_Info['sep_tab'],
                                                'quoted_rows': EPrime_Info['quoted_rows'],
                                                'hoffs':       EPrime_Info['hoffs'],
                                                'n_rows':      EPrime_Info['n_rows'],
                                                'n_cols':      EPrime_Info['n_cols'],

                                                'contents_ok': EPrime_Info['contents_ok'],
                                                'exp_in_file': EPrime_Info['exp_in_file'],
                                                'exper':       EPrime_Info['exper'],
                                                'fname_exp_match': EPrime_Info['fname_exp_match'],
                                                'exp_datime':  EPrime_Info['exp_datime'],
                                                'delay':       EPrime_Info['delay'],
                                                'exp_t0':      EPrime_Info['exp_t0'],

                                                'nruns':       EPrime_Info['nruns'],
                                                'run':         ri + 1,
                                                'run_t0':      EPrime_Info['exp_times'][ri],
                                                # 'run_t0':      EPrime_Info['run'][ri],
                                                # 'run_t0':      EPrime_Info['run'][ri]['run_t0'],

                                                'ok':          EPrime_Info['ok'],
                                                'naming_ok':   EPrime_Info['naming_ok'],
                                                'pGUIDmatch':  pGUIDmatch,
                                                'diagnos':     EPrime_Info['diagnos'],

                                                'msg':         EPrime_Info['msg'],
                                                'path': path }),  index=[j] )

                    Files = Files.append( record, ignore_index=True )
        else:
            pass

    EPrime_Info =  copy.deepcopy( Info_Prot )
    EPrime_Info['dir_found'] = 1

    EPrime_Data = pd.DataFrame()

    if len(Files) <= 0:
        # EPrime_Info['fname'] = fname
        EPrime_Info['ok'] = False
        EPrime_Info['msg'] = 'Unable to read or interpret files in directory. '
        # if Verbose or not optn2:
        if Verbose:
            print( EPrime_Info['msg'] )

        return Files, EPrime_Info, EPrime_Data

    else:
        EPrime_Info['file_found'] = 1

    # Reorder columns
    Files = Files[['file', 'encoding', 'sep_tab', 'quoted_rows', 'hoffs', 'n_rows', 'n_cols',
                   'ok', 'contents_ok', 'exp_in_file', 'exper', 'fname_exp_match', 'exp_datime', 'delay', 'exp_t0',
                   'nruns', 'run', 'run_t0',
                   'naming_ok', 'pGUIDmatch', 'diagnos', 'msg', 'path']]

    # Sort by ok, date & time, interpreted or not, prefered format
    Files = Files.sort_values( by=['ok','pGUIDmatch','exp_datime', 'exper', 'diagnos', 'naming_ok', 'fname_exp_match'],
                               ascending=[False, False, True, True, True, False, False] )

    Files = Files.reset_index( drop=True )

    # Keep values of pguid-match and naming_ok of best file's name found at this point, in case dir contains no valid files
    pGUIDmatch = Files.iloc[0]['pGUIDmatch']
    naming_ok   = Files.iloc[0]['naming_ok']

    if (Verbose or optn == 'ListFiles') and not quiet:
        print('Found files:')
        print( Files )

    if optn == 'PickFile':
        # Filter file list: valid files and, if specified, task.
        # Find prefered-encoding file for which exp_t0 is closest to given reference date & time
        # Calculate time difference in minutes and order files according to the absolute time time differences

        Files = Files[ Files['ok'] ]
        # Files = Files[ Files['contents_ok'] == 1 ]

        if len(Files) <= 0:
            EPrime_Info['ok'] = False
            EPrime_Info['msg'] = 'No valid files in directory. '
            EPrime_Info['pGUIDmatch'] = pGUIDmatch
            EPrime_Info['naming_ok']  = naming_ok
            # print_EPrime_info( EPrime_Info, header_show )
            return Files, EPrime_Info, EPrime_Data


        tan = 0   # Lenght of task argument; used below
        if task_arg:
            # Use only the number of characters in task_arg to match exper and task_arg
            tan = len(task_arg)
            Files = Files.loc[ [s[0:tan] == task_arg for s in Files['exper']] ]

        if len(Files) <= 0:
            EPrime_Info['ok'] = False
            EPrime_Info['msg'] = 'No files in directory satisfy request. '
            EPrime_Info['pGUIDmatch'] = pGUIDmatch
            EPrime_Info['naming_ok']  = naming_ok
            # print_EPrime_info( EPrime_Info, header_show )
            return Files, EPrime_Info, EPrime_Data


        if 13 <= len(ref_time_arg) and len(ref_time_arg) <= 16:
            ref_time = datetime.datetime.strptime( ref_time_arg, '%Y%m%d %H%M%S')

            if Verbose or not optn2:
                print()
                print('Reference time:', ref_time_arg, ':  ', ref_time )

            Files['tdiff']  = [ round( (s-ref_time).total_seconds()/60, 2 )  for s in Files['run_t0']]
            Files['tdiffa'] = [ abs(s) for s in Files['tdiff'] ]

            Files = Files.sort_values( by=['tdiffa', 'pGUIDmatch', 'exp_t0', 'exper', 'diagnos', 'naming_ok', 'fname_exp_match'],
                                       ascending=[True, False, True, True, True, False, False] )
            Files = Files.drop('tdiffa', axis=1)
        else:
            Files['tdiff']  = float('nan')
            Files = Files.sort_values( by=['pGUIDmatch', 'exp_t0', 'exper', 'diagnos', 'naming_ok', 'fname_exp_match'], ascending=[False, True, True, True, False, False] )

        if Verbose or not optn2:
            print('Files:')
            print( Files )
            print()
            print('Chosen file:')
            print( Files.iloc[0].to_frame().T )

        # Select top file in list for further processing
        fname = f_path + '/' + Files.iloc[0]['file']
        pGUIDmatch = Files.iloc[0]['pGUIDmatch']
        tdiff      = Files.iloc[0]['tdiff']
        run        = Files.iloc[0]['run']
        run_t0     = Files.iloc[0]['run_t0']

        EPrime_Info, EPrime_Data  =  EPrime_Info_and_Data_get( fname )

        EPrime_Info['dir_found']  = 1
        EPrime_Info['file_found'] = 1
        EPrime_Info['pGUIDmatch'] = pGUIDmatch
        EPrime_Info['naming_ok']  = naming_ok
        EPrime_Info['tdiff']  = tdiff
        EPrime_Info['run']    = run
        EPrime_Info['run_t0'] = run_t0

        if 'Info' in optn2:
            # EPrime_Info['dir_found']  = 1
            # EPrime_Info['file_found'] = 1
            # EPrime_Info['exper_ok']   = 1 if EPrime_Info['exper'] == task_arg else 0

            # print('EPrime_Info:')
            # print( EPrime_Info )

            if task_arg == EPrime_Info['exper'][0:tan]:
                if EPrime_Info['exper'] == 'nBack_Rec':
                    EPrime_Info['exper_ok'] = 1
                else:
                    if 0 < EPrime_Info['nruns'] and EPrime_Info['nruns'] <= nruns_max:
                        EPrime_Info['exper_ok'] = 1
                    else:
                        EPrime_Info['exper_ok'] = 0
            else:
                EPrime_Info['exper_ok'] = 0

            # if (EPrime_Info['exper'] == task_arg) and (0 < EPrime_Info['nruns'] and EPrime_Info['nruns'] <= nruns_max):
            #     EPrime_Info['exper_ok'] = 1
            # else:
            #     EPrime_Info['exper_ok'] = 0

            # print_EPrime_info( EPrime_Info, header_show )

        # else:
        #     print_EPrime_info( EPrime_Info )

        # if optn2 == 'ExportFile':
        #     ok, msg  =  ExportFile( EPrime_Data, fname_out )
        #     EPrime_Info['fout_ok']  = ok
        #     EPrime_Info['fout_msg'] = msg

    # else:
    #     print_EPrime_info( EPrime_Info, header_show )

    return Files, EPrime_Info, EPrime_Data
# ---------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------
def print_EPrime_info( EPrime_Info, header_show ):
    sep = ',  '

    if header_show:
        print( '  '.join(head_vars) )

    hn = len( head_vars )
    for j, var in enumerate( head_vars ):
        if var == 'pGUIDmatch':
            print( 1 if EPrime_Info['pGUIDmatch'] > 10.0 else 0, end='' )
        else:
            print( EPrime_Info[var], end='' )
        if (j+1) < hn:
            print( sep, end='' )
    print()
    return
# ---------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------------------
def eprime_file_find_and_interpret( fname, optn1, fname_out, ref_time_arg='', subj='', task_arg='', optn2='', header_show=False, quiet=False ):
    EPrime_Info = copy.deepcopy( Info_Prot )
    ok = False
    msg = ''

    # Determine if argument fname a directory or a file
    dir_found  = 1 if os.path.isdir(  fname ) else 0
    file_found = 1 if os.path.isfile( fname ) else 0

    if file_found:
        EPrime_Info['dir_found']  = 1
        EPrime_Info['file_found'] = 1

        if optn1 in ['Summary', 'Info', 'InfoCheckName', 'ExportFile']:
            try:
                EPrime_Info, EPrime_Data  =  EPrime_Info_and_Data_get( fname )

                if EPrime_Info['ok']  and  optn1 == 'InfoCheckName':
                    if not EPrime_Info['fname_exp_match']:
                        # Report mismatch by replacing exp_diagnos with its negative
                        if EPrime_Info['diagnos'] > 0:
                            EPrime_Info['diagnos'] = -64 - EPrime_Info['diagnos']
                        else:
                            EPrime_Info['diagnos'] = -64 + EPrime_Info['diagnos']
            except Exception as err:
                msg = 'Unable to read or interpret file: %s. ' % str(err)
                EPrime_Info['msg'] += msg
                if Verbose:
                    print( msg, fname )

            EPrime_Info['dir_found']  = 1
            EPrime_Info['file_found'] = 1

            if 'Info' in optn1:
                pass
            elif optn1 == 'ExportFile':
                if EPrime_Info['diagnos'] > 0:
                    ok, msg  =  ExportFile( EPrime_Data, fname_out )
                else:
                    ok = False
                    msg += 'Nothing to export. '
                EPrime_Info['fout_ok']  = ok
                EPrime_Info['fout_msg'] = msg
            else:
                EPrime_Info['msg'] = 'Error: invalid command'
        else:
            EPrime_Info['msg'] = 'Error: fname points to a file and command is not a valid file command'

    else:   # Path expected to be directory or invalid
        if dir_found:
            EPrime_Info['dir_found']  = 1
            EPrime_Info['file_found'] = 0   # Set below; for example after picking a file from inside a specified subdir

            if optn1 in ['ListFiles', 'PickFile']:

                Files, EPrime_Info, EPrime_Data  =  List_and_Pick_Files( fname, optn1, fname_out, subj, ref_time_arg, task_arg, optn2, header_show )

                if EPrime_Info['ok'] and (optn1 == 'ListFiles'):
                    EPrime_Info['diagnos'] = 0

                if optn2 == 'ExportFile':
                    if EPrime_Info['diagnos'] > 0:
                        ok, msg  =  ExportFile( EPrime_Data, fname_out )
                    else:
                        ok = False
                        msg += 'Nothing to export. '

                    EPrime_Info['fout_ok']  = ok
                    EPrime_Info['fout_msg'] = msg

            else:
                EPrime_Info['msg'] = 'Error: fname points to a directory and command is not a valid dir command'

        else:
            EPrime_Info['dir_found']  = 0
            EPrime_Info['diagnos'] = 0
            EPrime_Info['msg'] = 'Unable to find path: %s' % fname

    if not quiet:
        print_EPrime_info( EPrime_Info, header_show )

    return EPrime_Info, ok, msg
# ----------------------------------------------------------------------------------------------------------------------------------
# =============================================================================================================================================



# =============================================================================================================================================
if __name__ == "__main__":

    cmnd_syntx_ok, fname, optn1, fname_out, ref_time_arg, subj, task_arg, optn2, Verbose, header_show  =  command_line_get_variables()

    if cmnd_syntx_ok:
        EPrime_Info, ok, msg  =  eprime_file_find_and_interpret( fname, optn1, fname_out, ref_time_arg, subj, task_arg, optn2, header_show )
    else:
        sys.exit( 64 )   # 64 => command line usage error
# =============================================================================================================================================
