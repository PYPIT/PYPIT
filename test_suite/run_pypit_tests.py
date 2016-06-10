#!/usr/bin/env python
"""
This needs to be run in the TEST_SUITES directory.
Usually found in a Dropbox.

IMPORTANT: Output will be written to TST_PYPIT env variable
IMPORTANT: Path to PYPIT_CALIBS must be set

See README for more
"""
import sys, os, os.path
import pdb
import subprocess
import warnings 

#sys.path.append(os.getenv('PYPIT')+'/src/')
#import pypit, arload

# Point to all sub-folders 
walk = os.walk('./')
instruments = next(walk)[1]

pwd = os.getcwd()

# Loop on instruments
for instr in instruments:
    #if instr in ['Kast_blue','Kast_red']:  # For testing
    #    continue
    # Setups
    setups = next(os.walk(instr))[1]
    for setup in setups:
        # Look for redux file in PYPIT
        redfile = os.getenv('PYPIT')+'/test_suite/'+instr.lower()+'_'+setup.lower()+'.pypit'
        if not os.path.exists(redfile):
            warnings.warn('No redux file: {:s}'.format(redfile))
            warnings.warn('Not testing..')
            continue
        # Edit data directory 
        with open(redfile, 'r') as infile:
            lines = infile.readlines()
        for kk,iline in enumerate(lines):
            if 'data read' in iline:
                dpth = lines[kk+1]
                i0 = dpth.rfind('/')
                newdpth = ' '+pwd+'/'+instr+'/'+setup+dpth[i0:]
                lines[kk+1] = newdpth
            elif 'reduce useflat' in iline:
                cpth = lines[kk]
                i0 = cpth.rfind('/')
                newcpth = os.getenv('PYPIT_CALIBS')+'/'+cpth[i0:]
                lines[kk] = 'reduce useflat '+newcpth
        # Generate folder as need be
        idir = os.getenv('TST_PYPIT')+'/'+instr
        if not os.path.exists(idir):
            os.makedirs(idir)
        wdir = os.path.join(os.getenv('TST_PYPIT'),instr,setup)
        if not os.path.exists(wdir):
            os.makedirs(wdir)
        # Write to TST_PYPIT
        outfile = wdir+'/'+instr.lower()+'_'+setup.lower()+'.pypit'
        with open(outfile, 'w') as ofile:
            for iline in lines:
                ofile.writelines(iline)
        # Run       
        logfile = wdir+'/'+instr.lower()+'_'+setup.lower()+'.log'
        print('Running pypit on {:s} --- '.format(outfile))
        with open(logfile,'w') as f:
            subprocess.call(['python', os.getenv('PYPIT')+'/pypit/pypit.py', outfile], stderr=f, cwd=wdir)#, shell=True)
        print('Done running pypit on {:s} --- '.format(outfile))
        subprocess.call(['tail', logfile])
        # Need some merit of success..


# cd to Test Suite folder
#os.chdir(os.getenv('TST_PYPIT'))

