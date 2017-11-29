from __future__ import print_function
import os, sys

def print(*args, **kwargs):
    msg = ' '.join([str(a) for a in args])
    end = kwargs['end'] if 'end' in kwargs else os.linesep
    sys.stdout.write(msg+end)
    sys.stdout.flush()