from __future__ import print_function
import sys

def print(*args, **kwargs):
    sys.stdout.write(' '.join([str(a) for a in args])+(kwargs['end'] if 'end' in kwargs else ''))
    sys.stdout.flush()