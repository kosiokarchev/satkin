import sys

def print(*args, **kwargs):
    sys.stdout.write(' '.join([str(a) for a in args])+(kwargs['end'] if 'end' in kwargs else ''))
    if 'flush' in kwargs:
        sys.stdout.flush()