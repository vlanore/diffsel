import sys

def strip(str):
    if str[0]=='#':
        return str[1:]
    else:
        return str.strip()

if sys.stdout.isatty():
    class bcolors:
        HEADER = '\033[95m'
        OKBLUE = '\033[34m'
        OKGREEN = '\033[32m'
        WARNING = '\033[93m'
        CYAN = '\033[33m'
        FAIL = '\033[91m'
        ENDC = '\033[0m'
        BOLD = '\033[1m'
        UNDERLINE = '\033[4m'
else:
    class bcolors:
        HEADER = ''
        OKBLUE = ''
        OKGREEN = ''
        WARNING = ''
        CYAN = ''
        FAIL = ''
        ENDC = ''
        BOLD = ''
        UNDERLINE = ''

def error(string):
    return bcolors.FAIL+bcolors.BOLD+string+bcolors.ENDC

def param(myparam):
    return bcolors.OKBLUE+str(myparam)+bcolors.ENDC

def data(myparam):
    return bcolors.CYAN+str(myparam)+bcolors.ENDC

def step(string):
    return bcolors.BOLD+bcolors.HEADER+string+bcolors.ENDC

def success(string):
    return bcolors.BOLD+bcolors.OKGREEN+string+bcolors.ENDC

