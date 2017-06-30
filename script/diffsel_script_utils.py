import sys
import random

# String handling functions
def strip(str):
    if str[0]=='#':
        return str[1:]
    else:
        return str.strip()

# Color-related functions
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

def boldred(string):
    return bcolors.FAIL+bcolors.BOLD+string+bcolors.ENDC

def param(myparam):
    return bcolors.OKBLUE+str(myparam)+bcolors.ENDC

def data(myparam):
    return bcolors.CYAN+str(myparam)+bcolors.ENDC

def step(string):
    return bcolors.BOLD+bcolors.HEADER+string+bcolors.ENDC

def boldgreen(string):
    return bcolors.BOLD+bcolors.OKGREEN+string+bcolors.ENDC

def good(string):
    return "-- ("+bcolors.OKGREEN+"Good"+bcolors.ENDC+") "+string

def bad(string):
    return "-- ("+boldred("Bad")+") "+string

def success(string):
    return "-- ["+boldgreen("SUCCESS")+"] "+string

def failure(string):
    return "-- ["+boldred("FAILURE")+"] "+string


# Codon functions
bases = ["A", "C", "G", "T"]

def rand_codon():
    return random.choice(bases)+random.choice(bases)+random.choice(bases)

def selected_codon():
    aa1 = ["AAT", "AAC"]
    return random.choice(aa1)

def mutate(codon, proba=100):
    if random.randint(1,100) <= proba:
        result = list(codon)
        result[random.randint(0,2)] = random.choice(bases)
        print("Decided to mutate codon "+codon+" to "+"".join(result)+" with probability "+str(proba))
        return "".join(result)
    else:
        return codon
