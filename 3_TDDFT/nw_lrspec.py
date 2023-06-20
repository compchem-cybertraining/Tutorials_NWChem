#!/usr/bin/python

##
## nw_lrspec
##
## Parses NWChem output deck for linear response TDDFT roots and
## generates corresponding spectrum, with Lorentzian broadened peaks.
##
## All internal frequencies are in eV, but can output to nm or au, too.
##
## Should work with any version of python >= 2.4.3
##
## For online help run "nw_lrspec --help"
##
## Kenneth Lopata (kenneth.lopata@pnl.gov)
##

import sys
import array
import math
import inspect
from optparse import OptionParser


##
##  Global variables.
##
nr = 0                  # number of roots
re = array.array("d")   # energies in eV
rs = array.array("d")   # oscillator strengths
lnum  = 0               # current input line number

emin = 0.0
de = 0.0


##
## Identifiers for "tagging" the tddft output block and for extracting
## quantities.
##
start_tag = "NWChem TDDFT Module"
end_tag = "Excited state energy"

nr_tag = "No. of roots"
nr_delim_l = ":"
nr_delim_r = ""

os_tag = "Dipole Oscillator Strength"
os_delim_l = " "
os_delim_r = ""

## extract energy in eV between ( and e, e.g. "( 0.12345 eV )"
en_delim_l = "("
en_delim_r = "e"



def err_msg_stop():
    """Prints a message to stderr showing line we stopped at and
    exits."""

    global lnum

    c = inspect.stack()
    calledfrom = c[1][3]
    linenum = c[1][2]

    sys.stderr.write(__file__ + ": Error occured in function " + calledfrom + "().\n")
    sys.stderr.write(__file__ + ": Halted at line " + str(lnum) + " of input.\n")
    sys.exit(1)

    

def parse_get_line():
    """Reads in a line from stdout and returns it, then increments
    line number counter."""

    global lnum
    
    try:
        line = sys.stdin.readline()
        lnum = lnum + 1
        return line
    except:
        sys.stderr.write(__file__ + ": Failed to read line " + str(lnum) + ".\n")
        err_msg_stop()

        

def check_eof(line, target, dl, dr):
    """Checks that we didn't hit the end of the file unexpectedly."""
    if line == "":
        sys.stderr.write(__file__ + ': Failed to find "' + target + '" bounded between "' + dl + '" and "' + dr + '".\n')
        err_msg_stop()


        
def parse_jumpto_tag(target):
    """Scans through stdin until it finds [target] string,
    incrementing lnum as it goes. If it cannot find it, print an error
    message and stop."""

    global lnum

    while True:
        line = parse_get_line()
        check_eof(line, target, "", "")

        if line.rfind(target) != -1:
            return

        

def parse_get_val(target, delim_l, delim_r, dtype):
    """Scans through stdin looking for "foo D bar", where "foo" is the
    lable, "bar" is the value, and "D" is some delimiter. Exits if it
    cannot find it."""
    
    global lnum

    while True:
        line = parse_get_line()
        check_eof(line, target, delim_l, delim_r)

        dlloc = line.rfind(delim_l)               # location of left delimiter, if found
        drloc = line.rfind(delim_r)               # location of right delimiter, if found
        tloc = line.rfind(target)                 # location of tag string, if found
        
        if dlloc != -1 and drloc != -1 and tloc != -1:
            rside = line[dlloc+1:drloc].strip()         # from the left delimiter to the right delimiter

            if dtype == "int":
                try:
                    val = int(rside)
                    return val
                except:
                    sys.stderr.write(__file__+': Failed to convert int: "' + str(rside) + '".\n')
                    err_msg_stop()
                    
            elif dtype == "float":
                try:
                    val = float(rside)
                    return val
                except:
                    sys.stderr.write(__file__+': Failed to convert float: "' + str(rside) + '".\n')
                    err_msg_stop()

            elif dtype == "string":
                try:
                    val = str(rside)
                    return val
                except:
                    sys.stderr.write(__file__+': Failed to convert string: "' + str(rside) + '".\n')
                    err_msg_stop()
                    
            else:
                    sys.stderr.write(__file__+': Bad data type: "' + dtype + '".\n')
                    err_msg_stop()
                
    return val


def ev2au(e_ev):
    return (1.0 / 27.2114) * e_ev


def ev2nm(e_ev):
    return 27.2114 * 2.0 * 2.99792 * 2.41888 * math.pi / e_ev


def parse_stdin():
    """Parse NWChem output for for linear response TDDFT energies and
    oscillator strengths."""

    global lnum
    global nr
    global re, rs

    lnum = 0

    ## Jump to start of TDDFT section.
    parse_jumpto_tag(start_tag)

    
    ## Extract nr, the number of roots.
    nr = parse_get_val(nr_tag, nr_delim_l, nr_delim_r, "int")

    
    ##
    ## Extract all of the peaks.  The tag for the ir'th energy is
    ## "RootXXX", where XXX is the root number, e.g., "Root 35" or
    ## "Root241" (not no space when three digit number).
    ##
    for ir in range(1, nr + 1):

        ## Make the energy tag.
        if ir < 10:
            etag = "Root  " + str(ir)
        elif ir < 100:
            etag = "Root " + str(ir)
        else:
            etag = "Root" + str(ir)

        ener = parse_get_val(etag, en_delim_l, en_delim_r, "float")
        ostr = parse_get_val(os_tag, os_delim_l, os_delim_r, "float")

        ## Store the values.
        re.append(ener)
        rs.append(ostr)

    ## Jump to the end of the TDDFT section.
    parse_jumpto_tag(end_tag)

    

def dump_roots():
    """Dump list of roots to stdout (no broadening)."""
    global re, rs

    ## XXX: add header.
    ## XXX: perhaps add option for specifying delimiter?
    
    for i in range(nr):

        if options.units == "eV":
            ee = re[i]
        elif options.units == "au":
            ee = ev2au(re[i])
        elif options.units == "nm":
            ee = ev2nm(re[i])
        else:
            sys.stderr.write(__file__ + ': Bad unit system: "' + options.units + '".\n')
            err_msg_stop()

        ## space after % so that numbers are aligned even if negative
        print "% 10.8e%s%10.8e" %(ee, options.delim, rs[i])


def dump_spectrum_lorentzian():
    """Dump a spectrum of artifically Lorentzian broadened roots."""
    global options


    ##
    ## multiply by 0.5 as FWHM was supplied by gamma in lorenzian is
    ## actually HWHM.
    ##
    gamma = 0.5 * options.width


    ##
    ## L(w; w0, gamma) = gamma/pi * 1 / [(w-w0)^2 + gamma^2]
    ##
    prefac = gamma/math.pi*de


    ## Loop over all energies in spectrum.
    for ie in range(options.npoints):
        ee = emin + ie*de

        
        ## At this energy, add broadened contribution of each root.
        stot = 0.0

        for ii in range(nr):
            xx0 = ee-re[ii]
            stot = stot + prefac * rs[ii] / ( xx0*xx0 + gamma*gamma)   #Lorentzian


        ## Convert this energy to desired units.
        if options.units == "eV" or options.units == "ev":
            eout = ee
        elif options.units == "au":
            eout = ev2au(ee)
        elif options.units == "nm":
            eout= ev2nm(ee)
        else:
            sys.stderr.write(__file__ + ': Bad unit system: "' + options.units + '".\n')
            err_msg_stop()

        print "% 10.8e%s%10.8e" %(eout, options.delim, stot)


def grab_dump_nwchem_header():
    cchar="#"

    print "%s ====================================" %(cchar)
    print "%s NWChem linear response TDDFT parser " %(cchar)
    print "%s ====================================" %(cchar)
    print "%s" %(cchar)

    branch = parse_get_val("date", "=", "", "string")
    print "%s Date NWChem run : %s" %(cchar, branch)
    
    branch = parse_get_val("branch", "=", "", "string")
    print "%s NWChem branch   : %s" %(cchar, branch)

    branch = parse_get_val("input", "=", "", "string")
    print "%s Input deck      : %s" %(cchar, branch)
    

def dump_roots_header():
    
    global options

    cchar="#"
    
    print "%s" %(cchar)
    print "%s Linear response TDDFT roots." %(cchar)
    print "%s" %(cchar)

    if options.units == "nm":
        print "%s Wavelen [nm] %s Oscil. str." %(cchar,options.delim)
    else:
        print "%s  Energy [%s] %s Oscil. str." %(cchar,options.units,options.delim)
        
    print "%s----------------------------------" %(cchar)



def dump_spec_header():
    
    global options

    cchar="#"
    
    print "%s" %(cchar)
    print "%s Linear response TDDFT spectrum." %(cchar)
    print "%s Roots were artificially broadened by %10.5e eV." %(cchar, options.width)
    print "%s" %(cchar)

    if options.units == "nm":
        print "%s Wavelen [nm] %s Oscil. str." %(cchar,options.delim)
    else:
        print "%s  Energy [%s] %s Oscil. str." %(cchar,options.units,options.delim)
        
    print "%s----------------------------------" %(cchar)



def check_options():
    """Check options and replace with sane values if needed."""
    global options
    
    if options.units != "nm" and options.units != "eV" and options.units != "ev" and options.units != "au":
        sys.stderr.write(__file__ + ': Invalid unit type: "'+ options.units + '".\n')
        sys.exit(1)
        
    if options.npoints < 100:
        sys.stderr.write(__file__ + ": Increase number of points to at least 100 (asked for " + str(options.npoints) + ").\n")
        sys.exit(1)

    if options.width < 0.0:
        sys.stderr.write(__file__ + ": Peak width must be positive, supplied " + str(options.width) + ".\n")
        sys.exit(1)


        
def calc_options():
    global emin, de
    
    ## Spectrum will go between lowest root - 15*width to highest + 15*width.
    emin = min(re) - 15.0*options.width

    if emin < options.width:
        emin = options.width

    emax = max(re) + 15.0*options.width
    de = (emax - emin) / options.npoints

    
    ## Use width of at least two grid points 
    options.width = max (options.width, 2*de)

    
    
    
def main():

    global options

    
    ##
    ## Parse command line options.
    ##
    usage = "%prog [options]\n\n" + \
        "Reads NWChem output from stdin, parses for the linear response TDDFT roots,\n" + \
        "and prints absorption spectrum to stdout.  It will broaden peaks using a\n" + \
        "Lorentzian with FWHM of at least two energy/wavelength spacings.\n\n" + \
        "Example:\n" + \
        "\n\tnw_lrspec -b0.3 -p1000 -wnm < water.nwo > spectrum.dat\n\n" + \
        'Create absorption spectrum in nm named "spectrum.dat" from the NWChem\n' + \
        'output file "water.nwo" named spectrum.dat with peaks broadened by 0.3 eV\n' + \
        "and 1000 points in the spectrum.\n"
    parser = OptionParser(usage=usage)
    parser.add_option("-b", "--broad", type="float", dest="width", 
                      help="broaden peaks (FWHM) by WID eV (default 0.1 eV)", metavar="WID")
    parser.add_option("-p", "--points", type="int", dest="npoints",
                      help="create a spectrum with P points (default 500)", metavar="P")
    parser.add_option("-w", "--units", type="string", dest="units",
                      help="units for frequency:  eV (default), au, nm", metavar="UNT")
    parser.add_option("-d", "--delim", type="string", dest="delim",
                      help="use STR as output separator (four spaces default)", metavar="STR")
    parser.add_option("-x", "--extract-only", action="store_false", dest="genspec",
                      help="only extract roots; do not make spectrum")
    parser.add_option("-c", "--clean", action="store_false", dest="header",
                      help="clean output--data only, no header or comments")
    
    parser.set_defaults(width = 0.1, npoints = 500, units="eV", delim="    ", genspec=True, header=True)
    (options, args) = parser.parse_args()
    
    check_options()

    if options.header:
        grab_dump_nwchem_header()

    parse_stdin()
    calc_options()
    
    if options.genspec:

        if options.header:
            dump_spec_header()
            
        dump_spectrum_lorentzian()
        
    else:
        
        if options.header:
            dump_roots_header()
            
        dump_roots()


    
if __name__ == "__main__":
    main()

