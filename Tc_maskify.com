#! /bin/tcsh -f
#
#     Construct a mask that can be used to squash the contribution of the coordinate atoms
#     James Holton  5-29-23
#
set pdbfile = dry.pdb
set mapfile = 2fofc.map
set outfile = mask.map

set B = 50
set Tc_scale = 1

set SG = ""
set CELL = ""
set GRID
set AXIS

set tempfile = /dev/shm/${USER}/Tcmap_$$_
mkdir -p /dev/shm/${USER}

#========================================================================
#    Reading command-line arguments:

foreach arg ( $* )
    set key = `echo $arg | awk -F "=" '{print $1}'`
    set val = `echo $arg | awk -F "=" '{print $2}'`
    if("$arg" =~ *.pdb && "$arg" !~ *=*) then
        if(! -e "$arg") then
            echo "ERROR: $arg does not exist"
            exit 9
        endif
        set pdbfile = "$arg"
        continue
    endif
    if("$arg" =~ *.map && "$arg" !~ *=*) then
        set mapfile = "$arg"
    endif
    if( "$arg" =~ [PpCcIiFfRrHh][1-6]* ) then
        set SG = "$arg"
    endif
    if("$key" == mapfile) set mapfile = $val
    if("$key" == pdbfile) set pdbfile = $val
    if("$key" == outfile) set outfile = $val

    if("$key" == B) set B = $val
    if("$key" == Tc_scale) set Tc_scale = $val

    if("$key" == tempfile) set tempfile = $val
    if("$key" == debug) set DEBUG
end



# get cell dimensions 
set CELL = `awk '/CRYST1/{print $2,$3,$4,$5,$6,$7}' $pdbfile`
if("$SG" == "") then
    set pdbSG = `awk '/^CRYST/{print substr($0,56,12)}' $pdbfile | head -1`
    if("$pdbSG" == "R 32") set pdbSG = "R 3 2"
    if("$pdbSG" == "P 21") set pdbSG = "P 1 21 1"
    if("$pdbSG" == "R 3 2" && $CELL[6] == 120.00) set pdbSG = "H 3 2"
    set SGnum = `awk -v pdbSG="$pdbSG" -F "[\047]" 'pdbSG==$2{print;exit}' ${CLIBD}/symop.lib | awk '{print $1}'`
    set SG = `awk -v pdbSG="$pdbSG" -F "[\047]" 'pdbSG==$2{print;exit}' ${CLIBD}/symop.lib | awk '{print $4}'`
    if("$SG" == R3 && $CELL[6] == 120.00) set SG = H3
endif
if("$SG" == "") then
    echo "WARNING: cannot determine space group, using P1"
    set SG = P1
endif
if(! $?SGnum) set SGnum = `awk -v SG="$SG" 'SG==$4{print $1;exit}' ${CLIBD}/symop.lib`

# now get the grid of the map we should use
if(-e "$mapfile") then
    set GRID = `echo | mapdump mapin $mapfile | awk '/Grid sampling/{print "GRID", $(NF-2), $(NF-1), $NF; exit}'`
    set AXIS = `echo | mapdump mapin $mapfile | awk '/Fast, medium, slow/{print "AXIS",$(NF-2), $(NF-1), $NF}' | head -n 1`
endif

# convert coordinate atoms into an occupancy mask
# the function exp(-log(1-occ)/4.5623*Tc_ff(x,50)) goes from 1% at the core to 98.6% at 3A from each atom
# changing occupancy by -log(1-occ)/4.5623 modulates the central occupancy weight
echo $B |\
cat - $pdbfile |\
awk 'NR==1{B=$1;\
   ABC="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";\
   next} /^CRYST/{print;next} /^ATOM|^HETAT/{\
   occ=substr($0,55,6);xyz=substr($0,31,24);\
   w_min=1-occ;if(w_min<0.01)w_min=0.01\
   logocc=log(w_min)/log(0.01);\
   ++n;\
   c=substr(ABC ABC ABC,int(n/10000+1),1);\
   printf("ATOM  %5d TC   TC  %1s%4d    %24s%6.2f%6.2f%12s\n",n%100000,c,n%10000,xyz,logocc,B,"TC")}' |\
cat >! ${tempfile}sfallme.pdb

sfall xyzin ${tempfile}sfallme.pdb mapout ${tempfile}sfalled.map << EOF
mode atmmap
symm $SGnum
SFSG 1
$GRID
FORMFAC NGAUSS 5
EOF
mapmask mapin ${tempfile}sfalled.map mapout ${tempfile}negme.map << EOF
$AXIS
xyzlim asu
EOF

mapmask mapin ${tempfile}negme.map mapout ${tempfile}expme.map << EOF
scale factor -$Tc_scale
EOF

map_func.com -func exp  ${tempfile}expme.map -outfile $outfile

ls -l $outfile

if($?DEBUG || "$tempfile" == "") exit
rm -f ${tempfile}* >& /dev/null

