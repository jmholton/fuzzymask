#! /bin/tcsh -f
#
#   create a solvent mask by averaging over refmac-derived solvent masks
#   with PDB file coordinates randomly shifted using B and occupancy
#
#   put any extra refmac keywords into a file "refmac_opts.txt" before running
#
#   run with no arguments for online help
#
#   James Holton 12-18-18
#
#========================================================================
#    Setting defaults :

# bulk solvent density
set ksol = 0.35824

# solvent mask options
set vdwprobe = 1.4
set ionprobe = 0.8
set rshrink = 0.5

# scale for size of kicks
set kick_scale    = 1.0
# applied to non-water afer above scale
set drykick_scale = 1.0
# number of sub-kicks to perform
set kick_steps = 1
# weight given to thru-bond partners in post-jiggle averaging
set frac_thrubond = 0.9
# number of times to do thru-bond averaging
set ncyc_thrubond = 500
# weight given to shift magnitude re-scaling in post-jiggle optimization
set frac_magnforce = 1.1
# number of times to do magnitude re-scaling in post-jiggle optimization
set ncyc_magnforce = 500

# number of refmac cycles to refine in each sub-kick step
set ncyc_ideal = 5
set ncyc_data = 0

# number of CPUs to use
set seeds = 500
set user_CPUs  = auto

# grid spacing to use for map, default is to let refmac choose because refmac wont listen
set GRID = ""

# output map files
set mapout = "fuzzysolvent.map"
set coordmap = "fuzzycoord.map"
set carbonmap = "carbonized.map"

# output mtz file
set mtzout  = "Fparted.mtz"

# output pdb file
set multiconfpdb = "multiconf.pdb"

if($#argv != 0) goto Setup
# come back at after_Setup:

#========================================================================
#    help message if run with no arguments
Help:
cat << EOF
usage: $0 refme.pdb refme.mtz [ligand.cif] \
     [ksol=${ksol}] \
     [vdwprobe=$vdwprobe] [ionprobe=$ionprobe] [rshrink=$rshrink] \
     [drykick_scale=$drykick_scale] [kick_steps=$kick_steps]\
     [frac_thrubond=${frac_thrubond}] [ncyc_thrubond=${ncyc_thrubond}]\
     [frac_magnforce=${frac_magnforce}] [ncyc_magnforce=${ncyc_magnforce}]\
     [ncyc_ideal=${ncyc_ideal}] [ncyc_data=${ncyc_data}]\
     [seeds=${seeds}] [CPUs=auto] \
     [mapout=${mapout}] [coordmap=${coordmap}]\
     [output=${mtzout}]

where:
refme.pdb      PDB file being refined
refme.mtz      data being used in refinement
ligand.cif     is the geometry info for any ligands (must be one file)

ksol           electron density of bulk solvent far from any coordinate atoms in electrons/A^3
vdwprobe       van der Walls radius for solvent probe in refmac 
ionprobe       vdwprobe for ionic species in refmac 
rshrink        mask shrinkage radius in refmac 
kick_scale     scale-up rms jiggle, to compensate for minimization
drykick_scale  scale-up rms jiggle of non-water atoms, to compensate for minimization

frac_thrubond  weight of thru-bond positional averaging after random atomic jiggle (default: $frac_thrubond)
ncyc_thrubond  number of thru-bond positional averaging cycles to perform after jiggling (default: $ncyc_thrubond)
frac_magnforce scale for shift-magnitude re-normalization after thru-bond averaging cycle (default: $frac_magnforce)
ncyc_magnforce number of shift magnitude re-scaling cycles to perform (default: $ncyc_magnforce)
ncyc_ideal     number of geometry minimization refmac cycles to perform after jiggling (default: $ncyc_ideal)
ncyc_data      number of refmac cycles vs x-ray data to perform before taking solvent mask (default: $ncyc_data)

seeds          number of random maps to average (default: $seeds)
CPUs           number of CPUs to use in parallel (default: $user_CPUs)
$mapout map file of the fuzzy solvent to output in CCP4 format
$coordmap map file of the average fuzzy coordinate atoms
$mtzout    partial structure factors for solvent added to input mtz file

EOF

exit 9


# read command line, set defaults, deploy scripts
after_Setup:

cat << EOF
generating $seeds random masks with:
vdwprobe=$vdwprobe ionprobe=$ionprobe rshrink=$rshrink
kick_scale=$kick_scale
drykick_scale=$drykick_scale
kick_steps=$kick_steps
frac_thrubond=${frac_thrubond} ncyc_thrubond=${ncyc_thrubond}
frac_magnforce=${frac_magnforce} ncyc_magnforce=${ncyc_magnforce}
ncyc_ideal=${ncyc_ideal} ncyc_data=${ncyc_data}
ksol=${ksol}
mapout=${mapout} coordmap=${coordmap}
mtzout=${mtzout}
EOF

# make sure we have a temporary file location all cluster nodes and this node can see
set sharedtmp = ${tempfile}
if("$CLUSTER" != "") then
    # we are using a cluster make sure we can access results
    set sharedtmp = "./"
    if(-w "/scrapp/") mkdir -p /scrapp/${USER}/
    if(-w "/scrapp2/") mkdir -p /scrapp2/${USER}/
    if(-w "/scratch/") mkdir -p /scratch/${USER}/
    if(-w "/global/scratch}/") mkdir -p /global/scratch/${USER}/
    if(-w "/scrapp/${USER}/") set sharedtmp = "/scrapp/${USER}/fuzztemp_`hostname`_$$"
    if(-w "/scrapp2/${USER}/") set sharedtmp = "/scrapp2/${USER}/fuzztemp_`hostname`_$$"
    if(-w "/global/scratch/${USER}/") set sharedtmp = "/global/scratch/${USER}/fuzztemp_`hostname`_$$"
endif

set kick_step_list = `seq 1 $kick_steps`
#set kick_scale_step = `echo $kick_scale $kick_steps | awk '{print sqrt($1*$1/$2)}'`
#set drykick_scale_step = `echo $drykick_scale $kick_steps | awk '{print sqrt($1*$1/$2)}'`
set kick_scale_step = `echo $kick_scale $kick_steps | awk '{print $1/$2}'`
#set drykick_scale_step = `echo $drykick_scale $kick_steps | awk '{print $1/$2}'`
set drykick_scale_step = $drykick_scale

set refmacgrid = ""
if("$GRID" != "") then
    # refmac will not accept this!!!  How do we do it?
    set refmacgrid = "GRID $GRID"
endif

# create the refmac script we will run on each cpu
cat << EOF-script >! refmac_cpu.com
#! /bin/tcsh -f
#  SGE instructions
#\$ -S /bin/tcsh                    #-- the shell for the job
#\$ -o $pwd                         #-- output directory 
#\$ -e $pwd                         #-- error directory 
#\$ -cwd                            #-- job should start in your working directory
#\$ -r y                            #-- if a job crashes, restart
#\$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#\$ -l mem_free=16G                 #-- submits on nodes with enough free memory (required)
#\$ -l arch=linux-x64               #-- SGE resources (CPU type)
#\$ -l netapp=16G,scratch=16G       #-- SGE resources (home and scratch disks)
#\$ -l h_rt=00:10:00                #-- runtime limit (see above; this requests 24 hours)
#\$ -t 1-$seeds                     #-- the number of tasks
# SLURM instructions
#SBATCH 
#SBATCH --array=1-$seeds
#SBATCH 


set seed = \$1

if(\$?SGE_TASK_ID) then
    qstat -j \$JOB_ID
    set seed = \$SGE_TASK_ID
endif
if(\$?SLURM_ARRAY_TASK_ID) then
    set seed = \$SLURM_ARRAY_TASK_ID
endif

if(! \$?CCP4) then
    source ${CBIN}/ccp4.setup-csh
endif
set path = ( . `dirname $0` \$path )
set tempfile = ${tempfile}\$\$
mkdir -p `dirname $tempfile`

# prevent history substitution bombs
history -c
set savehist = ""
set histlit

if(\$?SGE_TASK_ID || \$?SLURM_ARRAY_TASK_ID) then
    hostname
    free -g
    echo "pid= \$\$"
    ps -fH
    df
endif
echo "seed = \$seed"

cp $pdbfile \${tempfile}seed\${seed}kickme.pdb
foreach itr ( $kick_step_list )

# mess it up
cat \${tempfile}seed\${seed}kickme.pdb |\
jigglepdb.awk -v shift=byB -v seed=\$seed \
  -v frac_thrubond=$frac_thrubond -v ncyc_thrubond=$ncyc_thrubond \
  -v frac_magnforce=$frac_magnforce -v ncyc_magnforce=$ncyc_magnforce \
  -v shift_scale=$kick_scale_step \
  -v dry_shift_scale=$drykick_scale_step |\
awk '! /^ATOM|^HETAT/{print;next} \
  substr(\$0,55)+0>0{print substr(\$0,1,16)" "substr(\$0,18)}' >! \${tempfile}seed\${seed}.pdb

# clean it up
if($ncyc_ideal != 0) then
    refmac5 xyzin \${tempfile}seed\${seed}.pdb $LIBSTUFF \
            xyzout \${tempfile}seed\${seed}minimized.pdb << EOF-refmac
#vdwrestraints 0
$otheropts
refi type ideal
ncyc $ncyc_ideal
EOF-refmac
else
    cp \${tempfile}seed\${seed}.pdb \${tempfile}seed\${seed}minimized.pdb
endif

# forget it all if it didnt minimize
if(! -e \${tempfile}seed\${seed}minimized.pdb) then
    echo "BAD" >! ${sharedtmp}seed\${seed}done.txt
    exit 9
endif

# map it out
refmac5 xyzin \${tempfile}seed\${seed}minimized.pdb \
        xyzout \${tempfile}seed\${seed}out.pdb \
        hklin $mtzfile  $LIBSTUFF \
    mskout \${tempfile}mask_seed\${seed}out.map \
    hklout \${tempfile}seed\${seed}out.mtz << EOF-refmac
#vdwrestraints 0
$otheropts
$refmacgrid
ncyc $ncyc_data
solvent vdwprobe $vdwprobe ionprobe $ionprobe rshrink $rshrink
solvent yes
EOF-refmac
if(\$status) then
    echo "what the frak! "
endif

if(! -s \${tempfile}seed\${seed}out.pdb) then
    echo "BAD" >! ${sharedtmp}seed\${seed}done.txt
    exit 9
endif


# get ready for next round
cp \${tempfile}seed\${seed}minimized.pdb \${tempfile}seed\${seed}kickme.pdb

end


if(\$?SGE_TASK_ID) then
    dmesg
    free
    qstat -j \$JOB_ID
endif

# minimize map size to speed up summation
mapmask mapin \${tempfile}mask_seed\${seed}out.map mapout ${sharedtmp}mask_seed\${seed}.map << EOF
xyzlim asu
EOF

if("$GRID" != "") set GRID = ( $GRID )
    set GRID = \`echo | mapdump mapin ${sharedtmp}mask_seed\${seed}.map | awk '/Grid sampling/{print \$(NF-2), \$(NF-1), \$NF; exit}'\`
    echo "using map grid: \$GRID"
endif

#cp \${tempfile}seed\${seed}minimized.pdb ${pwd}/seed\${seed}_minimized.pdb

if("$multiconfpdb" != "") then
    cp \${tempfile}seed\${seed}out.pdb ${sharedtmp}seed\${seed}out.pdb
endif

if("$coordmap" != "") then

if(0) then
# calcualte coordinate map from refmac Fs
fft hklin \${tempfile}seed\${seed}out.mtz \
   mapout \${tempfile}seed\${seed}.map << EOF
LABIN F1=FC PHI=PHIC
GRID \$GRID
EOF
endif

# covert coordinates to a map
cat \${tempfile}seed\${seed}out.pdb |\
awk '! /^ATOM|HETAT/{print;next} {l=\$0;Ee=substr(l,13,2);gsub(Ee," ","")}\
      Ee !~ /^[CNOHS]\$/{Ee=toupper(substr(\$NF,1,2))} {\
      printf("%s%2s%s\n",substr(l,1,12),Ee,substr(l,15,57))}' |\
cat >! \${tempfile}seed\${seed}sfallme.pdb

sfall xyzin \${tempfile}seed\${seed}sfallme.pdb \
   mapout \${tempfile}seed\${seed}.map << EOF
mode atmmap
SYMM $SG
cell $CELL
GRID \$GRID
SFSG 1
EOF

# minimize map size to speed up summation
mapmask mapin \${tempfile}seed\${seed}.map \
  mapout ${sharedtmp}coord_seed\${seed}.map << EOF
xyzlim asu
AXIS X Y Z
EOF

# carbonize the fuzzy coordiantes to check they converge to original model, if all non-hydrogen atoms were carbon
cat \${tempfile}seed\${seed}out.pdb |\
awk '! /^ATOM|HETAT/{print;next} \$NF!="H"{l=\$0\
      printf("%s C%s 5.00           C\n",substr(l,1,12),substr(l,15,47),substr(l,1,61))}' |\
cat >! \${tempfile}seed\${seed}carbonized.pdb

sfall xyzin \${tempfile}seed\${seed}carbonized.pdb \
   mapout \${tempfile}seed\${seed}.map << EOF
mode atmmap
SYMM $SG
cell $CELL
GRID \$GRID
SFSG 1
EOF
# minimize map size to speed up summation
mapmask mapin \${tempfile}seed\${seed}.map \
  mapout ${sharedtmp}carbonized_seed\${seed}.map << EOF
xyzlim asu
AXIS X Y Z
EOF

endif

# clean up
if(! $debug && -e ${sharedtmp}mask_seed\${seed}.map && "$CLUSTER" != "") then
    echo "CLEANING UP CLUSTER NODE "
    rm -f \${tempfile}seed\${seed}.pdb  >& /dev/null
    rm -f \${tempfile}seed\${seed}.map  >& /dev/null
    rm -f \${tempfile}seed\${seed}minimized.pdb  >& /dev/null
    rm -f \${tempfile}seed\${seed}carbonized.pdb  >& /dev/null
    rm -f \${tempfile}seed\${seed}out.pdb >& /dev/null
    rm -f \${tempfile}seed\${seed}out.map >& /dev/null
    rm -f \${tempfile}seed\${seed}out.mtz  >& /dev/null
endif

# signal that map is ready
touch ${sharedtmp}seed\${seed}done.txt

# did we clean up after ourselves
#ls -l ${tempfile}*

EOF-script
chmod a+x refmac_cpu.com


# launch refmac jobs in parallel
# launch jobs on a torq cluster, like pxproc
if("$CLUSTER" == "TORQ") then
    # use a TORQ queue
    if(! $quiet) echo "submitting jobs..."
    rm -f qsubs.log
    foreach seed ( `seq 1 $seeds` )
        echo -n "$seed " >> qsubs.log
        qsub -e seed${seed}_errors.log -o seed${seed}.log -d $pwd ./refmac_cpu.com -F "$seed" >> qsubs.log
    #    sleep 0.1
    end
    goto sum_maps
endif

# launch jobs on a SGE cluster, like the UCSF one
if("$CLUSTER" == "SGE") then
    # use a SGE queue
    if(! $quiet) echo "submitting SGE jobs..."
    qsub -cwd ./refmac_cpu.com 
    #    sleep 0.1
    goto sum_maps
endif

set jobs = 0
set lastjobs = 0
foreach seed ( `seq 1 $seeds` )
    if($quiet) then
        ( ./refmac_cpu.com $seed >! ${sharedtmp}seed${seed}.log & ) >& /dev/null
    else
        ./refmac_cpu.com $seed >! ${sharedtmp}seed${seed}.log &
    endif

    # now make sure we dont overload the box
    @ jobs = ( $jobs + 1 )
    while ( $jobs >= $CPUs )
        sleep 4
        set jobs = `ps -fu $USER | grep refmac_cpu.com | grep -v grep | wc -l`
        if($jobs != $lastjobs && ! $quiet) echo "$jobs jobs running..."
        set lastjobs = $jobs
    end
end
#wait
goto sum_maps



sum_maps:
set seed = 1
if(! $quiet) then
    echo "all jobs launched."
    if (! -e ${sharedtmp}seed${seed}done.txt) echo "waiting for ${seed} to finish..."
endif
while (! -e ${sharedtmp}seed${seed}done.txt)
   sleep 1
   ls -l ${sharedtmp}seed${seed}done.txt >& /dev/null
end
while (! -s ${sharedtmp}mask_seed${seed}.map)
   sleep 2
   ls -l ${sharedtmp}mask_seed${seed}.map >& /dev/null
end

# now start adding up the maps...
echo -n "" >! ${tempfile}rmsds.log
rm -f ${tempfile}sum.map
if(! $quiet) echo -n "summing maps: "
foreach seed ( `seq 1 $seeds` )

    # big rigamarol to make sure files are not stuck in NFS
    set timeleft = 100
    if(! $quiet) echo -n "$seed "
    while(! -e ${sharedtmp}seed${seed}done.txt && $timeleft)
        sleep 3
        @ timeleft = ( $timeleft - 1 )
        ls -l ${sharedtmp}seed${seed}done.txt >& /dev/null
        wait
    end
    set timeleft = 100
    while(! -s ${sharedtmp}mask_seed${seed}.map && $timeleft)
        echo "WARNING: map $seed should be here by now..."
        sleep 5
        @ timeleft = ( $timeleft - 1 )
        ls -l ${sharedtmp}mask_seed${seed}.map >& /dev/null
        wait
    end
    if(! $timeleft) then
        echo ""
        echo "trying to do $seed again..."
        ./refmac_cpu.com $seed >$! ${sharedtmp}seed${seed}.log
        if(-s ${sharedtmp}mask_seed${seed}.map) set timeleft = 1
    endif
    if(! $timeleft) then
        cat ${sharedtmp}seed${seed}.log
        cat refmac_cpu.com.*.${seed}
        set BAD = "timed out waiting for jobs to finish."
        goto exit
    endif
    # file exists and is readable


    # now actually add the maps
    if(! -e ${tempfile}sum.map) then
        cp -p ${sharedtmp}mask_seed${seed}.map ${tempfile}sum.map
    else
        echo maps add |\
        mapmask mapin1 ${sharedtmp}mask_seed${seed}.map mapin2 ${tempfile}sum.map \
           mapout ${tempfile}new.map >> $logfile
        if(! -e ${tempfile}new.map) then
            ls -l ${sharedtmp}mask_seed${seed}.map
            set BAD = "failed to make seed $seed "
            goto exit
        endif
        mv ${tempfile}new.map ${tempfile}sum.map
    endif

    # settle on a map grid ASAP
    if( "$GRID" == "" ) then
        set GRID = `echo | mapdump mapin ${tempfile}sum.map | awk '/Grid sampling/{print $(NF-2), $(NF-1), $NF; exit}'`
        if(! $quiet) echo "using map grid: $GRID"
    endif



    # same for coordinate map
    if( "$coordmap" != "" ) then
        if(! -e ${tempfile}coordsum.map) then
            cp -p ${sharedtmp}coord_seed${seed}.map ${tempfile}coordsum.map
        else
            echo maps add |\
            mapmask mapin1 ${sharedtmp}coord_seed${seed}.map mapin2 ${tempfile}coordsum.map \
               mapout ${tempfile}new.map >> $logfile        
            mv ${tempfile}new.map ${tempfile}coordsum.map
        endif
    endif


    # and carbonized map
    if( "$carbonmap" != "" ) then
        if(! -e ${tempfile}carbonsum.map) then
            cp -p ${sharedtmp}carbonized_seed${seed}.map ${tempfile}carbonsum.map
        else
        echo maps add |\
            mapmask mapin1 ${sharedtmp}carbonized_seed${seed}.map mapin2 ${tempfile}carbonsum.map \
               mapout ${tempfile}new.map >> $logfile        
           mv ${tempfile}new.map ${tempfile}carbonsum.map
        endif
    endif

    # and now pdb files
    if( "$multiconfpdb" != "" ) then
        if(! -e ${tempfile}multiconf.pdb) then
            egrep -v "^END" ${sharedtmp}seed${seed}out.pdb >! ${tempfile}multiconf.pdb
        else
            # just append atoms
            egrep "^ATOM|^HETAT" ${sharedtmp}seed${seed}out.pdb >> ${tempfile}multiconf.pdb
        endif
    endif


    # clean up
    if(! $debug) then
        rm -f ${sharedtmp}mask_seed${seed}.map ${sharedtmp}seed${seed}done.txt
        rm -f ${sharedtmp}coord_seed${seed}.map
        rm -f ${sharedtmp}carbonized_seed${seed}.map
    endif
    # sumarize rms geometry deviations
    if("$CLUSTER" == "SGE") then
        mv refmac_cpu.com.*.$seed ${sharedtmp}seed${seed}.log
    endif
    while (! -e ${sharedtmp}seed${seed}.log)
        echo "waiting for $seed log"
        sleep 1
        if("$CLUSTER" == "SGE")  mv refmac_cpu.com.*.$seed ${sharedtmp}seed${seed}.log
        find ${sharedtmp}seed${seed}.log
    end

    awk '/rmsBOND/,/Final/ && NF>1' ${sharedtmp}seed${seed}.log |\
    awk 'NF>1 && ! /[a-z]/' |\
    tail -n 1 >> ${tempfile}rmsds.log
end
echo ""
# wait for SMP jobs
wait

if(1 || $debug) then
    rmsd2B ${sharedtmp}*out.pdb | egrep "^ATOM|^HETAT" >! rmsd2B.pdb

    awk '$NF!="H"' $pdbfile rmsd2B.pdb | rmsd -v debug=1 >! rmsds.txt
endif

# examine RMS deviations?
echo -n "rms bonds, angles = "
awk '{print $7,$9}' ${tempfile}rmsds.log |\
awk '{sum1+=$1*$1;sum2+=$2*$2} \
   END{print sqrt(sum1/NR),sqrt(sum2/NR)}'

cp -p ${tempfile}rmsds.log rmsds.log

# wait for queued jobs?

if(! $debug) then
    rm -f ${tempfile}*seed*.log >& /dev/null
    rm -f ${tempfile}*seed*.pdb >& /dev/null
    rm -f seed*.log >& /dev/null
    rm -f qsubs.log >& /dev/null
    rm -f refmac_cpu.com.* >& /dev/null
endif

# save as readable name
if("$multiconfpdb" != "") then
    mv ${tempfile}multiconf.pdb $multiconfpdb
endif

# final check.  Did it work
if(! -e ${tempfile}sum.map) then
    set BAD = "map summation failed"
    goto exit
endif

# put on the appropriate scale, and expand to cell
set max = `echo | mapdump mapin ${tempfile}sum.map | awk '/Maximum density/{print $NF}'`
#set max = $seeds
set scale = `echo $ksol $max | awk '$2+0==0{print $1;exit} {print $1/$2}'`

set axis = "Z X Y"

reaxis:
# reorganize map for SFALL
mapmask mapin ${tempfile}sum.map mapout ${tempfile}ksol.map << EOF >> $logfile
scale factor $scale 0
xyzlim cell
AXIS $axis
EOF
if("$mapout" != "") then
    mapmask mapin ${tempfile}ksol.map mapout $mapout << EOF >> $logfile
    xyzlim asu
    AXIS X Y Z
EOF
    echo -n "fuzzy solvent map: "
    echo | mapdump mapin $mapout | awk '/density/{gsub("imum","");printf("%s %s ",$1,$NF)}'
    echo "$mapout"
endif


if( "$coordmap" == "" ) goto skipcoordmap

# put on the appropriate scale, and expand to cell
set scale = `echo $seeds | awk '{print 1/$1}'`

# reorganize map for SFALL
mapmask mapin ${tempfile}coordsum.map mapout ${tempfile}coordavg.map << EOF >> $logfile
scale factor $scale 0
xyzlim cell
AXIS $axis
EOF
mapmask mapin ${tempfile}coordavg.map mapout $coordmap << EOF >> $logfile
xyzlim asu
axis X Y Z
EOF
echo -n "fuzzy coord   map: " 
echo | mapdump mapin $coordmap | awk '/density/{gsub("imum","");printf("%s %s ",$1,$NF)}'
echo "$coordmap"


# do the same for the original coord
cat $pdbfile |\
awk '! /^ATOM|HETAT/{print;next} {Ee=substr($0,13,2);gsub(Ee," ","")}\
      Ee !~ /^[CNOHS]$/{Ee=toupper(substr($NF,1,2))} {l=$0;\
      printf("%s%2s%s\n",substr(l,1,12),Ee,substr(l,15,57))}' |\
cat >! ${tempfile}sfallme.pdb
sfall xyzin ${tempfile}sfallme.pdb mapout ${tempfile}sfalled.map << EOF >> $logfile
mode atmmap
SYMM $SG
cell $CELL
GRID $GRID
SFSG 1
EOF
mapmask mapin ${tempfile}sfalled.map mapout ${tempfile}orig_coord.map << EOF >> $logfile
xyzlim cell
AXIS $axis
EOF
cp ${tempfile}orig_coord.map orig_coord.map


echo scale factor -1 |\
mapmask mapin ${tempfile}orig_coord.map mapout ${tempfile}neg.map >> $logfile
echo maps add |\
mapmask mapin1 ${tempfile}coordavg.map mapin2 ${tempfile}neg.map mapout coord_diff.map >> $logfile
echo -n "fz-orig coord map: "
echo | mapdump mapin coord_diff.map | awk '/density/{gsub("imum","");printf("%s %s ",$1,$NF)}'
echo "coord_diff.map"


# compare fuzzy coordinate map to original model
echo correlate section |\
overlapmap mapin1 ${tempfile}orig_coord.map \
           mapin2 ${tempfile}coordavg.map |\
awk '/Total corr/{print "CC between original and fuzzy coordinate maps:", $NF}'


# compare fuzzy coordinate map to original model as structure factors
sfall mapin ${tempfile}orig_coord.map hklout ${tempfile}orig_coord.mtz << EOF >> $logfile
mode sfcalc mapin
SFSG 1
resolution $reso
EOF
sfall mapin ${tempfile}coordavg.map hklout ${tempfile}fuzzed.mtz << EOF >> $logfile
mode sfcalc mapin
SFSG 1
resolution $reso
EOF
rm -f ${tempfile}cadded.mtz
sftools << EOF >> $logfile
read ${tempfile}orig_coord.mtz col FC
read ${tempfile}fuzzed.mtz col FC
set labels
Forig
Ffuzz
calc Q col SIGF = 0.1
write ${tempfile}cadded.mtz
quit
y
EOF
scaleit hklin ${tempfile}cadded.mtz hklout ${tempfile}scaled.mtz << EOF | tee ${tempfile}scaleit.log >> $logfile
labin FP=Forig SIGFP=SIGF FPH1=Ffuzz SIGFPH1=SIGF
EOF
awk '/THE TOTALS/{print "fuzz vs orig coord-maps R factor:", substr($0,50,6)}' ${tempfile}scaleit.log
set scale = `awk '$1=="Derivative" && ! /itle/{print $3}' ${tempfile}scaleit.log | tail -1`
set B     = `awk '/equivalent iso/{print $NF}' ${tempfile}scaleit.log | tail -1`
echo "scale= $scale B= $B"

skipcoordmap:




if( "$carbonmap" == "" ) goto skipcarbonmap

# put on the appropriate scale, and expand to cell
set scale = `echo $seeds | awk '{print 1/$1}'`

# reorganize map for SFALL
mapmask mapin ${tempfile}carbonsum.map mapout ${tempfile}carbonavg.map << EOF >> $logfile
scale factor $scale 0
xyzlim cell
AXIS $axis
EOF
mapmask mapin ${tempfile}carbonavg.map mapout $carbonmap << EOF >> $logfile
xyzlim asu
axis X Y Z
EOF
echo -n "fuzzy carbon  map: " 
echo | mapdump mapin $carbonmap | awk '/density/{gsub("imum","");printf("%s %s ",$1,$NF)}'
echo "$carbonmap"


# carbonize the original coord
cat $pdbfile |\
awk '! /^ATOM|HETAT/{print;next} $NF!="H"{l=$0\
      printf("%s C%s 5.00           C\n",substr(l,1,12),substr(l,15,47),substr(l,1,61))}' |\
cat >! ${tempfile}carbonized.pdb
sfall xyzin ${tempfile}carbonized.pdb mapout ${tempfile}sfalled.map << EOF >> $logfile
mode atmmap
SYMM $SG
cell $CELL
GRID $GRID
SFSG 1
EOF
mapmask mapin ${tempfile}sfalled.map mapout ${tempfile}orig_carbon.map << EOF >> $logfile
xyzlim cell
AXIS $axis
EOF
cp ${tempfile}orig_carbon.map orig_carbon.map


echo scale factor -1 |\
mapmask mapin ${tempfile}orig_carbon.map mapout ${tempfile}neg.map >> $logfile
echo maps add |\
mapmask mapin1 ${tempfile}carbonavg.map mapin2 ${tempfile}neg.map mapout carbon_diff.map >> $logfile
echo -n "fz-orig carbon map: "
echo | mapdump mapin carbon_diff.map | awk '/density/{gsub("imum","");printf("%s %s ",$1,$NF)}'
echo "carbon_diff.map"


# compare fuzzy carbonized map to original model
echo correlate section |\
overlapmap mapin1 ${tempfile}orig_carbon.map \
           mapin2 ${tempfile}carbonavg.map |\
awk '/Total corr/{print "CC between original and fuzzy carbonized coordinate maps:", $NF}'


# compare fuzzy coordinate map to original model as structure factors
sfall mapin ${tempfile}orig_carbon.map hklout ${tempfile}orig_carbon.mtz << EOF >> $logfile
mode sfcalc mapin
SFSG 1
resolution $reso
EOF
sfall mapin ${tempfile}carbonavg.map hklout ${tempfile}fuzzed.mtz << EOF >> $logfile
mode sfcalc mapin
SFSG 1
resolution $reso
EOF
rm -f ${tempfile}cadded.mtz
sftools << EOF >> $logfile
read ${tempfile}orig_carbon.mtz col FC
read ${tempfile}fuzzed.mtz col FC
set labels
Forig
Ffuzz
calc Q col SIGF = 0.1
write ${tempfile}cadded.mtz
quit
y
EOF
scaleit hklin ${tempfile}cadded.mtz hklout ${tempfile}scaled.mtz << EOF | tee ${tempfile}scaleit.log >> $logfile
labin FP=Forig SIGFP=SIGF FPH1=Ffuzz SIGFPH1=SIGF
EOF
awk '/THE TOTALS/{print "fuzz vs orig carbonized-maps R factor:", substr($0,50,6)}' ${tempfile}scaleit.log
set scale = `awk '$1=="Derivative" && ! /itle/{print $3}' ${tempfile}scaleit.log | tail -1`
set B     = `awk '/equivalent iso/{print $NF}' ${tempfile}scaleit.log | tail -1`
echo "scale= $scale B= $B"

skipcarbonmap:














sfall:
# first one will crash?
echo "sfall..."
sfall mapin ${tempfile}ksol.map hklout ${tempfile}sfalled.mtz << EOF | tee -a ${tempfile}crash.log >> $logfile
mode sfcalc mapin
SFSG 1
resolution $reso
EOF

# try to recover from SFALL crash?
if(! -e ${tempfile}sfalled.mtz && ! $?newaxis) then
    # get the axis order that SFALL wants
    set newaxis = `awk '/Check Iuvw/{gsub("1","X");gsub("2","Y");gsub("3","Z");print $6,$7,$8}' ${tempfile}crash.log`
    if("$newaxis" != "$axis") then
        set axis = ( $newaxis )
        goto reaxis
    endif
endif

if(! -e ${tempfile}sfalled.mtz) then
    cat ${tempfile}crash.log
    set BAD = "SFALL failed"
    goto exit
endif

# combine with refinement mtz
echo "combining $mtzfile with Fpart PHIpart"
echo "" |\
mtzdump hklin $mtzfile |\
awk 'NF>10' | awk '$(NF-1)~/^[FDQIJPWGKLMAR]$/{++n;print $NF" "}' |\
egrep -v "part" |\
awk '{++n;print "E"n"="$1}' >! ${tempfile}tokens.txt
set tokens = `cat ${tempfile}tokens.txt`

cad hklin1 $mtzfile hklin2 ${tempfile}sfalled.mtz \
    hklout ${tempfile}cadded.mtz << EOF >> $logfile
labin file 1 $tokens
labin file 2 E1=FC E2=PHIC
labou file 2 E1=Fpart E2=PHIpart
EOF
if(! -e ${tempfile}cadded.mtz) then
    set BAD = "failed to make $mtzout"
    goto exit
endif
mv ${tempfile}cadded.mtz $mtzout
if(-e $mtzout) then
    echo "output mtz file: $mtzout"
else
    set BAD = "failed to make $mtzout"
    goto exit
endif




# compare to F in mtz file?






final_message:
cat << EOF

add this to refmac input:
LABIN FP=FP SIGFP=SIGFP FREE=FreeR_flag  FPART1=Fpart PHIP1=PHIpart
SCPART 1
SOLVENT NO
EOF


exit:

if($?BAD) then
    echo "ERROR: $BAD"
    exit 9
endif

if($debug) exit

rm -f ${tempfile}* >& /dev/null

exit




















#========================================================================
#    initial setup routines down here and out of the way
Setup:

# default file names
set pdbfile = ""
set mtzfile = ""
set libfile = ""

set quiet = 0
set debug = 0

# abort if we cant run
if(! $?CCP4_SCR) then
    set BAD = "CCP4 is not set up."
    goto exit
endif

# pick temp filename location
set logfile = /dev/null
mkdir -p ${CCP4_SCR} >&! /dev/null
set tempfile = ${CCP4_SCR}/fuzzymask$$temp
if(-w /dev/shm/ ) then
    mkdir -p /dev/shm/${USER}
    set tempfile = /dev/shm/${USER}/fuzzymask$$temp
endif
if(-w /scratch/ ) then
    set tempfile = /scratch/${USER}fuzzymask$$temp
    mkdir -p /scratch/${USER}
endif

# some platforms dont have these?
if(! $?USER) then
    setenv USER `whoami`
endif

#========================================================================
#    Reading command-line arguments:

foreach arg ( $* )
    if("$arg" =~ *.pdb) then
        if(! -e "$arg") then
            echo "ERROR: $arg does not exist"
            exit 9
        endif
        set pdbfile = "$arg"
        continue
    endif
    if("$arg" =~ *.map && "$arg" !~ *=*) then
        set mapout = "$arg"
    endif
    if("$arg" =~ *.mtz && "$arg" !~ *=*) then
        if(! -e "$arg") then
            echo "ERROR: $arg does not exist"
            exit 9
        endif
        set mtzfile = "$arg"
    endif
    if("$arg" =~ *.cif || "$arg" =~ *.lib) then
        if(! -e "$arg") then
            echo "ERROR: $arg does not exist"
            exit 9
        endif
        set libfile = "$arg"
    endif
    if("$arg" == debug) then
        set debug = 1
    endif
    if("$arg" == quiet) then
        set quiet = 1
    endif
    if("$arg" =~ ksol=*) set ksol = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ CPUs=*) set user_CPUs = `echo $arg | awk -F "=" '{print $2}'`
    if("$arg" =~ CLUSTER=*) set CLUSTER = `echo $arg | awk -F "=" '{print $2}'`
    if("$arg" =~ seeds=*) set seeds = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ frac_thrubond=*) set frac_thrubond = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ ncyc_thrubond=*) set ncyc_thrubond = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ frac_magnforce=*) set frac_magnforce = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ ncyc_magnforce=*) set ncyc_magnforce = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ ncyc_ideal=*) set ncyc_ideal = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ ncyc_data=*) set ncyc_data = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ vdwprobe=*) set vdwprobe = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ ionprobe=*) set ionprobe = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ rshrink=*) set rshrink = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ kick_steps=*) set kick_steps = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ kick_scale=*) set kick_scale = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ drykick_scale=*) set drykick_scale = `echo $arg | awk -F "=" '{print $2+0}'`
    if("$arg" =~ mtzout=*mtz) set mtzout = `echo $arg | awk -F "=" '{print $2}'`
    if("$arg" =~ output=*mtz) set mtzout = `echo $arg | awk -F "=" '{print $2}'`
    if("$arg" =~ mapout=*) set mapout = `echo $arg | awk -F "=" '{print $2}'`
    if("$arg" =~ coordmap=*) set coordmap = `echo $arg | awk -F "=" '{print $2}'`
    if("$arg" =~ pdbfile=*) set pdbfile = `echo $arg | awk -F "=" '{print $2}'`
    # this does not work
    if("$arg" =~ grid=*) set GRID = `echo $arg | awk -F "=" '{print $2}' | awk -F "," '{print $1,$2,$3}'`
end

# bug out here if we cant run
if(! -e "$pdbfile") then
    echo "please specify an existing PDB file."
    goto Help
endif
if(! -e "$mtzfile") then
    echo "please specify an existing MTZ file."
    goto Help
endif

if("$coordmap" !~ *.map && "$coordmap" !~ *.ccp4) then
    set coordmap = ""
endif

# examine the MTZ file
echo | mtzdump hklin $mtzfile >&! ${tempfile}mtzdump.txt
set reso = `awk '/Resolution Range/{getline;getline;print $6}' ${tempfile}mtzdump.txt`
rm -f ${tempfile}mtzdump.txt
if("$reso" == "") then
    echo "cannot read $mtzfile as mtz"
    goto Help
endif


# get cell dimensions and space group
set CELL = `awk '/CRYST1/{print $2,$3,$4,$5,$6,$7}' $pdbfile`
set pdbSG = `awk '/^CRYST/{print substr($0,56,12)}' $pdbfile | head -1`
if("$pdbSG" == "R 32") set pdbSG = "R 3 2"
if("$pdbSG" == "P 21") set pdbSG = "P 1 21 1"
if("$pdbSG" == "R 3 2" && $CELL[6] == 120.00) set pdbSG = "H 3 2"
set SG = `awk -v pdbSG="$pdbSG" -F "[\047]" 'pdbSG==$2{print;exit}' ${CLIBD}/symop.lib | awk '{print $4}'`
if("$SG" == R3 && $CELL[6] == 120.00) set SG = H3
if("$SG" == "") set SG = P1



# decide on a grid spacing
set default_res = `awk '/^REMARK   2 RESOLUTION/{print $4}' $pdbfile`
if("$default_res" == "") then
    set minB = `awk '/^ATOM/ || /^HETATM/{print substr($0, 61, 6)+0}' $pdbfile | sort -n | awk '{v[++n]=$1} END{print v[int(n/2)]}'`
    set default_res = `echo $minB | awk '$1>0{print 3*sqrt($1/80)}'`
endif
if("$default_res" == "") set default_res = 1.5

# reduce default resolution if it will overload sftools
set megapoints = `echo $CELL $default_res | awk '{printf "%d", 27*($1 * $2 * $3)/($NF*$NF*$NF)/900000}'`
if("$megapoints" > 100) then
    set default_res = `echo $CELL 100000000 | awk '{print ($1*$2*$3*27/$NF)^(1/3)}'`
endif
if("$reso" == "") then
    set reso = $default_res
endif




# system setup
if(! $?CLUSTER) set CLUSTER = ""
set pwd = `pwd`
set uname = `uname`

if("$CLUSTER" != "" && "$CLUSTER" != "none") goto cluster_detected

# test for torq cluster
cat << EOF >! test$$.csh
#! /bin/tcsh -f
#\$ -r n                            #-- if job crashes, do not restart
#\$ -l mem_free=1M                  #-- submits on nodes with enough free memory (required)
#\$ -l arch=linux-x64               #-- SGE resources (CPU type)
#\$ -l netapp=1M,scratch=1M         #-- SGE resources (home and scratch disks)
#\$ -l h_rt=00:00:01                #-- runtime limit
touch \$1
EOF
chmod a+x test$$.csh

# test for TORQ cluster
set queued = 0
qsub -e /dev/null -o /dev/null -d $pwd test$$.csh -F "testtorq$$.txt" >& /dev/null
if(! $status) set queued = 1
set timeout = 100
while ( $queued && $timeout )
    @ timeout = ( $timeout - 1 )
    sleep 0.1
    set queued = `qstat |& awk '$(NF-1)~/[RQ]/{print}' | wc -l`
end
if(-e ${pwd}/testtorq$$.txt) then
    echo "detected torq-based cluster"
    set CLUSTER = TORQ
    goto cluster_detected
endif

# test for SGE
qsub -cwd test$$.csh testsge$$.txt >&! ${tempfile}jid.txt
if(! $status) set queued = 1
set jid = `awk '$3+0>0{print $3}' ${tempfile}jid.txt | tail -n 1`
if("$jid" != "") echo "testing SGE cluster, to skip set environment CLUSTER=SGE"
set timeout = 300
while ( ! -e ${pwd}/testsge$$.txt && $queued && $timeout )
    @ timeout = ( $timeout - 1 )
    set queued = `qstat -j $jid |& awk 'NF>1 && ! /jobs do not exist/' | wc -l`
    sleep 1
end
if("$timeout" == "0" && ! -e ${pwd}/testsge$$.txt) echo "SGE cluster not working."
if(-e ${pwd}/testsge$$.txt) then
    echo "detected SGE-based cluster"
    set CLUSTER = SGE
    goto cluster_detected
endif

cluster_detected:
rm -f test$$.csh* testtorq$$.txt testsge$$.txt >& /dev/null
rm -f ${tempfile}jid.txt >& /dev/null



if("$uname" == "Linux") then
    set CPUs = `awk '/^processor/' /proc/cpuinfo | wc -l`
    set chips = `awk '/^physical/' /proc/cpuinfo | sort -u | wc -l`
    set cores_per_chip = `awk '/^core/' /proc/cpuinfo | sort -u | wc -l`
    set cores = `echo $chips $cores_per_chip | awk '{print $1*$2}'`
    set threads_per_chip = `awk '/^siblings/{print $NF;exit}' /proc/cpuinfo`
    set threads_per_core = `echo $threads_per_chip $cores_per_chip | awk '{threads=int($1/$2+0.001)} threads==0{threads=1} {print threads}'`
#    echo "found $CPUs CPUs on $chips chips with $cores_per_chip cores each ($threads_per_core threads/core)"

    set freeCPUs = `w | cat - /proc/cpuinfo - | awk '/^processor/{++p} /load aver/{l=$(NF-2)+0} END{print int(p-l+0.5)}'`
#    echo "found $freeCPUs free CPUs"

    # guestimate optimal memory usage
    set freemem = `free -m | awk '/^Mem:/{print $4+$7;exit}'`
    set totalmem = `free -m | awk '/^Mem:/{print $2;exit}'`
#    echo "$totalmem MB RAM total"
endif
if("$uname" == "Darwin") then
    # for some reason on macs: don't have wget
    alias wget 'curl -o `basename \!:1` \!:1'

    set cores = `sysctl hw.physicalcpu | awk '/^hw./{print $NF}'`
    set CPUs = `sysctl hw.logicalcpu | awk '/^hw./{print $NF}'`
    set threads = `echo $CPUs $cores | awk '{threads=int($1/$2+0.001)} threads==0{threads=1} {print threads}'`
    set chips = 1
    set cores_per_chip = `echo $cores $chips | awk '{print $1/$2}'`
    set threads_per_chip = `echo $threads $chips | awk '{print $1/$2}'`
    set threads_per_core = `echo $threads $cores | awk '{print $1/$2}'`
#    echo "found $CPUs CPUs in $cores cores ($threads threads/core)"

    set freeCPUs = `w | awk -v p=$CPUs '/load aver/{l=$(NF-2)+0} END{print int(p-l+0.5)}'`
#    echo "found $freeCPUs free CPUs"

    # guestimate optimal memory usage
    set totalmem = `sysctl hw.memsize | awk '/hw.memsize/ && $NF+0>0{print $NF/1024/1024;exit}'`
#    echo "$totalmem MB RAM available"
endif
if(! $?CPUs) then
    ehco "WARNING: unknown platform! "
endif
# allow user to override
if("$user_CPUs" == "cores") set user_CPUs = $cores
if("$user_CPUs" == "free") set user_CPUs = $freeCPUs
if("$user_CPUs" == "all") set user_CPUs = $CPUs
if("$user_CPUs" != "auto") set CPUs = $user_CPUs
if("$CLUSTER" == "none") set CLUSTER = ""
if("$CLUSTER" == "") echo "will use $CPUs CPUs"



# change things for debug mode
if($debug) then
    set tempfile = ./tempfile
    set logfile = fuzzymask_debug.log
endif




# see if there is a user-provided parameter file
set otheropts
if(-e refmac_opts.txt) then
    awk '/^occupa/{next} {print}' refmac_opts.txt >! ${tempfile}refmac_opts.txt
    set otheropts = "@${tempfile}refmac_opts.txt"
endif
set LIBSTUFF
if(-e "$libfile") set LIBSTUFF = "LIBIN $libfile"
if(-e atomsf.lib) set LIBSTUFF = "$LIBSTUFF ATOMSF ./atomsf.lib"








# requires jigglepdb.awk
set path = ( `dirname $0` . $path )


# deploy scripts?
set test = `echo | jigglepdb.awk | awk '/jiggled by dXYZ/{print 1}'`
if("$test" == "1") goto after_Setup

echo "deploying jigglepdb.awk ..."
cat << EOF-script >! jigglepdb.awk
#! `which awk` -f
#! /bin/awk -f
#
#
#        Jiggles a pdb file's coordinates by some random value
#        run like this:
#
#        jigglepdb.awk -v seed=2343 -v shift=1.0 old.pdb >! jiggled.pdb
#         (use a different seed when you want a different output file)
#
BEGIN {

    if(! shift)  shift = 0.5
    if(! Bshift) Bshift = shift
    if(shift == "byB") Bshift = 0
    if(shift == "Lorentz") Bshift = 0
    if(! shift_scale) shift_scale = 1
    if(! dry_shift_scale) dry_shift_scale = 1
    pshift = shift
    shift_opt = shift
    if(pshift == "byB") pshift = "sqrt(B/8)/pi"
    if(pshift == "LorentzB") pshift = "Lorentzian B"
    if(seed) srand(seed+0)
    if(! keepocc) keepocc=0
    if(! distribution) distribution="gaussian";
    if(! frac_thrubond) frac_thrubond=0.5
    if(! ncyc_thrubond) ncyc_thrubond=0
    if(! frac_magnforce) frac_magnforce=5.0
    if(! ncyc_magnforce) ncyc_magnforce=ncyc_thrubond/2

    pi=4*atan2(1,1);

    # random number between 1 and 0 to select conformer choices
    global_confsel=rand();

    print "REMARK jiggled by dXYZ=", pshift, "dB=", Bshift
    print "REMARK shift_scale=",shift_scale,"dry_shift_scale=",dry_shift_scale
    print "REMARK frac_thrubond=",frac_thrubond,"ncyc_thrubond=",ncyc_thrubond
    print "REMARK frac_magnforce=",frac_magnforce,"ncyc_magnforce=",ncyc_magnforce
    print "REMARK random number seed: " seed+0
}

# count all lines
{++n;line[n]=\$0}

/^ATOM|^HETAT/{
    if(debug) print tolower(\$0)

#######################################################################################
#    electrons = substr(\$0, 67,6)
#    XPLORSegid = substr(\$0, 73, 4)            # XPLOR-style segment ID
#    split(XPLORSegid, a)
#    XPLORSegid = a[1];
    Element = substr(\$0, 67)

#    Atomnum= substr(\$0,  7, 5)+0
    if(Element !~ /^[A-Z]/) Element= substr(\$0, 13, 2);
    Greek= substr(\$0, 15, 2);
    split(Element Greek, a)
    Atomtype[n]   = a[1];
    if(length(a[1])==4 && Element ~ /^H/)Element="H";
    gsub(" ","",Element);
    Ee[n] = Element;
    #main[n]=(Atom[n] ~ /^[NCO]\$/ || Atom[n] ~ /^C[AB]\$/);
    prefix[n] = substr(\$0,  1,30)
    Conf[n]   = substr(\$0, 17, 1)                # conformer letter
    Restyp[n] = substr(\$0, 18, 3)
    Segid[n]  = substr(\$0, 22, 1)            # O/Brookhaven-style segment ID
    Resnum[n] = substr(\$0, 23, 4)+0
    X[n]      = substr(\$0, 31, 8)+0
    Y[n]      = substr(\$0, 39, 8)+0
    Z[n]      = substr(\$0, 47, 8)+0
    Occ[n]    = substr(\$0, 55, 6)+0
    Bfac[n]   = substr(\$0, 61, 6)+0
    rest[n]   = substr(\$0, 67)
#    ATOM   = toupper(substr(\$0, 1, 6))
#######################################################################################
}

END{

  # min/max bond lengths
  if(min_bond_d == "") min_bond_d = 1
  if(max_bond_d == "") max_bond_d = 2


  for(i=1;i<=n;++i)
  {
    if(prefix[i] !~ /^ATOM|^HETAT/ ) continue;

    if(shift_opt=="byB" || shift_opt=="LorentzB"){
        # switch on "thermal" shift magnitudes
        shift=sqrt(Bfac[i]/8)/pi*sqrt(3);

        # kick them more than byB?
        if(shift_scale != 1){
            shift *= shift_scale;
        }

        # kick them more if they are not water
        if(Restyp[i] != "HOH" && dry_shift_scale != 1){
            shift *= dry_shift_scale;
        }

        # randomly "skip" conformers with occ<1
        if(Occ[i]+0<1){
            # remember all occupancies
            if(conf_hi[Conf[i],Segid[i],Resnum[i]]==""){
                conf_lo[Conf[i],Segid[i],Resnum[i]]=cum_occ[Segid[i],Resnum[i]]+0;
                cum_occ[Segid[i],Resnum[i]]+=Occ[i];
                conf_hi[Conf[i],Segid[i],Resnum[i]]=cum_occ[Segid[i],Resnum[i]];
            }
        }
    }
    if(shift_opt == "LorentzB")
    {
        distribution = "Lorentz"
    }
    if(distribution == "Lorentz")
    {
        dx = lorentzrand(shift/sqrt(3));
        dy = lorentzrand(shift/sqrt(3));
        dz = lorentzrand(shift/sqrt(3));
    }
    
    if(distribution == "gaussian" || distribution == "Gauss")
    {
        dx = gaussrand(shift/sqrt(3));
        dy = gaussrand(shift/sqrt(3));
        dz = gaussrand(shift/sqrt(3));
    }
    if(distribution == "uniform")
    {
        dR=2
        while(dR>1)
        {
            dx = (2*rand()-1);
            dy = (2*rand()-1);
            dz = (2*rand()-1);
            dR = sqrt(dx^2+dy^2+dz^2);
        }
        dx *= shift;
        dy *= shift;
        dz *= shift;
    }

    dX[i] = dx;
    dY[i] = dy;
    dZ[i] = dz;

    # pick a random shift on B-factor
    if(Bshift+0>0) Bfac[i] += gaussrand(Bshift)
    if(Oshift+0>0) Occ[i] += gaussrand(Oshift)
    
    # use same occopancy for given conformer
    if(! keepocc && conf_hi[Conf[i],Segid[i],Resnum[i]]!=""){
        # use same random number for all conformer choices
        confsel = global_confsel;
        # unless occupancies do not add up
        if(Conf[i]==" "){
            # save this for later?
            confsel = rand();
        }
        Occ[i] = 0;
        # atom only exists if it falls in the chosen interval
        lo=conf_lo[Conf[i],Segid[i],Resnum[i]];
        hi=conf_hi[Conf[i],Segid[i],Resnum[i]];
        if(lo < confsel && confsel <= hi) Occ[i]=1;
    }
  }


    if(frac_thrubond != 0 && ncyc_thrubond != 0)
    {
        print "REMARK measuring original bond lengths"
        # now, find all the bonds in the unperturbed structure that don't involve zero-occupancy members
        for(i=1;i<=n;++i){
            #skip zero occupancy
            if( Occ[i]==0 ) continue;
            for(j=1;j<i;++j){
                # skip zero occupancy
                if( Occ[j]==0 ) continue;
                # skip nonsencial conformer relationships?
#                if( Conf[j] != Conf[i] && ! ( Conf[i] == " " || Conf[j] == " " || main[j] && main[i]) ) continue;
                mind=min_bond_d;maxd=max_bond_d;
                # hydrogen bonding lengths are shorter
                if( Ee[i] == "H" || Ee[j] == "H" ){
                    mind=0.5;maxd=1.5;
                    if(Ee[i]==Ee[j])maxd=0;
                };
                # recognize disulfides
                if( Ee[i]=="S" && Ee[j] == "S" && Restyp[i]=="CYS" && Restyp[j]=="CYS"){mind=1.5;maxd=2.5};
                # measure distance
                if(X[i]>X[j]+maxd || X[i]<X[j]-maxd) continue;
                if(Y[i]>Y[j]+maxd || Y[i]<Y[j]-maxd) continue;
                if(Z[i]>Z[j]+maxd || Z[i]<Z[j]-maxd) continue;
                d=sqrt((X[i]-X[j])^2+(Y[i]-Y[j])^2+(Z[i]-Z[j])^2);
                if(d>mind && d<maxd) {
                    # distance falls within range
                    newbond=1;
                    for(k=1;k<=nbonds[i]+0;++k){
                        # woops, this bond already exists
                        if(bond[i,k]==j){
                            newbond=0;
                            break;
                        }
                    }
                    if(newbond) {
                        # legitimate new bond, increment the list
                        ++nbonds[i];
                        bond[i,nbonds[i]]=j;
                        bondlen[i,j]=d;
                    }
                    newbond=1;
                    # check if reverse bond is already there
                    for(k=1;k<=nbonds[j]+0;++k){
                        if(bond[j,k]==i){
                            newbond=0;
                            break;
                        }
                    }
                    if(newbond) {
                        # also account for reverse bond
                        ++nbonds[j];
                        bond[j,nbonds[j]]=i;
                        bondlen[j,i]=d;
                    }
                }
            }
        }

        # now go through and force ideal bonds that we know about
        for(i=1;i<=n;++i){
            if(Occ[i]==0) continue;
        }

        # measure and store all currently applied displacement magnitudes
        for(i=1;i<=n;++i){
            master_dXYZ[i]=sqrt(dX[i]**2+dY[i]**2+dZ[i]**2);
        }

        # iterative application of thru-bond smoothing
        for(l=1;l<=ncyc_thrubond;++l)
        {
            print "REMARK thru-bond averaging cycle",l
            # loop over all atoms
            for(i=1;i<=n;++i){
                if(Occ[i]==0) continue;
                realbonds=ddX[i]=ddY[i]=ddZ[i]=0;
                # take every atom bonded to this atom
                for(u=1;u<=nbonds[i];++u){
                j = bond[i,u];
                if(i==j) continue;
                if(Occ[j]==0)continue;
                ++realbonds;
                ddX[i] += dX[j];
                ddY[i] += dY[j];
                ddZ[i] += dZ[j];
            }
            if(! realbonds) continue;
                ddX[i] /= realbonds;
                ddY[i] /= realbonds;
                ddZ[i] /= realbonds;
            }
            # second pass to update delta-positions
            for(i=1;i<=n;++i){
                if(Occ[i]==0) continue;
                bw=frac_thrubond;
                if(Ee[i]=="H") bw=1;
                dX[i] = (1-bw)*dX[i] + bw*ddX[i];
                dY[i] = (1-bw)*dY[i] + bw*ddY[i];
                dZ[i] = (1-bw)*dZ[i] + bw*ddZ[i];
            }
            # third pass to rescale shift magnitudes
            if(l<=ncyc_magnforce)
            {
                mw=(1-(l-1)/(ncyc_magnforce-1))
                ms=frac_magnforce
                if(l>ncyc_magnforce)mw=0;
                for(i=1;i<=n;++i){
                    if(Occ[i]==0) continue;
                    if(Ee[i]=="H") mw=0;
                    mag = sqrt(dX[i]**2+dY[i]**2+dZ[i]**2);
                    if(mag<=0.0) {
                        continue;
                        # make something up?
                        mag=1e-6;
                        dX[i]=mag*(rand()-0.5);
                        dY[i]=mag*(rand()-0.5);
                        dZ[i]=mag*(rand()-0.5);
                    }
                    # get difference between current and scaled version of originally prescribed shift magnitude
                    dmag = (ms*master_dXYZ[i]-mag);
                    # scale thet shift to be more like original magnitude
                    scale = (mag+mw*dmag)/mag;
                    dX[i] = scale*dX[i];
                    dY[i] = scale*dY[i];
                    dZ[i] = scale*dZ[i];
#if(i==456) print "GOTHERE1",dX[i],dY[i],dZ[i],"    ",master_dXYZ[i],mag,"    ",dmag,scale,mw
                }
            }
        }
    }


    for(i=1;i<=n;++i)  
    {
        if(prefix[i] !~ /^ATOM|^HETAT/ )
        {
            print line[i];
            continue;
        }
      
        X[i] += dX[i];
        Y[i] += dY[i];
        Z[i] += dZ[i];
      
        # now print out the new atom
        printf("%s%8.3f%8.3f%8.3f %5.2f%6.2f%s\\n",prefix[i],X[i],Y[i],Z[i],Occ[i],Bfac[i],rest[i]);        
    }
}



#######################################################################################
# function for producing a random number on a gaussian distribution
function gaussrand(sigma){
    if(! sigma) sigma=1
    rsq=0
    while((rsq >= 1)||(rsq == 0))
    {
        x=2.0*rand()-1.0
        y=2.0*rand()-1.0
        rsq=x*x+y*y
    }
    fac = sqrt(-2.0*log(rsq)/rsq);
    return sigma*x*fac
}

# function for producing a random number on a Lorentzian distribution
function lorentzrand(fwhm){
    if(! fwhm) fwhm=1

    return fwhm/2*tan(pi*(rand()-0.5))
}

function tan(x){
    return sin(x)/cos(x)
}
EOF-script
chmod a+x jigglepdb.awk


goto after_Setup



