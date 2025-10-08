#! /bin/tcsh -f
#
# remove bulk solvent from Fobs  - ala "squeeze"   - James Holton  9-17-25
#
# requires Tc_maskify.com map_func.com map_scaleB.com
#

set pdbfile = noH.pdb
set mtzfile = refmacout.mtz
set outfile = refme_minusol.mtz
set FP = "auto"
set FreeR_flag = "auto"
set reso = "auto"
set FC_subtract = FC_ALL,PHIC_ALL
set FC_addback = FC,PHIC

set squish_scale = -0.5


# B factor for enlarging Shannon voxels
set ShannonB = 0

# B factor to use when passing sharp probability mask through fft filter
# setting to zero means no filter, very small value means just back-and-forth fft
set filterB = 1e-6
# number of times to pass mask through fft filter
set recycles = 5

# select buffer zone around significant difference peaks
set radius = 0.5
set bevel = 0.5

# probabilities considered definitely there vs unlikely
# protein map will be range-stretched to make this range equal 0-1
set protein_highprob = 0.9
set protein_lowprob = 0.1

# probabilities considered definitely there vs unlikely
#set highprob = 0.87
#set lowprob = 0.145
set highprob = 0.9
set lowprob = 0.1

#set squish_frac = 1.0
#set probfilter_frac = 1.0

set tempfile = ./tempfile
set debug = 0

set logfile = details.log

set path = ( $path `dirname $0` )

# read the command line to update variables and other settings
foreach Arg ( $* )
    set arg = `echo $Arg | awk '{print tolower($0)}'`
    set assign = `echo $Arg | awk '{print ( /=/ )}'`
    set Key = `echo $Arg | awk -F "=" '{print $1}'`
    set Val = `echo $Arg | awk '{print substr($0,index($0,"=")+1)}'`
    set Csv = `echo $Val | awk 'BEGIN{RS=","} {print}'`
    set key = `echo $Key | awk '{print tolower($1)}'`
    set num = `echo $Val | awk '{print $1+0}'`
    set int = `echo $Val | awk '{print int($1+0)}'`

    if( $assign ) then
      # re-set any existing variables
      set test = `set | awk -F "\t" '{print $1}' | egrep "^${Key}"'$' | wc -l`
      if ( $test ) then
          set $Key = $Val
          echo "$Key = $Val"
          continue
      endif
      # synonyms
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) set pdbfile = "$Arg"
      if("$Arg" =~ *.mtz ) set mtzfile = "$Arg"
    endif
    if("$Key" == "debug") set debug = "1"
end
# shorthand for temp stuff
set t = ${tempfile}

if(! -e "$pdbfile") then
    set BAD = "need pdb file: $pdbfile"
    goto exit
endif
if(! -e "$mtzfile") then
    set BAD = "need refmac output mtz file: $mtzfile"
    goto exit
endif

foreach dependency ( map_scaleB.com map_func.com Tc_maskify.com addup_maps_runme.com )
   echo -n "using: "
   which $dependency
   if( $status ) then
       set BAD = "need $dependency in "'$'"path"
       goto exit
   endif
end

echo header | mtzdump hklin $mtzfile |\
   tee ${t}mtzdump.txt |\
   awk '/^ H K L /{for(i=1;i<=NF;++i)label[i]=$i}\
        /^ H H H /{for(i=1;i<=NF;++i)print $i,label[i]}' |\
   cat >! ${t}labels.txt
set SGnum = `awk '/Space group =/{print $NF+0}' ${t}mtzdump.txt | tail -n 1`
set symops = `awk -v n=$SGnum '$1==n{print $3;exit}' ${CLIBD}/symop.lib`
set CELL = `awk '/Cell Dimensions/{getline;getline;print $1+0,$2+0,$3+0,$4+0,$5+0,$6+0;exit}' ${t}mtzdump.txt`
echo $CELL |\
awk 'NF==6{DTR=atan2(1,1)/45; A=cos(DTR*$4); B=cos(DTR*$5); G=cos(DTR*$6); \
 skew = 1 + 2*A*B*G - A*A - B*B - G*G ; if(skew < 0) skew = -skew;\
 printf "%.3f\n", $1*$2*$3*sqrt(skew)}' |\
cat >! ${t}volume
set CELLvolume = `cat ${t}volume`
rm -f ${t}volume

if( "$reso" == "auto" ) then
   set reso = `awk '/Resolution Range/{getline;getline;print $6}' ${t}mtzdump.txt`
   echo "resolution in $mtzfile is $reso"
endif
if( "$FP" == "auto" ) then
   cat  ${t}labels.txt |\
   egrep -v FC |\
   awk '$1=="F"{print $2}' >! ${t}Fs.txt
   set FP = `head -n 1 ${t}Fs.txt`
   if( "$FP" == "" ) then
      set BAD = "no Fobs in $mtzfile"
      goto exit
   endif
   set SIGFP = `grep "Q $FP" ${t}labels.txt | awk '{print $2;exit}'`
   if( "$SIGFP" == "" ) then
      set SIGFP = `egrep "^Q " ${t}labels.txt | awk '{print $2;exit}'`
   endif
   if( "$SIGFP" == "" ) then
      set BAD = "no sigFobs in $mtzfile"
      goto exit
   endif
   echo "selecting Fobs = $FP $SIGFP"
endif
if( "$FreeR_flag" == "auto" ) then
   cat  ${t}labels.txt |\
   awk '$1=="I" && tolower($2)~/free/{print $2}' >! ${t}free.txt
   set FreeR_flag = `head -n 1 ${t}free.txt`
   if( "$FreeR_flag" == "" ) then
      set FreeR_flag = `egrep "^I " ${t}labels.txt | awk '{print $2;exit}'`
   endif
   if( "$FreeR_flag" == "" ) then
      set BAD = "no FreeR_flag in $mtzfile"
      goto exit
   endif
   echo "selecting FreeR_flag = $FreeR_flag"
endif
set hires = `echo $reso | awk '{print (($1**3)*0.8)**0.3333}'`

foreach FC ( $FC_subtract $FC_addback )
  echo $FC |\
   awk -F "," '{print $1,$2}' |\
   cat - ${t}labels.txt |\
   awk 'NR==1{F=$1;P=$2;next}\
     {++seen[$2]}\
     END{print seen[F],seen[P]}' >! ${t}test.txt
  set test = `cat ${t}test.txt`
  if( "$test" != "1 1" ) then
    set BAD = "no $FC in $mtzfile"
    goto exit
  endif
end


echo "making $FP - $FC_subtract difference map"
set FC = `echo $FC_subtract | awk -F "," '{print $1}'`
set PHIC = `echo $FC_subtract | awk -F "," '{print $2}'`
fft hklin $mtzfile mapout ffted.map << EOF >! $logfile
labin F1=$FP F2=$FC PHI=$PHIC
reso $hires
EOF
mapmask mapin ffted.map mapout fofc.map << EOF >> $logfile
xyzlim asu
axis X Y Z
EOF
fft hklin $mtzfile mapout ffted.map << EOF >! $logfile
labin F1=$FP PHI=$PHIC
reso $hires
EOF
mapmask mapin ffted.map mapout Fobs.map << EOF >> $logfile
xyzlim asu
axis X Y Z
EOF

echo | mapdump mapin fofc.map >! ${t}mapdump.txt


foreach C ( C C_ALL C_ALL_LS )
fft hklin $mtzfile mapout ffted.map << EOF >! $logfile
labin F1=F$C PHI=PHI$C
reso $hires
EOF
mapmask mapin ffted.map mapout F${C}.map << EOF >> $logfile
xyzlim asu
axis X Y Z
EOF
end


# ' "probability that something is there" '
# ' prob(rho) = sign(rho)*pow(erf(abs(rho/sigma(rho))/sqrt(2)),V/2/d^3) '




# number of Shannon voxels expected
set exponent = `echo $CELLvolume $symops $ShannonB | awk '{V=$1;n=$2;B=$3+$4;pi=atan2(1,1)*4;print V/n/2/(sqrt((B+9.484)*log(2))/pi)**3}'`
echo "estimating $exponent Shannon voxels in ASU"

echo "taking absolute value of Fo-Fc scaled to sigma=0.7071"
mapmask mapin fofc.map mapout fofc_sigsqrt2.map << EOF >> $logfile
scale sigma 0.70701 0
EOF
map_func.com -func abs  fofc_sigsqrt2.map -output abs.map >> $logfile
echo "erfpow"
map_func.com -func erfpow -param $exponent  abs.map -output erfpow.map >> $logfile

echo | mapdump mapin erfpow.map >! ${t}mapdump.txt
set xyzlim = `awk '/stop points on col/{gsub(":"," ");print $9,$10,$11,$12,$13,$14; exit}' ${t}mapdump.txt`
set paddedlim = `echo $xyzlim | awk '{print $1-5,$2+5,$3-5,$4+5,$5-5,$6+5}'`

mapmask mapin erfpow.map mapout padded.map << EOF >> $logfile
xyzlim $paddedlim
EOF
echo "feather"
map_func.com -func feather  padded.map -output padfeather.map >> $logfile
mapmask mapin padfeather.map mapout feathered.map << EOF >> $logfile
xyzlim asu
EOF

foreach itr ( `seq 1 $recycles` )
    echo "recombining feathered mask with ffted version of itself: $itr of $recycles"
    map_scaleB.com feathered.map B=$filterB output=filtered.map clip=0 >> $logfile

    map_func.com -func max feathered.map filtered.map -output maxmix.map >> $logfile

    cp maxmix.map feathered.map
end

# compute scale and B needed to generate specified radius and feather
echo "selecting area within $radius A of selected voxels, beveled by $bevel A"
echo $radius $bevel |\
awk '{radius=$1;feather=$2;pi=4*atan2(1,1);\
  rat = ((radius+feather/2)/radius)**2;\
  s = (0.977322 - 0.159175*rat)/(rat - 1);\
  B = -4*pi**2*radius**2/log(-log(0.5)*10**-s);\
  print -(10**s),B}' >! ${t}sB.txt
set bevel_scale = `awk '{print $1}' ${t}sB.txt`
set bevel_B     = `awk '{print $2}' ${t}sB.txt`

map_scaleB.com maxmix.map B=$bevel_B scale=$bevel_scale output=expme.map
map_func.com -func exp expme.map -output notmask.map >> $logfile
map_func.com -func negate -param 1 notmask.map -output beveled.map >> $logfile

echo "stretching range of significance mask: $lowprob : $highprob -> 0 : 1"
rm -f temp.map stretched.map loclip.map clipped.map > /dev/null
set mult = `echo $highprob $lowprob | awk '{print 1./($1-$2)}'`
set offs = `echo $highprob $lowprob | awk '{print -$2/($1-$2)}'`
# started out using mult=1.5 adn offs=-0.2
map_func.com -func mult -param $mult beveled.map -output temp.map >> $logfile
map_func.com -func add -param $offs temp.map -output stretched.map >> $logfile
map_func.com -func max -param 0 stretched.map -output loclip.map >> $logfile
map_func.com -func min -param 1 loclip.map -output clipped.map >> $logfile

cp clipped.map significance.map
# significance map reflects probability that the fo-fc feature is noise

# now we need to mask out the protein
echo "counting conformers"
cat $pdbfile |\
awk '! /^ATOM|^HETAT/{next}\
 {conf=substr($0,17,1);occ=substr($0,55,6);\
  if(conf==" ")next;\
  ++count[conf];sum[conf]+=occ}\
 ! order[conf]{++c;order[conf]=c}\
 END{for(conf in count){\
   print order[conf],conf,sum[conf]/count[conf],count[conf]}}' |\
sort -g |\
tee ${t}conf_occ.txt 


set masks = ""
foreach cnfnum ( `awk '{print $1}' ${t}conf_occ.txt` )
set cnfchar = `egrep "^${cnfnum} " ${t}conf_occ.txt | awk '{print $2}'`

echo "masking around conf $cnfchar "

egrep "^${cnfnum} " ${t}conf_occ.txt |\
cat - $pdbfile |\
awk 'NR==1{selected=$2;avgocc=$3+1e-6;next}\
  ! /^ATOM|^HETAT/{print;next}\
    {conf=substr($0,17,1);occ=substr($0,55,6);\
     pre=substr($0,1,54);post=substr($0,61)}\
 conf==" "{print;next}\
 conf!=selected{next}\
    {occ=occ/avgocc}\
  occ>1{occ=1} occ<0{occ=0}\
  {printf("%s%6.2f%s\n",pre,occ,post)}' |\
cat >! maskme_${cnfnum}.pdb
Tc_maskify.com fofc.map maskme_${cnfnum}.pdb outfile=mask_${cnfnum}.map >! ${t}maskify_${cnfnum}.log
set masks = ( $masks mask_${cnfnum}.map )
end

echo "averaging $#masks masks..."
addup_maps_runme.com $masks >&! ${t}addup.log

set scale = `echo $#masks | awk '{print 1/$1}'`

mapmask mapin sum.map mapout avg.map << EOF >> $logfile
scale factor $scale
axis X Y Z
EOF
cp avg.map not_protein_rough.map
#map_func.com -func negate not_protein.map -outfile squish_solvent_rough.map

echo "stretching range of protein mask: $protein_lowprob : $protein_highprob -> 0 : 1"
rm -f temp.map stretched.map loclip.map clipped.map > /dev/null
set mult = `echo $protein_highprob $protein_lowprob | awk '{print 1./($1-$2)}'`
set offs = `echo $protein_highprob $protein_lowprob | awk '{print -$2/($1-$2)}'`
map_func.com -func mult -param $mult not_protein_rough.map -output temp.map >> $logfile
map_func.com -func add -param $offs temp.map -output stretched.map >> $logfile
map_func.com -func max -param 0 stretched.map -output loclip.map >> $logfile
map_func.com -func min -param 1 loclip.map -output clipped.map >> $logfile

cp clipped.map not_protein.map
# not protein is a map that is 0  around the protein and 1 far from it



echo maps mult |\
mapmask mapin1 not_protein.map mapin2 significance.map mapout significant_solvent_mask.map >> $logfile

map_func.com -func negate significant_solvent_mask.map -outfile squish_solvent_mask.map

map_func.com -func multiply fofc.map squish_solvent_mask.map -outfile fofc_squish_sol.map

map_func.com -func multiply fofc.map significant_solvent_mask.map -outfile fofc_sig_sol.map

gemmi map2sf -v --dmin=$reso fofc_squish_sol.map gemmi.mtz dF PHIdF 
echo labin file 1 all | cad hklin1 gemmi.mtz hklout fofc_squish_sol.mtz  >> $logfile
gemmi map2sf -v --dmin=$reso fofc_sig_sol.map gemmi.mtz sF PHIsF 
cad hklin1 gemmi.mtz hklout fofc_sig_sol.mtz << EOF >> $logfile
labin file 1 all
scale file 1 $squish_scale
EOF

echo "adding $FC_addback to modified fofc to get new FP"
set FC = `echo $FC_addback | awk -F "," '{print $1}'`
set PHIC = `echo $FC_addback | awk -F "," '{print $2}'`
cad hklin1 $mtzfile hklout cadded.mtz << EOF >> $logfile
labin file 1 E1=$FP E2=$SIGFP E3=$FC E4=$PHIC
labou file 1 E1=FP E2=SIGFP E3=FC E4=PHIC
EOF

rm -f new.mtz
sftools << EOF >> $logfile
read fofc_squish_sol.mtz
read fofc_sig_sol.mtz
read cadded.mtz col FP SIGFP FC PHIC
calc ( COL Fms PHI ) = ( COL FC PHIC ) ( COL dF PHIdF ) +
calc ( COL Fms PHI ) = ( COL Fms PHI ) ( COL sF PHIsF ) +
calc F col FP = col Fms
absent col FP if col SIGFP absent
select col SIGFP = PRESENT
purge nodata yes
select all
write new.mtz col FP SIGFP PHI
quit
y
EOF
echo "cadding in $FreeR_flag for final $outfile"
cad hklin1 new.mtz hklin2 $mtzfile hklout $outfile << EOF >> $logfile
labin file 1 E1=FP E2=SIGFP
labin file 2 E1=$FreeR_flag
EOF

fft hklin new.mtz mapout ffted.map << EOF >! $logfile
labin F1=FP PHI=PHI
reso $hires
EOF
mapmask mapin ffted.map mapout new_FP.map << EOF >> $logfile
xyzlim asu
axis X Y Z
EOF


exit:

if("$tempfile" == "") set  tempfile = "./"
set tempdir = `dirname $tempfile`
if(! $debug && ! ( "$tempdir" == "." || "$tempdir" == "" ) ) then
    rm -f ${tempfile}*
endif

if($?BAD) then
   echo "ERROR: $BAD"
   exit 9
endif
exit






rm -f new.mtz
sftools << EOF >> $logfile
read fofc_resid.mtz
read refmacout.mtz col FC PHIC
calc ( COL FP PHI ) = ( COL FC PHIC ) ( COL dF PHIdF ) +
absent col FP if col SIGFP absent
select col SIGFP = PRESENT
purge nodata yes
select all
write new.mtz col FP PHI
quit
y
EOF
fft hklin new.mtz mapout ffted.map << EOF >! $logfile
labin F1=FP PHI=PHI
EOF
mapmask mapin ffted.map mapout new_FP.map << EOF >> $logfile
xyzlim asu
axis X Y Z
EOF

