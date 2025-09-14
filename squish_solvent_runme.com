#! /bin/tcsh -f
#
# remove bulk solvent from Fobs
#
# requires Tc_maskify.com map_func.com 
#

set pdbfile = refmacout.pdb
set mtzfile = refmacout.mtz
set outfile = refme_minusol.mtz
set FP = FP
set mask_scale = 50
set reso = "auto"

set tempfile = ./tempfile
set debug = 0


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



set hires = `echo $reso | awk '{print $1*0.9}'`

fft hklin $mtzfile mapout ffted.map << EOF
labin F1=FP F2=FC_ALL PHI=PHIC_ALL
reso $hires
EOF
echo xyzlim asu | mapmask mapin ffted.map mapout fofc.map 

cat $pdbfile |\
awk -v scale=$mask_scale '! /^ATOM|^HETAT/{print;next}\
   {occ=substr($0,55,6);pre=substr($0,1,54);post=substr($0,61);\
    newocc=occ*scale;}\
   newocc>1{newocc=1}\
   {printf("%s%6.2f%s\n",pre,newocc,post)}' |\
cat >! maskme.pdb
Tc_maskify.com fofc.map maskme.pdb

map_func.com -func negate mask.map -outfile solsquash_mask.map

map_func.com -func multiply fofc.map solsquash_mask.map -outfile fofc_sqsh.map

gemmi map2sf -v --dmin=$reso fofc_sqsh.map gemmi.mtz dF PHIdF 
echo labin file 1 all | cad hklin1 gemmi.mtz hklout fofc_sqsh.mtz
echo labin file 1 all | cad hklin1 $mtzfile hklout cadded.mtz

rm -f new.mtz
sftools << EOF
read fofc_sqsh.mtz
read cadded.mtz col FP SIGFP FC PHIC
calc ( COL FP PHI ) = ( COL FC PHIC ) ( COL dF PHIdF ) +
absent col FP if col SIGFP absent
select col SIGFP = PRESENT
purge nodata yes
select all
write new.mtz col FP SIGFP
quit
y
EOF
cad hklin1 new.mtz hklin2 cadded.mtz hklout $outfile << EOF
labin file 1 all
labin file 2 E1=FreeR_flag
EOF

