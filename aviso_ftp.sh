#!/bin/sh

HOST='ftp://ftp.aviso.altimetry.fr'
USER='uw_hennon' # username required by this ftp site
PASSWD='ovsd56dfsd' # HYCOM prefers an email address, but really you can enter anything
DIR='/global/delayed-time/grids/madt/all-sat-merged/uv'

echo "What year are you downloading?"
read year

 suffix="20140106.nc.gz"
#suffix="20140704.nc.gz"

 startmonth=1
 endmonth=12
 
# download all the things!
for month in $(seq $startmonth $endmonth) 
do
	for day in 1 9 17 25
	do
		dirpath="${HOST}${DIR}/${year}/" # directory of the file
		fname=$(printf 'dt_global_allsat_madt_uv_%s%02d%02d_%s' $year $month $day $suffix)
		fullpath="${dirpath}${fname}" # full path of file
		echo "$fullpath" # sanity check: print the name of the file being downloaded
		wget --connect-timeout=21600 -nv -nc -c --user=$USER --password=$PASSWD $fullpath
	done
done

