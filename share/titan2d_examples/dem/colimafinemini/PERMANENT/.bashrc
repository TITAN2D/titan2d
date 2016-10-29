test -r ~/.alias && . ~/.alias
PS1='GRASS 7.0.4 (colimafinemini):\w > '
grass_prompt() {
	LOCATION="`g.gisenv get=GISDBASE,LOCATION_NAME,MAPSET separator='/'`"
	if test -d "$LOCATION/grid3/G3D_MASK" && test -f "$LOCATION/cell/MASK" ; then
		echo [2D and 3D raster MASKs present]
	elif test -f "$LOCATION/cell/MASK" ; then
		echo [Raster MASK present]
	elif test -d "$LOCATION/grid3/G3D_MASK" ; then
		echo [3D raster MASK present]
	fi
}
PROMPT_COMMAND=grass_prompt
export GRASS_GNUPLOT="gnuplot -persist"
export GRASS_PROJSHARE=/usr/share/proj
export GRASS_ADDON_BASE=/home/mikola/.grass7/addons
export GRASS_HTML_BROWSER=xdg-open
export GRASS_PYTHON=python
export GRASS_PAGER=pager
export GRASS_VERSION=7.0.4
export PATH="/usr/lib/grass70/bin:/usr/lib/grass70/scripts:/home/mikola/.grass7/addons/bin:/home/mikola/.grass7/addons/scripts:/var/www/akrr/bin:/usr/local/ParaView-4.4.0-Qt4-Linux-64bit/bin:/usr/local/openmpi-1.10.1-icc/bin:/opt/intel/compilers_and_libraries_2016.1.150/linux/bin/intel64:/opt/intel/compilers_and_libraries_2016.1.150/linux/mpi/intel64/bin:/opt/intel/debugger_2016/gdb/intel64_mic/bin:/var/www/akrr/bin:/usr/local/ParaView-4.4.0-Qt4-Linux-64bit/bin:/usr/local/openmpi-1.10.1-icc/bin:/opt/intel/compilers_and_libraries_2016.1.150/linux/bin/intel64:/opt/intel/compilers_and_libraries_2016.1.150/linux/mpi/intel64/bin:/opt/intel/debugger_2016/gdb/intel64_mic/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games"
export HOME="/home/mikola"
