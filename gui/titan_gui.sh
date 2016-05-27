#! /bin/sh

# Script to start the Titan2D GUI.  Enter: ./run_titan

# The location of the Titan2D application binary must be included in the PATH environment variable.  
# (Note: On VHub, the location of the Titan application binary is set via HUBzero invoke processing when the Titan2D tool GUI is launched).  Examples:
# VHub:
# use titan2d-v4.0.0
# rush:
# module load titan
# otherwise:
# Add the path to the Titan2D application binary to the PATH environment variable.
# Example, in your .bashrc file:
# export PATH=$PATH:"/projects/academic/gmfg/renettej/titan2d_github/titan2d/bin"

# Note: on VHub, the TOOLDIR environment variable is set via HUBzero invoke processing when the Titan2D tool GUI is launched.  Need TOOLDIR defined to create KML files.  Setting TOOLDIR to the parent directory.
export TOOLDIR="$(dirname "$(pwd)")"
echo $TOOLDIR

# Note: The Titan2D GUI's Job Submission Tab's Hub-Submit Run Style option is available on VHub only.
# Note: On VHub, only 1 node and 1 cpu are allowed per job.
# Default, assumes running Titan2d GUI on VHub
export E_VHUB="false"

# Specify the start directory for the Titan2D GUI's directory selectors.
# Default, $Home
# Example:
# export E_INPUTDIR="/projects/academic/gmfg/renettej/titan_java_gui"

fullpath="$(readlink -f $0)"

export PATH="$(echo "$fullpath" | sed "s?/bin/titan_gui.sh?/bin?"):$PATH"

libpath="$(echo "$fullpath" | sed "s?/bin/titan_gui.sh?/lib/titan_java_gui?")"

# Help files are installed to $(docdir)
helppath="$(echo "$fullpath" | sed "s?/bin/titan_gui.sh?/share/doc/titan2d?")"

export CLASSPATH=${CLASSPATH}:$libpath/titan_gui.jar:$libpath/jh.jar:$libpath/derby.jar:$helppath

java titan.gui.Titan
