#!/bin/sh

# First, mount to nightraid (This protocol works best/only? from
# Macs).  This requires you to have a DSA/RSA Key Set up already. 

echo "Password to mount to nightraid?" 
stty -echo
read PASS
stty echo

echo "     MSG: Mount to nightraid..."
mkdir /Volumes/nightraid_dbernat
mount_afp afp://dbernat:$PASS@nightraid.astro.cornell.edu/dbernat /Volumes/nightraid_dbernat

logname=`date "+%y%m%d-%H%M%S"` 
logname=$logname".dwblog"
echo "     MSG: Logging IDL output in file: " $logname

#
# Configure the local machine to use IDL and IDL Code on Nightraid
#

echo "     MSG: Adding to PATH the idl directory on nightraid..."
export LM_LICENSE_FILE=1700@idl-server.astro.cornell.edu
export IDL_DIR=/Volumes/nightraid_dbernat/itt/idl704
export PATH=$PATH:$IDL_DIR/bin/

echo "     MSG: Adding /Volumes/nightraid_dbernat/Research/idl-code to local IDL_PATH"
export IDL_PATH=$IDL_PATH:+/Volumes/nightraid_dbernat/Research/:'<IDL_DEFAULT>'

#
# Launch Cron
#

echo "     MSG: Launching IDL Cron..."
nohup idl -e "cron_emp_f_dist, WORKINGDIR='/Volumes/nightraid_dbernat/Research/redux/'" > $logname &
#nohup idl -e "print, indgen(10)"

#
# Display output of log to screen to follow
#

# tail -n 20 -f $logname