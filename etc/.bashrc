# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions
module load gcc
module load slurm 
umask 0007

# Marty's additions
module load julia


# Marty's additions
export JULIA_APP_BOOTSTRAP=/home/mhimmelstein/CodeandSampleData/conf/taylor_work.conf

