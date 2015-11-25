
# Setting PATH for Python 3.4
# The orginal version is saved in .bash_profile.pysave


# Source global definitions
if [ -f /etc/bashrc ]; then
    . /etc/bashrc
fi




# Aliases
alias top='htop'
alias ..='cd ..'
alias rm='rm -i'
alias c='clear'
alias vi='vim'
alias l='ls -Glvh'
alias rmtgz='ls | grep -v *.tgz | xargs rm'  # Removes all the files in the directory that are not .tgz files
alias tree='tree -h'
alias du='du -h'

# SSH aliases
alias sshHPC='ssh je714@login.cx1.hpc.ic.ac.uk'
alias sshGuanine='ssh je714@guanine.ch.ic.ac.uk'
alias sshBertha='ssh je714@bertha.ch.ic.ac.uk'
alias sshTitan='ssh je714@ch-ch-igould-titanx1.ch.ic.ac.uk'

# Conda environments
alias condapy27='source activate py27'
alias condaroot='source activate root'

# Quick cd to working directories
alias phospho='cd /Users/je714/Troponin/IAN_Troponin/completehowarthcut/phospho/hmr_runs'
alias salted='cd /Users/je714/Troponin/IAN_Troponin/completehowarthcut/salted'

# VMD and Chimera aliases
vmdappdir='/Applications/VMD1.9.2.app/Contents'
alias vmd='"$vmdappdir/Resources/VMD.app/Contents/MacOS/VMD" $*'
alias myvmd='vmd -e ~/Scripts/StateFile -size 1920 1080'
alias chimera='/Applications/Chimera.app/Contents/MacOS/chimera'

# Connections
export CX1="login.cx1.hpc.ic.ac.uk"
export TITAN='ch-ch-igould-titanx1.ch.ic.ac.uk'

# Command line prompt
PS1='$(whoami)@$(hostname):$(pwd)$ '


# MacPorts Installer addition on 2015-07-27_at_16:06:16: adding an appropriate PATH variable for use with MacPorts.
export PATH="/opt/local/bin:/opt/local/sbin:/opt/bin:$PATH"

# added by Anaconda3 2.3.0 installer
export PATH="/Users/je714/anaconda/bin:$PATH"

# MacPorts Installer addition on 2015-10-28_at_10:10:59: adding an appropriate PATH variable for use with MacPorts.
export PATH="/opt/local/bin:/opt/local/sbin:$PATH"
# Finished adapting your PATH environment variable for use with MacPorts.

# Add the Scripts repository (git@github.com:jeiros/Scripts.git) to the PATH 
export PATH="/Users/je714/Scripts:/Users/je714/Scripts/AnalysisMDTraj:/Users/je714/Scripts/AnalysisCpptraj:$PATH"

# AMBER installation
export AMBERHOME="/usr/local/amber15"
test -f /usr/local/amber15/amber.sh && source /usr/local/amber15/amber.sh


export SCRIPTS='/Users/je714/Scripts'
