
# Setting PATH for Python 3.4
# The orginal version is saved in .bash_profile.pysave
alias c='clear'
alias sshImp='ssh je714@155.198.225.40'
alias vi='vim'
alias l='ls -Glvh'
alias sshHPC='ssh je714@login.cx1.hpc.ic.ac.uk'
alias myvmd='/Applications/VMD1.9.2.app/Contents/MacOS/startup.command -size 1200 750 -e /Users/Juan/Desktop/StateFile_macos.vmd'
export IMP="155.198.225.40"
PATH="/Library/Frameworks/Python.framework/Versions/3.4/bin:${PATH}"
export PATH
PS1='$(whoami)@$(hostname):$(pwd)$ '
export SCRIPTS='/Users/Juan/Desktop/Scripts'
export HD='/Volumes/My\ Passport/'
alias sshMarce='ssh Juan@188.166.104.176'

##
# Your previous /Users/Juan/.bash_profile file was backed up as /Users/Juan/.bash_profile.macports-saved_2014-12-12_at_15:53:17
##

# MacPorts Installer addition on 2014-12-12_at_15:53:17: adding an appropriate PATH variable for use with MacPorts.
export PATH="/opt/local/bin:/opt/local/sbin:/opt/bin:$SCRIPTS:$PATH"
# Finished adapting your PATH environment variable for use with MacPorts.

export CLICOLOR=1
export LSCOLORS=ExFxCxDxBxegedabagacad

# added by Anaconda3 2.2.0 installer
export PATH="/Users/Juan/anaconda/bin:$PATH"
export AMBERHOME="/usr/local/amber15"
source $AMBERHOME/amber.sh
##
# Your previous /Users/Juan/.bash_profile file was backed up as /Users/Juan/.bash_profile.macports-saved_2015-07-27_at_16:06:16
##

# MacPorts Installer addition on 2015-07-27_at_16:06:16: adding an appropriate PATH variable for use with MacPorts.
export PATH="/opt/local/bin:/opt/local/sbin:$PATH"
# Finished adapting your PATH environment variable for use with MacPorts.

