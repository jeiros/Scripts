# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# enable color support of ls and also add handy aliases
if [ -x /usr/bin/dircolors ]; then
    test -r ~/.dircolors && eval "$(dircolors -b ~/.dircolors)" || eval "$(dircolors -b)"
    alias ls='ls --color=auto'
    #alias dir='dir --color=auto'
    #alias vdir='vdir --color=auto'
 
    alias grep='grep --color=auto'
    alias fgrep='fgrep --color=auto'
    alias egrep='egrep --color=auto'
fi

# User specific aliases and functions
alias cpptraj2='/home/je714/cpptraj-master/bin/cpptraj'
alias c='clear'
alias vi='vim'
alias l='ls -XGlh'
alias sshHPC='ssh login.cx1.hpc.ic.ac.uk'
alias dropbox='/home/je714/dropbox.py'
alias myvmd='/usr/local/bin/vmd -size 1920 1200 -e /home/je714/StateFile'
alias ..='cd ../'
alias ....='cd ../../'
alias phospho='cd /home/je714/Troponin/IAN_Troponin/completehowarthcut/phospho/hmr_runs'
alias condapy27='source activate py27'
alias condaroot='source activate root'
export msmbuilder='/home/je714/anaconda3/lib/python3.4/site-packages/msmbuilder-3.3.0-py3.4-linux-x86_64.egg/msmbuilder'


export SCRIPTS=/home/je714/Scripts/
export AMBERHOME=/usr/local/amber15/
export PATH="$AMBERHOME/bin:/usr/local/amber/bin/:/usr/local/bin/:$SCRIPTS:$SCRIPTS/MSManalysis/:/usr/hs/bin:$PATH"
source $AMBERHOME/amber.sh
export LD_LIBRARY_PATH="/home/je714/pytraj/cpptraj/lib/:${LD_LIBRARY_PATH}"
export CX1="login.cx1.hpc.ic.ac.uk"
PS1='$(whoami):$(pwd)$ '

# added by Anaconda3 2.1.0 installer
export PATH="/home/je714/anaconda3/bin:$PATH"
