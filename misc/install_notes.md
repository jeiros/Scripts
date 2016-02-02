|_   _| |          | ___ \         (_)         
  | | | |__   ___  | |_/ / __ _ ___ _  ___ ___ 
  | | | '_ \ / _ \ | ___ \/ _` / __| |/ __/ __|
  | | | | | |  __/ | |_/ / (_| \__ \ | (__\__ \
  \_/ |_| |_|\___| \____/ \__,_|___/_|\___|___/
                                                 
The BASH interpreter in OS X loads up .bash_profile, not .bashrc!! Care with that.
The home (~) directory in OS X is /Users/yourusername/, unlike Linux which is /home/yourusername/
Read **all** of this before starting to do anything                                             

# Registering your Mac
Either ICT does it for you and they install their build, or you do it yourself (and follow this guide)

    Follow these steps in ORDER:
    1)Email Nick Davies (n.davies@imperial.ac.uk) and he will give you a tag asset number
        1.1) Once you got it, plug your iMac to the ethernet and follow the self-registration process
    2) Register your machine through ICT's website https://imperial.service-now.com/ict/request.do
    3) I requested "Add or change the registration details of a computer"
        3.1) You'll need your MAC adress for that
        3.2) Give it a name (ch-whatever)
        3.3) Your machine's host name will then be ch-whatever.ch.ic.ac.uk
            So you can ssh from another machine by doing: ssh yourusername@ch-whatever.ch.ic.ac.uk

# Antivirus
Install the symantec antivirus protection software from the software shop
https://www.imperial.ac.uk/ict/services/software/shop/index.asp
you'll receive an e-mail with the instructions to follow

# Microsoft Office 365
Follow these instructions:
http://www.imperial.ac.uk/admin-services/ict/shop/software/microsoft-office-365/install-office-365/mac/

    *Care* --> Point 2 is important. DON'T use your e-mail address, it's username@ic.ac.uk
            otherwise Microsoft Office won't recognize you as an Imperial College member.


# Installing Anaconda Python distribution

Go to this website https://www.continuum.io/downloads
and install the PYTHON 3.5 OSX version. Follow the instructions.

http://conda.pydata.org/docs/intro.html 




# Installing AmberTools
Follow Jason Swail's guide --> http://jswails.wikidot.com/mac-os-x

Basically: 
    1) Instal Xcode
    2) Enable command-line tools
    3) Download macports
    4) Use macports to install the compilers for  AmberTools to work


The gcc version 4.8 (known bug in here: https://trac.macports.org/ticket/48471) failed for me so I installed the 4.9
so the commands were:

    sudo port install gcc49
    sudo port install mpich-gcc49
    sudo port select --set gcc mp-gcc49
    sudo port select --set mpi mpich-gcc49-fortran

Now onto actually installing the AmberTools:

    1) Download them from here --> http://ambermd.org/AmberTools15-get.html
    2) Follow the instructions in page 23 of the manual
        2.1) Everything worked for the most simple installation (./configure gnu)

                All the tests PASSED
                    1382 file comparisons passed
                       0 file comparisons failed
                       0 tests experienced errors
                Test log file saved as /usr/local/amber15/logs/test_at_serial/2015-10-28_15-26-42.log
                No test diffs to save!

# ~~~~~Installing VMD
Here are all of the downloadables --> http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD

I chose the latest version (1.9.2) and the MacOS X OpenGL (32-bit Intel x86) link
A .dmg file is downloaded. Open it and then drag the VMD icon to the Applications folder. Done!

It won't let you just open it yet: Apple blocks applications from non-verified distributors (a.k.a non-Apple)
so you have to go to System Preferences > Security & Privacy and click on the "General" tab. There, click on the
lower-left lock icon, and then select the "Anywhere" option under the "Allow apps downloaded from:" section. 
Now you can use vmd.  

How to use vmd from the command line: 
    1) Remove the space from the name of the Application (i.e. VMD1.9.2)
    2) Add this to your .bash_profile
        vmdappdir='/Applications/VMD1.9.2.app/Contents'
        alias vmd='"$vmdappdir/Resources/VMD.app/Contents/MacOS/VMD" $*'

NOTE: VMD only has 32-bit version for Mac OS X. This sucks because you'll only be able to open
trajectories that are half the size of your RAM (i.e. < 4GB.)

# ~~~~~~Installing UCSF Chimera
Here are all the downloadables --> https://www.cgl.ucsf.edu/chimera/download.html
Same as VMD. Click on .dmg file, drag Chimera icon to the Applications folder.

The binary is in /Applications/Chimera.app/Contents/MacOS/chimera
    either create an alias in your .bash_profile to it or do a symbolic link, either option works.
    


 _____ _   _                     _          __  __ 
|  _  | | | |                   | |        / _|/ _|
| | | | |_| |__   ___ _ __   ___| |_ _   _| |_| |_ 
| | | | __| '_ \ / _ \ '__| / __| __| | | |  _|  _|
\ \_/ / |_| | | |  __/ |    \__ \ |_| |_| | | | |  
 \___/ \__|_| |_|\___|_|    |___/\__|\__,_|_| |_|  
                                                   
                                                   

# ~~~~~Generic Note 1

The new OS X El Capitan includes a new "rootless" mode that makes certain system directories
read-only even for admins. This is a bit annoying, because one of these is /usr/ (and thus everything inside /usr/ can't be touched
the only subdir that is excluded is /usr/local).
This can be disabled following these instructions but I haven't tried it myself.
http://apple.stackexchange.com/questions/196224/unix-ln-s-command-not-permitted-in-osx-el-capitan-beta3 

# ~~~~~Generic Note 2

Another "feature" that has been added is the creation of bash sessions. This creates a .bash_sessions/ directory
in your home that keeps storing files. It can be annoying since these are read
every time you launch a terminal and in the end it can make it slower to start. This can be disabled if you create
a .bash_sessions_disable file in your home directory, and delete de .bash_sessions/ directory with rm -rf ~/.bash_sessions

# ~~~~~Generic Note 3

Not OS X specific but I find it useful to be able to ssh withouth entering passwords all the time. 
This is handy when you are working with several machines (yours, the HPC, your laptop, other local machines with GPUs...)
This link explains how to do it. http://www.linuxproblem.org/art_9.html
If properly set up, you'll be able to ssh and scp across diferent machines without having to enter the password every time.
Super useful when you add this to scripts



Written by: Juan Eiros
j.eiros-zamora14@imperial.ac.uk
I will keep adding stuff that I find useful

