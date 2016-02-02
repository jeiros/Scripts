# The Basics
This guide is written for not-very-experienced unix users. Sorry if some of the stuff comes across as obvious. Also, some installation procedures
are a personal choice so please feel free to change them around. 

Nevertheless, please read **ALL** of this before starting to do anything.
                                      

**Note of care for Linux users:** The home `~/` directory in OS X is `/Users/yourusername/`, unlike Linux which is `/home/yourusername/`.
The BASH interpreter in OS X loads up `~/.bash_profile`, not `~.bashrc`.

## Registering your Mac
Either ICT does it for you and they install their build, or you do it yourself (and follow this guide).

Follow these steps:

1. Email [Nick Davies](mailto:n.davies@imperial.ac.uk) and he will give you a tag asset number

2. Register your machine through [ICT's website](https://imperial.service-now.com/ict/request.do)

 * Once you got it, plug your iMac to the ethernet and follow the self-registration process

3. I requested "Add or change the registration details of a computer"

 * You'll need your MAC adress for that. [Here](http://www.iclarified.com/30929/how-to-find-your-mac-address-in-mac-os-x)'s how you find it. 
 * Give it a name (*ch-whatever*)
 * Your machine's host name will then be `ch-whatever.ch.ic.ac.uk`

Then you'll be able to ssh from another machine by doing: `ssh yourusername@ch-whatever.ch.ic.ac.uk`


## Antivirus
Install the symantec antivirus protection software from the [software shop](https://www.imperial.ac.uk/ict/services/software/shop/index.asp).
You'll receive an e-mail with the instructions to follow.

## Microsoft Office 365
Follow these [instructions](http://www.imperial.ac.uk/admin-services/ict/shop/software/microsoft-office-365/install-office-365/mac/).

**Care** :heavy_exclamation_mark:
Point 2 is important. DON'T use your e-mail address, it's *username@ic.ac.uk*. Otherwise Microsoft Office won't recognize you as an Imperial College member.


## Installing Anaconda Python distribution

Go to this [website](https://www.continuum.io/downloads)
and install the PYTHON 3.5 OSX version. Follow the instructions.

[Here](http://conda.pydata.org/docs/intro.html) is a nice introduction to what is Anaconda and what can you do with it.
[This](http://conda.pydata.org/docs/_downloads/conda-cheatsheet.pdf) cheat sheet is also very useful for quick
reference.




## Installing AmberTools
Follow Jason Swail's [guide](http://jswails.wikidot.com/mac-os-x). 

It's specially important in here to also *read* the **ENTIRE** guide before starting
to blindly type any of the commands. That'll help you prevent errors.

In essence, these are the necessary things for AmberTools to work on a Mac:

1. Instal Xcode

2. Enable command-line tools

3. Download macports

4. Use macports to install the compilers for AmberTools to work


The gcc version 4.8 (known bug in [here](https://trac.macports.org/ticket/48471) failed for me so I installed the 4.9 version. 
The commands were:
```
sudo port install gcc49
sudo port install mpich-gcc49
sudo port select --set gcc mp-gcc49
sudo port select --set mpi mpich-gcc49-fortran
```
Once you complete all of the above steps without errors, you can actually start the installation of AmberTools:
1. Download them from [here](http://ambermd.org/AmberTools15-get.html)

2. Follow the instructions in page 23 of the [manual](http://ambermd.org/doc12/Amber15.pdf).

 * Everything worked for me using the most simple installation (`./configure gnu`)

If you do `make test` in your `$AMBERHOME` directory and everything is OK, you should see the following (after a lot of output):

```
All the tests PASSED
    1382 file comparisons passed
    0 file comparisons failed
    0 tests experienced errors
    Test log file saved as /usr/local/amber15/logs/test_at_serial/2015-10-28_15-26-42.log
No test diffs to save!
```

## Installing VMD
[Here](http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD) are all of the downloadables.

I chose the latest version (1.9.2) and the MacOS X OpenGL (32-bit Intel x86) link.

A `.dmg` file is downloaded. Open it and then drag the VMD icon to the Applications folder. Done!

It won't let you just open it yet: Apple blocks applications from non-verified distributors (a.k.a non-Apple)
so you have to go to *System Preferences > Security & Privacy* and click on the *General* tab. There, click on the
lower-left lock icon, and then select the *Anywhere* option under the *Allow apps downloaded from:* section. 
Now you can use vmd.  

### How to use vmd from the command line
Using VMD from the command line is very convenient as you can directly load the topology and trajectories files, or several PDBs
very quickly. If you open VMD manually from the Dock it points every time to your home directory and you loose a lot of time manually
navigating to wherever the files you want to visualize are. 


1. Remove the space from the name of the Application (i.e. VMD1.9.2). This is best done from the Finder (spaces in the command line are a pain to work with).

2. Add this to your `~/.bash_profile`:
```
vmdappdir='/Applications/VMD1.9.2.app/Contents'
alias vmd='"$vmdappdir/Resources/VMD.app/Contents/MacOS/VMD" $*'
```
**NOTE:** VMD only has 32-bit version for Mac OS X. This is bad because you'll only be able to open
trajectories that are half the size of your RAM. More info on this issue [here](http://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/26606.html).

**I've come up with a solution to this** :heavy_exclamation_mark:
Copy or download [this](https://github.com/jeiros/Scripts/blob/master/load_big_trajs.sh) script and save it somewhere in your 
`$PATH` (for example, `usr/local/bin`). Then make it an executable with `chmod +x load_big_trajs.sh`.
**Take a look at line 34 in the script and change the path to your VMD executable if it's different!**
If you've followed the previous steps it should be the same as the one in the script, though. Now you'll be able to load multiple trajectory files 
from the command line like so: `load_big_trajs.sh topology.prmtop trajectories*.nc` (provided you have a sensible naming scheme of your trajectories
and they show up sequentially. Check this by doing `ls trajectories*.nc`).

<p align="center">
    <img src="https://github.com/jeiros/Scripts/blob/master/misc/thumbsup.gif"/>
</p>


##Installing UCSF Chimera
[Here](https://www.cgl.ucsf.edu/chimera/download.html) are all the downloadables.

Same procedure as VMD. Click on .dmg file, drag Chimera icon to the Applications folder.

The binary is in `/Applications/Chimera.app/Contents/MacOS/chimera`.

### How to use chimera from the command line
Either create an alias in your `.bash_profile`:
```
alias chimera = '/Applications/Chimera.app/Contents/MacOS/chimera'
```
Or do a symbolic link to somewhere in your `$PATH`:
```
sudo ln -s /Applications/Chimera.app/Contents/MacOS/chimera /usr/local/bin/chimera
```
Either option should work. I *think* `/usr/local/bin` is in the default $PATH, check it by typing `echo $PATH` and looking for it. If not create it (with sudo)
and append it to the $PATH in your `~/.bash_profile`.
    


# Extra things

## Generic Note 1

The new OS X El Capitan includes a new "rootless" mode that makes certain system directories
read-only even for admins. This is a bit annoying, because one of these is `/usr/` (and thus everything inside it can't be touched). 
The exception to this is `/usr/local`).
This can be disabled following [these](http://apple.stackexchange.com/questions/196224/unix-ln-s-command-not-permitted-in-osx-el-capitan-beta3) instructions but I haven't tried it myself.


## Generic Note 2

Another "feature" that has been added is the creation of bash sessions. This creates a `~/.bash_sessions/` directory
in your home that keeps storing files. It can be annoying since these are read
every time you launch a terminal and in the end it can make it slower to start. This can be disabled if you do:
```
touch ~/.bash_sessions_disable
rm -rf ~/.bash_sessions/
```

## Generic Note 3

Not OS X specific but I find it useful to be able to ssh withouth entering passwords all the time. 
This is handy when you are working with several machines (yours, the HPC, your laptop, other local machines with GPUs...)
[This](http://www.linuxproblem.org/art_9.html) link explains how to do it. 
If properly set up, you'll be able to ssh and scp across diferent machines without having to enter the password every time.
Super useful when you add this to scripts.



Written by: [Juan Eiros](mailto:je714@ic.ac.uk)

I will keep adding stuff that I find useful.

