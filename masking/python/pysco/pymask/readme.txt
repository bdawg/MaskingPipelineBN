Hi, this is the readme to pymask, Ben's closure phase fitting package! This currently runs emcee and MultiNest to do Bayesian fitting on already-reduced closure phase data. This is a work in progress! In the near future I'll be adding support for closure phase covariances, statistically-independent closure phases and contrast detection limit calculation.

These are all included in a subfolder that should handle most of the nitty-gritty so that it can be called with very simple commands, such as are used in cpanalysis.py. 

-------------------------------------------------------------------------------------------------

Install/history documentation:

I'd like to point out that this depends on the wonderful Bayesian fitting packages emcee and MultiNest (including the latter's PyMultiNest wrapper) and therefore lots of credit goes to the authors of those packages and you should cite them on any work using their code. In addition to this, the code is very much modelled after and contains a few functions verbatim from pysco, a package I worked on with Frantz Martinache and Mike Ireland, so thanks to both of them for their contributions to that. It also includes oifits.py borrowed from Paul Boley at http://www.mpia-hd.mpg.de/homes/boley/oifits/. 

Anyway, the issue in getting pymask to work will be in installing those packages!

Here's how you do it:

emcee: this one is easy. Go to http://dan.iel.fm/emcee/ and download it, or use sudo easy_install emcee or sudo pip install emcee. 

MultiNest: this is a world of pain.

First go to http://ccpforge.cse.rl.ac.uk/gf/project/multinest/, sign up and download multinest v 2.18, which is the version I've used in putting together pymask. I'm sure other versions will work, but there could be dragons.

There's a great blog about how to install MultiNest at http://www.astro.gla.ac.uk/~matthew/blog/?p=342. Personally, I opted for the gfortran install rather than the Intel compilers, which are probably faster but seemed like a lot of effort. For my gfortran install, I changed the top of the Makefile in the main directory to read as:

#FC = mpif90 -DMPI
FC = gfortran
CC = mpicc
CXX = mpiCC
#FFLAGS +=  -w -O3
FFLAGS +=  -w -O3 -ffree-line-length-none -fPIC
CFLAGS += -O3 -DMPI -fPIC

LAPACKLIB = -llapack

NESTLIBDIR = ./

export FC CC CXX FFLAGS CFLAGS LAPACKLIB

Importantly, doing this recquires gfortran to have the appropriate libs so you need to run sudo apt-get install gfortran libblas-dev liblapack-dev. 

Anyway, that should about make it work, and if not, check the blog. 

The next issue is getting PyMultiNest to work. I went to https://github.com/JohannesBuchner/PyMultiNest and downloaded it, and did python setup.py install from the command line. Follow the installation instructions at http://johannesbuchner.github.io/PyMultiNest/install.html#installing-the-python-library as follows:

1) add to .bashrc 
export MULTINEST=/my/multinest/directory

2) go to that directory and type at the command line: make libnest3.so - this creates a dynamic library, so you have to compile it on your own machine as it has relative directory structures.

3) go to the pymultinest directory and type at the command line: make -C multinest_bridge libcnest.so

4) go to .bashrc and add 
export LD_LIBRARY_PATH=$MULTINEST:/my/pymultinest/directory/multinest_bridge

5) test the libraries as 

$ python -c 'import pymultinest'
$ python -c 'import pyapemost'

in the appropriate folders.

I wish you luck!
