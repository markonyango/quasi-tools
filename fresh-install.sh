#!/bin/bash
sudo apt-get install build-essential
sudo apt-get install git subversion
sudo apt-get install m4

# We also need some security measures
sudo apt-get install clamav ufw denyhosts psad rkhunter chkrootkit logwatch

sudo ufw allow ssh
sudo ufw enable

sudo rkhunter --update
sudo rkhunter --propupd

# Install the GMP/MPFR/MPC library
wget https://gmplib.org/download/gmp/gmp-6.1.2.tar.bz2
tar -xjf gmp-6.1.2.tar.bz2
cd gmp-6.1.2
./configure --prefix=$HOME
make
make check
make install
cd ~

wget http://www.mpfr.org/mpfr-current/mpfr-3.1.5.tar.bz2
tar -xjf mpfr-3.1.5.tar.bz2
cd mpfr-3.1.5
./configure --prefix=$HOME --with-gmp-include=$HOME/include --with-gmp-lib=$HOME/lib --enable-thread-safe
make
make check
make install
cd ~

wget ftp://ftp.gnu.org/gnu/mpc/mpc-1.0.3.tar.gz
tar -xzf mpc-1.0.3.tar.gz
cd mpc-1.0.3
./configure --prefix=$HOME --with-gmp-include=$HOME/include --with-gmp-lib=$HOME/lib --with-mpfr-include=$HOME --with-mpfr-lib=$HOME/lib --enable-valgrind-tests
make
make check
make install
cd ~

# Import new libraries
sudo touch /etc/ld.so.conf.d/local.conf
sudo echo '$HOME/lib' > /etc/ld.so.conf.d/local.conf
sudo ldconfig


# Install the R Environment
sudo apt-get install r-base
sudo apt-get install libxml2 libxml2-dev

mkdir ~/R-library
sudo sed -i 's/^R_LIBS_USER=~\/R-library/' /etc/R/Renviron
sudo sed -i 's/\/usr\/local\/lib\/R\/site-library/~\/R-library/' /etc/R/Renviron
echo 'R_LIBS_USER="~/R-library"' > ~/.Renviron

R --vanilla --silent --interactive <<RSCRIPT

source("http://bioconductor.org/biocLite.R")
biocLite(ask=FALSE, suppressUpdates=T)
biocLite(pkgs=c("DESeq","baySeq","edgeR","hexbin"),ask=FALSE, suppressUpdates=T)
quit()
RSCRIPT

# Cloning GIT repository
git clone https://github.com/markonyango/quasi-tools.git
git clone https://github.com/markonyango/quasi-express.git

# Install QUASI-tools
sh ./quasi-tools/install.sh

# Install QUASI-express
cd quasi-express
npm install

