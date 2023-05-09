## JAGS and JAGS-WIENER installation steps for OS

#### 1. Install JAGS:

```shell
sudo apt-get update
sudo apt-get install jags
```

#### 2. download JAGS-WIENER

```shell
cd ~
mkdir install
cd install
http://downloads.sourceforge.net/project/jags-wiener/JAGS-WIENER-MODULE-1.1.tar.gz
tar -zxvf JAGS-WIENER-MODULE-1.1.tar.gz
cd JAGS-WIENER-MODULE-1.1
```

#### 3. Install JAGS-WIENER

Configure and compile the source code for your system with the following commands in a terminal window:

```shell
./configure && make
```
When this is done,install the libraries on your system with the following command, which usually requires root privileges: 

```shell
sudo make install
```

#### 4. Test installation (run commands within the R terminal):

```R
library(R2jags)

load.module("wiener")
load.module("dic")
list.modules()
```

A successful installation of JAGS should load the modules: "basemod" "bugs" "dic". A successful installation of JAGS-WIENER that is found by JAGS will load "wiener"
