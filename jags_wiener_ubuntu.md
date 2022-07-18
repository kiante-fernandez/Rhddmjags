## JAGS and JAGS-WIENER installation steps for OS

#### 1. Install JAGS:

```bash
sudo apt-get update
sudo apt-get install jags
```

#### 2. download and install JAGS-WIENER

```bash
cd ~
mkdir install
cd install
wget http://downloads.sourceforge.net/project/jags-wiener/JAGS-WIENER-MODULE-1.1.tar.gz
tar -zxvf JAGS-WIENER-MODULE-1.1.tar.gz
cd JAGS-WIENER-MODULE-1.1
./configure --prefix=/usr/lib/x86_64-linux-gnu
make
sudo make install
```

Note that the prefix in the 3rd to last line above may change based on the install location of JAGS. This is especially likely if this error occurs:

```
DWiener.h:5:10: fatal error: distribution/ScalarDist.h: No such file or directory
#include <distribution/ScalarDist.h>
```

Try this command to find the correct prefix in the terminal:

```bash
find /usr -type f -name dic.la
```

For instance, if the result is ”/usr/lib/local/JAGS/modules-4/dic.la” Then I would use the following prefix in the 3rd to last line:

```bash
./configure --prefix=/usr/lib/local
```

#### 5. Test installation (last three commands will be within the R terminal):

```bash
jags
```
```
load.module("wiener")
load.module("dic")
list.modules()
```

A successful installation of JAGS should load the modules: basemod “ok”, bugs “ok”, and dic “ok”. A successful installation of JAGS-WIENER that is found by JAGS will load wiener “ok”.

