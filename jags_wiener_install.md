# JAGS and JAGS-WIENER Installation Guide for macOS
This guide provides step-by-step instructions for installing JAGS and the JAGS-WIENER module on macOS systems.

## Check and Install Homebrew
First, check if Homebrew is installed:
```bash
which brew
```

If Homebrew is not installed (no output from the above command), install it:
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

After installation, if you're using an Apple Silicon Mac (M1/M2/M3), add Homebrew to your PATH:
```bash
echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> ~/.zprofile
eval "$(/opt/homebrew/bin/brew shellenv)"
```

Verify Homebrew installation:
```bash
brew --version
```

## Prerequisites
Now install the required build tools:
```bash
brew install autoconf automake libtool
```

## 1. Install JAGS

Install JAGS using Homebrew:

```bash
brew install jags
```

## 2. Install R Packages

Open R and install the required packages:

```r
install.packages(c("rjags", "R2jags"))
```

If you encounter any package installation errors, try reinstalling with:

```r
install.packages("rjags", type="source", INSTALL_opts = "--no-compress-data")
```

## 3. Install JAGS-WIENER

### Clone the Repository

```bash
cd ~
git clone https://github.com/yeagle/jags-wiener.git
cd jags-wiener
```

### Check and Create JAGS Modules Directory

```bash
# Check if modules directory exists
if [ ! -d "/opt/homebrew/lib/JAGS/modules-4" ]; then
    echo "Creating JAGS modules directory..."
    sudo mkdir -p /opt/homebrew/lib/JAGS/modules-4
else
    echo "JAGS modules directory already exists"
fi

# Verify JAGS installation location
brew list jags

# Verify JAGS include files
ls /opt/homebrew/Cellar/jags/4.3.2/include/JAGS/distribution/

# If ScalarDist.h is not found, locate it:
find $(brew --prefix) -name "ScalarDist.h"
```

### Configure and Install

```bash
# Clean any previous attempts
make clean

# Create auxiliary files
autoreconf -fvi

# Configure with correct paths (adjust version if needed based on your JAGS installation)
./configure \
  CPPFLAGS="-I/opt/homebrew/Cellar/jags/4.3.2/include/JAGS" \
  LDFLAGS="-L/opt/homebrew/Cellar/jags/4.3.2/lib" \
  --prefix=/opt/homebrew \
  --with-jags-module=/opt/homebrew/lib/JAGS/modules-4

# Build and install
make
sudo make install
```

## 4. Verify Installation

### Check Module Location and Permissions

```bash
# Check if module was installed
ls -l /opt/homebrew/lib/JAGS/modules-4/wiener.so

# If not found, search for it
find /opt/homebrew -name "wiener.so"

# Check permissions
ls -l /opt/homebrew/lib/JAGS/modules-4/
```

### Test in R

Open R and run:

```r
library(R2jags)
load.module("wiener")
load.module("dic")
list.modules()
```

A successful installation should show:
- "basemod"
- "bugs"
- "dic"
- "wiener"

## Troubleshooting

### If Module Not Found

If you get the error "File not found: /opt/homebrew/lib/JAGS/modules-4/wiener.so":

1. Check module directory structure:
```bash
ls -la /opt/homebrew/lib/JAGS/
ls -la /opt/homebrew/lib/JAGS/modules-4/
```

2. Locate the built module:
```bash
find /opt/homebrew -name "wiener.so"
```

3. If found in a different location, create a symbolic link:
```bash
sudo ln -s /path/to/found/wiener.so /opt/homebrew/lib/JAGS/modules-4/wiener.so
```

4. Check file permissions:
```bash
sudo chmod 755 /opt/homebrew/lib/JAGS/modules-4/wiener.so
```

### If Build Fails

If you encounter build errors, try completely removing and reinstalling both JAGS and JAGS-WIENER:

```bash
# Remove current installations
brew uninstall jags
cd jags-wiener
sudo make uninstall

# Clean up any remaining files
sudo rm -rf /opt/homebrew/lib/JAGS/modules-4/wiener.so

# Reinstall JAGS
brew install jags

# Then follow the JAGS-WIENER installation steps again
```

