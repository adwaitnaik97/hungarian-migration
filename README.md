# Hungarian-migration
Refactor and migration of the [Hungarian algorithm](https://github.com/mcximing/hungarian-algorithm-cpp) to C++17. This is lincensed under the [BSD-2-Clause](https://github.com/adwaitnaik97/hungarian-migration/blob/main/LICENSE) License.

# Installation

1. Install CMake 3.25.2 or higher to compile.

### Step 1: Remove existing CMake(if necessary)

```bash
sudo apt remove --purge cmake
```

### Step 2: Download CMake 3.25.2

```bash
cd /tmp
wget https://github.com/Kitware/CMake/releases/download/v3.25.2/cmake-3.25.2.tar.gz
```

### Step 3: Extract the Archive

```bash
tar -xvzf cmake-3.25.2.tar.gz`
cd cmake-3.25.2
```

### Step 4: Build and install CMake

```bash
./bootstrap
make
sudo make install
```
### Step 5: Clean-up the files in /tmp

```bash 
cd /tmp
rm -rf cmake-3.25.2 cmake-3.25.2.tar.gz
```

### Step 6: Verify the installation

```bash
cmake --version #cmake version 3.25.2 
```

2. Install `gcc/g++` version `11`

### Step 1: Add the Toolchain PPA

```bash
sudo apt update
sudo apt install -y software-properties-common
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt update
```

### Step 2: Install GCC and G++

It is recommended to install gcc/g++ version 11

```bash
sudo apt install -y gcc-11 g++-11
```

### Step 3: Verify the installation

```bash
gcc-11 --version #gcc (Ubuntu 11.4.0-2ubuntu1~18.04) 11.4.0
g++-11 --version #g++ (Ubuntu 11.4.0-2ubuntu1~18.04) 11.4.0
```
### Step 4: Set as the default compiler

```bash
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 100
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-11 100
```

3. Install Google Test (gtest) from source.

### Step 1: Install prerequisites

```bash
sudo apt update
sudo apt install -y build-essential cmake git
```

### Step 2: Clone the Google Test repository

```bash
cd /tmp
git clone https://github.com/google/googletest.git
```

### Step 3: Create a build directory

```bash
cd googletest
mkdir build
cd build
```

### Step 4: Build Google Test

```bash
cmake ..
make
```

### Step 5: Install Google Test

```bash
sudo make install
```

### Step 6: Clean-up the files in /tmp

```bash
cd /tmp
rm -rf googletest
```

### Step 7: Verify the installation

```bash
ls /usr/local/lib | grep gtest # Should list gtest libraries
```

# Compilation & Usage

### Step 1: Compile using Cmake

```bash
mkdir build
cd build
cmake .. -DBUILD_TESTS=ON
make
```

### Step 2: Testing

```bash
ctest --verbose
```



