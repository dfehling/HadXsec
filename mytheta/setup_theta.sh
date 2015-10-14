#!/bin/bash
# A script to checkout the latest version of theta into a subfolder named theta and compile it

#Make a directory to store everything in. This could be made configurable I guess
echo "Creating theta directory"
mkdir theta

#Check out the latest recommended version of theta
echo "Checking out theta code"
svn co https://ekptrac.physik.uni-karlsruhe.de/public/theta/tags/testing theta --username=guest

#Change to the theta directory and compile with normal options
echo "Building theta code, this may take a few minutes"
cd theta
make -j 4

echo "Job complete"
