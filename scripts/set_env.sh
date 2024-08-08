#!/usr/bin/sh

# set environment variables
export ROSE_DEPENDENCIES_ROOT=$(find /glade/work/al*  -maxdepth 1 -type d -name  "ROSE")
export ROSE_INSTALL_PATH=$(find /glade/work/jd*  -maxdepth 1 -type d -name  "ROSE")/install
export PROSE_REPO_PATH=$(git rev-parse --show-superproject-working-tree --show-toplevel | head -1)
export LD_LIBRARY_PATH="${ROSE_INSTALL_PATH}"/lib:"${ROSE_DEPENDENCIES_ROOT}"/jdk1.8.0_241/jre/lib/amd64/server:"${ROSE_DEPENDENCIES_ROOT}"/flex-2.6.4/install/lib:"${ROSE_DEPENDENCIES_ROOT}"/boost/1_67_0/install/lib:"${ROSE_DEPENDENCIES_ROOT}"/gcc/7.4.0/install/lib64:"${LD_LIBRARY_PATH}"
export PATH="${PROSE_REPO_PATH}"/scripts:"${ROSE_INSTALL_PATH}"/bin:"${ROSE_DEPENDENCIES_ROOT}"/jdk1.8.0_241/bin/:"${ROSE_DEPENDENCIES_ROOT}"/jdk1.8.0_241/jre/bin/:"${ROSE_DEPENDENCIES_ROOT}"/automake-1.16.2/install/bin/:"${ROSE_DEPENDENCIES_ROOT}"/flex-2.6.4/install/bin:"${ROSE_DEPENDENCIES_ROOT}"/gettext-0.19.7/install/bin:"${ROSE_DEPENDENCIES_ROOT}"/gcc/7.4.0/install/bin:"${PATH}"
export JRE_HOME="${ROSE_DEPENDENCIES_ROOT}"/jdk1.8.0_241/jre
export JAVA_BINDIR="${ROSE_DEPENDENCIES_ROOT}"/jdk1.8.0_241/bin
export JAVA_HOME="${ROSE_DEPENDENCIES_ROOT}"/jdk1.8.0_241
export SDK_HOME="${ROSE_DEPENDENCIES_ROOT}"/jdk1.8.0_241
export JDK_HOME="${ROSE_DEPENDENCIES_ROOT}"/jdk1.8.0_241
export JAVA_ROOT="${ROSE_DEPENDENCIES_ROOT}"/jdk1.8.0_241
export ROSE_EXE_PATH="${ROSE_INSTALL_PATH}"/bin
export PROSE_PLUGIN_PATH="${PROSE_REPO_PATH}"/src/cpp

export TMPDIR=$(find /glade/d*/scratch/ -maxdepth 1 -name ${USER})/temp
mkdir -p $TMPDIR

module purge
module reset
module load conda

conda activate $(find /glade/work/jd*/conda-envs  -maxdepth 1 -type d -name  "my-env")