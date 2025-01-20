#!/bin/bash
# run this script to set up a xtb environment
# requirements: $QCxMSHOME is set to `pwd`
if [ -z "${QCxMS2HOME}" ]; then
   QCxMS2HOME="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../../"
fi

# set up path for QCxMS2, using the xtb directory and the users home directory
QCxMS2HOME=${QCxMS2HOME}/@datadir@:${HOME}

# to include the documentation we include our man pages in the users manpath
MANPATH=${MANPATH}:${QCxMS2HOME}/@mandir@

# finally we have to make the binaries
PATH=${PATH}:${QCxMS2HOME}/@bindir@

# enable package config for QCxMS 
PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:${QCxMS2HOME}/@libdir@/pkgconfig

export PATH QCxMS2PATH MANPATH PKG_CONFIG_PATH
