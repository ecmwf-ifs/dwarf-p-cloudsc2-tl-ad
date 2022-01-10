#!/bin/bash

MODULES=( )
# MODULES=( daint-gpu cray-python/3.6.5.7 cudatoolkit )
PYTHON=python3.9
PIP_UPGRADE=0
VENV=venv
FRESH_INSTALL=0

function install()
{
  # activate environment
  source $VENV/bin/activate

  # upgrade pip and setuptools
  if [ "$PIP_UPGRADE" -gt 0 ]; then
    pip install --upgrade pip setuptools
  fi

  # install cloudsc2py
  pip install -e .

  # install gt sources
  python -m gt4py.gt_src_manager install

  # install development packages
  pip install -r requirements_dev.txt

  # deactivate environment
  deactivate

  # On OSX only: change matplotlib backend from macosx to TkAgg
  if [[ "$OSTYPE" == "darwin"* ]]; then
    cat $VENV/lib/$PYTHON/site-packages/matplotlib/mpl-data/matplotlibrc | \
      sed -e 's/^backend.*: macosx/backend : TkAgg/g' > /tmp/.matplotlibrc && \
      cp /tmp/.matplotlibrc $VENV/lib/$PYTHON/site-packages/matplotlib/mpl-data/matplotlibrc && \
      rm /tmp/.matplotlibrc
  fi
}

for MODULE in "${MODULES[@]}"
do
  module load $MODULE
done

if [ "$FRESH_INSTALL" -gt 0 ]
then
  echo -e "Creating new environment..."
  rm -rf $VENV
  $PYTHON -m venv $VENV
fi

install || deactivate

echo -e ""
echo -e "Command to activate environment:"
echo -e "\t\$ source $VENV/bin/activate"
echo -e ""
echo -e "Command to deactivate environment:"
echo -e "\t\$ deactivate"
echo -e ""
