#/bin/bash

if [ -e "$HOME/.ppi_feature_trees.cfg" ]
then
  CONFIG_FILE="$HOME/.ppi_feature_trees.cfg"
else
  GITHOME=$(git rev-parse --show-toplevel)
  CONFIG_FILE=$GITHOME/Configs/$USER.cfg
fi
