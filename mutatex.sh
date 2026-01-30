#!/usr/bin/env bash

mutatex run \
  --structure complex.pdb \
  --chain A B \
  --engine foldx \
  --interface-only
