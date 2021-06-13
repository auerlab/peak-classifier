#!/bin/sh -e

awk '$1 ~ /^[0-9]+$/' filtered.bed > autosomes.bed
