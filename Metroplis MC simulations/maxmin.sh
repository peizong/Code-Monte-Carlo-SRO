#!/usr/bin/env bash

[ -z "$1" ] && { echo "Needs a limit as first argument." >&2; exit 1; }
read number
echo Min: $(( $number  < $1 ? $number : $1 ))
echo Max: $(( $number  > $1 ? $number : $1 ))
