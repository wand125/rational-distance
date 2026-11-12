#!/bin/sh
# exec keeps this harness plus its one child to at most two processes.
set -eu
case "$0" in
    /*) script_path=$0 ;;
    *) script_path=$PWD/$0 ;;
esac
exec python3 "${script_path%/*}/accept.py" "$@"
