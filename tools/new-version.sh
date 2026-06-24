#!/usr/bin/env sh
if [ $# -ne 1 ]; then
    echo "usage: $0 version"
    exit 1
fi
NEWV=$1
D=$(date +'%Y-%m-%d')
sed -i -e "s/VERSION \([[:digit:]]\+[\.]\?\)\{3\} /VERSION ${NEWV} /g" CMakeLists.txt
sed -i -e "s/version: \([[:digit:]]\+[\.]\?\)\{3\}/version: ${NEWV}/g" CODE
sed -i -e "s/release-date: .*$/release-date: ${D}/g" CODE
sed -i -e "s/Version: \([[:digit:]]\+[\.]\?\)\{3\}.*$/Version: ${NEWV} (${D})/g" README.md
sed -i -e "s/version = \"\([[:digit:]]\+[\.]\?\)\{3\}\"/version = \"${NEWV}\"/g" fpm.toml
