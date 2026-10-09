#!/usr/bin/env sh
if [ $# -ne 1 ]; then
    echo "usage $0 VERSION"
    exit 1
fi
VERSION=$1
MAJOR=$(echo $VERSION | cut -d. -f1)
MINOR=$(echo $VERSION | cut -d. -f2)
PATCH=$(echo $VERSION | cut -d. -f3)

DATE=$(date +"%Y-%m-%d")
sed -i -e "s/\(PROJECT(qrupdate-ng VERSION\) [0-9\.]* \(LANGUAGES Fortran)\)/\1 ${VERSION} \2/g" CMakeLists.txt
sed -i -e "s/version: [0-9\.]*/version: ${VERSION}/g" CODE
sed -i -e "s/release-date: [0-9-]*/release-date: ${DATE}/g" CODE
sed -i -e "s/QRUPDATE_VERSION \"[0-9\.]*\"/QRUPDATE_VERSION \"${VERSION}\"/g" src/fpm_include/qrupdate_config.h
sed -i -e "s/QRUPDATE_VERSION_MAJOR [0-9]*/QRUPDATE_VERSION_MAJOR ${MAJOR}/g" src/fpm_include/qrupdate_config.h
sed -i -e "s/QRUPDATE_VERSION_MINOR [0-9]*/QRUPDATE_VERSION_MINOR ${MINOR}/g" src/fpm_include/qrupdate_config.h
sed -i -e "s/QRUPDATE_VERSION_PATCH [0-9]*/QRUPDATE_VERSION_PATCH ${PATCH}/g" src/fpm_include/qrupdate_config.h
sed -i -e "s/Version: [0-9\.]* ([0-9-]*)/Version: ${VERSION} (${DATE})/g" README.md
sed -i -e "s/version = \"[0-9\.]*\"/version = \"${VERSION}\"/g" fpm.toml
