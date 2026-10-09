#!/usr/bin/env sh

#!/usr/bin/env sh
if [ $# -ne 1 ]; then
    echo "usage $0 VERSION"
    exit 1
fi
VERSION=$1


git tag -m "Version ${VERSION}" v${VERSION}
git push origin v${VERSION}
git archive --format tar.gz --prefix qrupdate-ng-${VERSION}/ -o ../qrupdate-ng-${VERSION}.tar.gz v${VERSION}

