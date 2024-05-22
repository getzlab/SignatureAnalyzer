#!/bin/bash

# method name
METHOD_NAME=signature_analyzer

# docker build parameters (most likely, the defaults are OK)
EXTRA_DOCKER_BUILD_ARGS=""
DOCKER_BUILD_PATH="."

#
# DO NOT EDIT PAST THIS POINT UNLESS YOU HAVE VERY GOOD REASON TO DO SO
#

# automatically create method version number
NCOMMIT=$(git rev-list --count HEAD)
BRANCH=$(git rev-parse --abbrev-ref HEAD)

if [[ $BRANCH == "master" ]]; then
	BRANCH=""
else
	BRANCH+="_"
fi
VERSION=${BRANCH}v${NCOMMIT}

# build method docker
docker build -t getzlab/$METHOD_NAME:${VERSION} ${EXTRA_DOCKER_BUILD_ARGS} ${DOCKER_BUILD_PATH}

# push method docker
docker tag getzlab/$METHOD_NAME:${VERSION} gcr.io/broad-getzlab-workflows/$METHOD_NAME:${VERSION}
docker push gcr.io/broad-getzlab-workflows/$METHOD_NAME:${VERSION}
