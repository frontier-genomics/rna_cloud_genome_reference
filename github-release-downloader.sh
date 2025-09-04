#!/bin/bash

# usage
# ./github-release-downloader.sh <GITHUB_ORG> <GITHUB_REPO> <RELEASE> <SEARCH_PATTERN> <OUTPUT_PATH>
# ./github-release-downloader.sh gohugoio hugo "1.2.3" "_linux-amd64.tar.gz" ./hugo.tar.gz

declare -r GITHUB_AUTH_HEADER="Authorization: Bearer ${GH_PAT}"
declare -r GITHUB_ORG=${GITHUB_ORG:-$1}
declare -r GITHUB_REPO=${GITHUB_REPO:-$2}
declare -r RELEASE=${RELEASE:-$3}
declare -r SEARCH_PATTERN=${SEARCH_PATTERN:-$4}
declare -r OUTPUT_PATH=${OUTPUT_PATH:-$5}

# Variables
declare RELEASES
declare ASSETS
declare ASSET

echo "Downloading..."

echo "   GITHUB_ORG       : ${GITHUB_ORG}"
echo "   GITHUB_REPO      : ${GITHUB_REPO}"
echo "   RELEASE          : ${RELEASE}"
echo "   SEARCH_PATTERN   : ${SEARCH_PATTERN}"
echo "   OUTPUT_PATH      : ${OUTPUT_PATH}"

# Retrieve
RELEASES=$(curl -H "$GITHUB_AUTH_HEADER" -sL https://api.github.com/repos/${GITHUB_ORG}/${GITHUB_REPO}/releases)
ASSETS=$(echo ${RELEASES} | jq -r ".[] | select(.name == \"$RELEASE\") | .assets")
ASSET=$(echo $ASSETS | jq -r ".[] | select(.name | contains(\"$SEARCH_PATTERN\"))")
if [ -z "$ASSET" ]; then
  echo "❌ No asset found for pattern: $SEARCH_PATTERN"
  exit 1
fi
ASSET_NAME=$(echo ${ASSET} | jq -r ".name" | head -n 1)
ASSET_URL=$(echo ${ASSET} | jq -r ".url" | head -n 1)
OUTPUT_FILE_PATH="${OUTPUT_PATH}/${ASSET_NAME}"

echo "   NAME               : ${ASSET_NAME}"
echo "   URL                : ${ASSET_URL}"
echo "   TO                 : ${OUTPUT_FILE_PATH}"
curl -L -H "$GITHUB_AUTH_HEADER" -H "Accept: application/octet-stream" -o "${OUTPUT_FILE_PATH}" "${ASSET_URL}"
echo "Downloaded: ${OUTPUT_FILE_PATH}"
