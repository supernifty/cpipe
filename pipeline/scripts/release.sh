#!/bin/bash

function err() {
    echo 
    echo "ERROR: $1"
    echo
    exit 1
}

echo
echo "================================================="
echo
echo " Melbourne Genomics Release Script"
echo
echo "================================================="
echo "
What this script does:

  1  allows you to enter a version number
  2  allows you to enter release notes for the release, based
     on changelog from git
  3  creates a tag for the version in the git repository
  4  pushes tag to main repository
"

read -p "Please enter the version number for the release: "

VERSION=$REPLY
echo
echo "Version is $VERSION"
echo
echo "
    Press enter now to open release notes editor. The editor will open with the 
    git changelog and past release notes in separate panes of the editor.
"

read

PREV_VERSION=`cat version.txt`

echo "==================================================

Melbourne Genomics Demonstration Project Pipeline

Version $VERSION Release Notes

==================================================

" > notes.tmp.txt

git log ${PREV_VERSION}..HEAD | grep -v '^Author' | grep -v '^commit' | grep -v '^Date' >> notes.tmp.txt

#vim -O3 notes.tmp.txt docs/ReleaseNotes-${PREV_VERSION}.txt git.tmp.txt

touch .notes.tmp.txt.timestamp

vi notes.tmp.txt

if [ .notes.tmp.txt.timestamp -nt notes.tmp.txt ];
then
        echo "Release notes not saved - aborting ...."
        echo
        exit 1
fi

cp notes.tmp.txt docs/ReleaseNotes-$VERSION.txt || err "Unable to copy release notes to docs directory"

git add docs/ReleaseNotes-$VERSION.txt || err "Unable to add new release notes to git repo"
echo $VERSION > version.txt || err "Unable to update version file"
git commit -m "Release notes for version $VERSION" version.txt docs/ReleaseNotes-$VERSION.txt || err "Unable to commit to git repository"
git tag $VERSION || err "Unable to create tag for version $VERSION"
git push || err "Unable to push code to origin repo"
git push origin $VERSION || err "Unable to push tag $VERSION to origin repo"

rm .notes.tmp.txt.timestamp notes.tmp.txt 

echo
echo "Done."
echo

