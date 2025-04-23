#!/bin/bash

eval "$(ssh-agent -s)"

ssh-add ~/.ssh/git_key

git fetch origin
git reset --hard origin/main
