#!/bin/sh


if [ -f addnoise ]
then
  rm addnoise
fi

goimports -w .
go fmt
go build

if [ -f addnoise ]
then
  ./addnoise

#  feh -d -g+1000+0 s.png
  feh -d -g+1000+0 c.png
fi
