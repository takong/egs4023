#!/bin/bash

unzip egs4023-master.zip

cd egs4023-master

tar -xjvf resource.tar.bz2

bash link_som_exec.sh

bash runExtreme.sh

ls -lsh -R 
