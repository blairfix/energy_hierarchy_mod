#!/bin/bash

# move executables from source directory to executables directory

find ./source_code/ -maxdepth 1 -iname '*.o' -exec mv -t ./executables/ {} +

find ./source_code/ -maxdepth 1 -type f ! -name "*.*" -exec mv -t ./executables/ {} \+
