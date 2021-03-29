#!/bin/bash

find . -type f -regex '.*\.\(hh\|cc\|cpp\|h\|C\)' -exec clang-format -i "{}" \;

