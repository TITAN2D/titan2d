#!/bin/bash

aclocal
autoheader
autoconf
automake --add-missing --copy
