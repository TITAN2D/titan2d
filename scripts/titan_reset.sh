#! /bin/bash

for f in *; do
   case $f in
       *data)
           ;;
       funky*)
           ;;
       *lhs*)
           ;;
       *.pl)
           ;;
       titan_reset*)
           ;;
       *)
          if [ -f $f ]; then
             rm -f $f
          fi
          ;;
   esac
done
cp ../titan .
