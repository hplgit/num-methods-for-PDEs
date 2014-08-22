#!/bin/sh
# Combine 2 or 4 pdf images to one
if [ $# -eq 2 ]; then
  nup=2x1
else
  nup=2x2
fi
pdftk $@ output tmp.pdf
pdfnup --nup $nup tmp.pdf
mv -f tmp-nup.pdf tmp.pdf
pdfcrop tmp.pdf
mv -f tmp-crop.pdf tmp.pdf
/bin/ls tmp.pdf
