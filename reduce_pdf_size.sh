#!/bin/bash
set -eu
# Use as ./reduce_pdf_size.sh in.pdf out.pdf

if [[ $# -lt 2 ]]; then
    printf "Provide two arguments\n"
    printf "Use as ./reduce_pdf_size.sh in.pdf out.pdf\n"
    exit 1
fi

in=$1
out=$2

ghostscript -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 \
	-dPDFSETTINGS=/ebook -dNOPAUSE -dQUIET -dBATCH \
	-sOutputFile=${out} ${in}

printf "Done\n"