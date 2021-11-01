#!/bin/bash
find . -type d -name 'Undete*' -exec sh -c '
for file do
mv "$file" "${file// /_}"
done
' exec-sh {} +
