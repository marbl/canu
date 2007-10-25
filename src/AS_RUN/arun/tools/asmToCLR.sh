#!/bin/sh
 grep -A 5 "{AFG" $1 | grep -E "acc|clr" | cut -d':' -f2 | tr -d '(' | tr -d ')' | paste -s -d" \n" | tr ' ' ',' | awk -F, '{print $1,$2,$3}'
