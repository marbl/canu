#if [ -f toProcess.txt ]; then
#
###########################################################################
#
# This file is part of Celera Assembler, a software program that 
# assembles whole-genome shotgun reads into contigs and scaffolds.
# Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received (LICENSE.txt) a copy of the GNU General Public 
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
###########################################################################
#
# $Id: filterExplain.sh,v 1.4 2005-03-22 19:48:58 jason_miller Exp $
#
#  rm toProcess.txt
#fi

Offset=${1}

if [ -z ${Offset} ] || [ ${Offset} == "bell-style" ]; then
  echo "Need offset of first result value in lines"
  exit
fi


if [ ! -f toProcess.txt ]; then 
  for file in `ls 0[0-9][0-9]_all.csv`; do
    chr=${file%%_*}
    tail +9 ${file} | sed 's/,//g' | gawk -v o=${Offset} -v c=${chr} \
      '{ \
        printf("%s ", c); \
        for(i=1;i<=NF;i++) \
          printf("%s ", $i); \
        printf("\n"); \
      }' >> toProcess.txt
  done
fi

gawk -v o=${Offset} \
  'BEGIN{ \
    twok=0;tenk=0;fiftyk=0;bacEnd=0;twokOrTenk=0;ttbe=0; \
    twokOrTenkStretchedCompressed=0; twokOrTenkStretchedCompressedSat=0;\
    twokSat=0;tenkSat=0;fiftykSat=0;bacEndSat=0;twokOrTenkSat=0;ttbeSat=0; \
    stretched=0;compressed=0;inversion=0;transposition=0; \
    stretchedOnly=0;compressedOnly=0;inversionOnly=0;transpositionOnly=0; \
    mixedBad=0;mixedBadOnly=0; \
    totalBadless=0;totalGoodless=0;nada=0; \
    tt5be=0;tt5scbe=0;tt5beSat=0;tt5scbeSat=0; \
    } \
    { \
      if($(o)==0&&$(o+4)==0&&$(o+12)==0&&$(o+16)==0) \
      { \
        totalGoodless++; \
        if($(o+20)>0) \
        { \
          if($(o+28)==0&&$(o+36)==0&&$(o+44)==0) \
          { \
            compressedOnly++; \
          } \
          else \
          { \
            mixedBadOnly++; \
          } \
        } \
        else \
        { \
          if($(o+28)>0) \
          { \
            if($(o+36)==0&&$(o+44)==0) \
            { \
              inversionOnly++; \
            } \
            else \
            { \
              mixedBadOnly++; \
            } \
          } \
          else \
          { \
            if($(o+36)>0) \
            { \
              if($(o+44)==0) \
              { \
                stretchedOnly++; \
              } \
              else \
              { \
                mixedBadOnly++; \
              } \
            } \
            else \
            { \
              if($(o+44)>0) \
              { \
                transpositionOnly++; \
              } \
              else \
              { \
                nada++; \
              } \
            } \
          } \
        } \
      } \
      else \
      { \
        if($(o)>0) tenk++; \
        if($(o+4)>0) twok++; \
        if($(o)>0||$(o+4)>0||$(o+8)>0||$(o+16)>0) twokOrTenkStretchedCompressed++; \
        if($(o)>0||$(o+4)>0||$(o+8)>0||$(o+12)>0||$(o+16)>0) tt5scbe++; \
        if($(o+12)>0) fiftyk++; \
        if($(o+16)>0) bacEnd++; \
        if($(o)>0||$(o+4)>0) twokOrTenk++; \
        if($(o)>0||$(o+4)>0||$(o+16)>0) ttbe++; \
        if($(o)>0||$(o+4)>0||$(o+12)>0||$(o+16)>0) tt5be++; \
        if($(o+20)>0) compressed++; \
        if($(o+28)>0) inversion++; \
        if($(o+36)>0) stretched++; \
        if($(o+44)>0) transposition++; \
        if($(o+20)==0&&$(o+28)==0&&$(o+36)==0&&$(o+44)==0) \
        { \
          totalBadless++; \
          if($(o)>0) tenkSat++; \
          if($(o+4)>0) twokSat++; \
          if($(o+12)>0) fiftykSat++; \
          if($(o+16)>0) bacEndSat++; \
          if($(o)>0||$(o+4)>0) twokOrTenkSat++; \
          if($(o)>0||$(o+4)>0||$(o+16)>0) ttbeSat++; \
          if($(o)>0||$(o+4)>0||$(o+12)>0||$(o+16)>0) tt5beSat++; \
          if($(o)>0||$(o+4)>0||$(o+8)>0||$(o+16)>0) twokOrTenkStretchedCompressedSat++; \
          if($(o)>0||$(o+4)>0||$(o+8)>0||$(o+12)>0||$(o+16)>0) tt5scbeSat++; \
        } \
      } \
    } \
    END{ \
      print NR, "total intervals."; \
      print "Satisfied clones, ignoring confirmed unsatisfied intervals:"; \
      print twok, "intervals spanned by one or more satisfied 2k clones"; \
      print tenk, "intervals spanned by one or more satisfied 10k clones"; \
      print fiftyk, "intervals spanned by one or more satisfied 50k clones"; \
      print bacEnd, "intervals spanned by one or more satisfied BACs"; \
      print "-"; \
      print twokOrTenk, "intervals spanned by one or more satisfied 2k or 10k clones"; \
      print ttbe, "intervals spanned by one or more satisfied 2k or 10k clones or a BAC"; \
      print twokOrTenkStretchedCompressed, "intervals spanned one or more satisfied 2k, 10k, BAC, or stretched/compressed 2k or 10k clones"; \
      print tt5be, "intervals spanned by one or more satisfied 2k, 10k, 50k clones or a BAC"; \
      print tt5scbe, "intervals spanned by one or more satisfied 2k, 10k, 50k, BAC, or stretched/compressed 2k or 10k clones"; \
      print ""; \
      print "Satisfied clones and not confirmed unsatisfied intervals:"; \
      print totalBadless, "intervals with no confirmed unsatisfied intervals"; \
      print twokSat, "intervals spanned by one or more satisfied 2k clones"; \
      print tenkSat, "intervals spanned by one or more satisfied 10k clones"; \
      print fiftykSat, "intervals spanned by one or more satisfied 50k clones"; \
      print bacEndSat, "intervals spanned by one or more satisfied BACs"; \
      print "-"; \
      print twokOrTenkSat, "intervals spanned by one or more satisfied 2k or 10k clones"; \
      print ttbeSat, "intervals spanned by one or more satisfied 2k or 10k clones or a BAC"; \
      print twokOrTenkStretchedCompressedSat, "intervals spanned one or more satisfied 2k, 10k, BAC, or stretched/compressed 2k or 10k clones"; \
      print tt5beSat, "intervals spanned by one or more satisfied 2k, 10k, 50k clones or a BAC"; \
      print tt5scbeSat, "intervals spanned by one or more satisfied 2k, 10k, 50k, BAC, or stretched/compressed 2k or 10k clones"; \
      print ""; \
      print "Confirmed unsatisfied intervals, ignoring satisfied clones:"; \
      print compressed, "intervals spanned by compressed"; \
      print stretched, "intervals spanned by stretched"; \
      print inversion, "intervals spanned by inversions"; \
      print transposition, "intervals spanned by transpositions"; \
      print mixedBad, "intervals spanned by mixed bad"; \
      print "";
      print "Confirmed unsatisfied intervals and not satisfied clones:"; \
      print totalGoodless, "intervals with no spanning satisfied clones"; \
      print compressedOnly, "intervals spanned only by compressed"; \
      print stretchedOnly, "intervals spanned only by stretched"; \
      print inversionOnly, "intervals spanned only by inversions"; \
      print transpositionOnly, "intervals spanned only by transpositions"; \
      print mixedBadOnly, "intervals spanned only by mixed unsatisfied"; \
      print ""; \
      print nada, "intervals spanned by nothing"; \
    }' toProcess.txt

#head -8 000_all.csv > header.txt

#if [ -f temp.txt ]; then
#  rm temp.txt
#fi

#for file in `ls *_all.csv`; do
#  tail +9 ${file} >> temp.txt
#done

#sort -k 8n temp.txt > sortedAll.csv
