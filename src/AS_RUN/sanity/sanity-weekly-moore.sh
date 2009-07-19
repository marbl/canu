#!/bin/sh

date=$1

email="bwalenz@jcvi.org jmiller@jcvi.org skoren@jcvi.org gsutton@jcvi.org jjohnson@jcvi.org"
email="bwalenz@jcvi.org jmiller@jcvi.org skoren@jcvi.org"
email="bwalenz@jcvi.org"

#  Weekly Moore
#
perl sanity.pl assemble $date moore \
  SPECFILES/moore/M001.spec \
  SPECFILES/moore/M002.spec \
  SPECFILES/moore/M003.spec \
  SPECFILES/moore/M005.spec \
  SPECFILES/moore/M006.spec \
  SPECFILES/moore/M007.spec \
  SPECFILES/moore/M008.spec \
  SPECFILES/moore/M010.spec \
  SPECFILES/moore/M011.spec \
  SPECFILES/moore/M012.spec \
  SPECFILES/moore/M013.spec \
  SPECFILES/moore/M014.spec \
  SPECFILES/moore/M015.spec \
  SPECFILES/moore/M016.spec \
  SPECFILES/moore/M017.spec \
  SPECFILES/moore/M018.spec \
  SPECFILES/moore/M019.spec \
  SPECFILES/moore/M020.spec \
  SPECFILES/moore/M022.spec \
  SPECFILES/moore/M023.spec \
  SPECFILES/moore/M024.spec \
  SPECFILES/moore/M025.spec \
  SPECFILES/moore/M026.spec \
  SPECFILES/moore/M027.spec \
  SPECFILES/moore/M028.spec \
  SPECFILES/moore/M029.spec \
  SPECFILES/moore/M030.spec \
  SPECFILES/moore/M032.spec \
  SPECFILES/moore/M033.spec \
  SPECFILES/moore/M035.spec \
  SPECFILES/moore/M036.spec \
  SPECFILES/moore/M037.spec \
  SPECFILES/moore/M038.spec \
  SPECFILES/moore/M039.spec \
  SPECFILES/moore/M040.spec \
  SPECFILES/moore/M043.spec \
  SPECFILES/moore/M046.spec \
  SPECFILES/moore/M047.spec \
  SPECFILES/moore/M048.spec \
  SPECFILES/moore/M049.spec \
  SPECFILES/moore/M051.spec \
  SPECFILES/moore/M053.spec \
  SPECFILES/moore/M054.spec \
  SPECFILES/moore/M055.spec \
  SPECFILES/moore/M056.spec \
  SPECFILES/moore/M057.spec \
  SPECFILES/moore/M058.spec \
  SPECFILES/moore/M059.spec \
  SPECFILES/moore/M061.spec \
  SPECFILES/moore/M062.spec \
  SPECFILES/moore/M063.spec \
  SPECFILES/moore/M065.spec \
  SPECFILES/moore/M066.spec \
  SPECFILES/moore/M068.spec \
  SPECFILES/moore/M069.spec \
  SPECFILES/moore/M070.spec \
  SPECFILES/moore/M071.spec \
  SPECFILES/moore/M072.spec \
  SPECFILES/moore/M073.spec \
  SPECFILES/moore/M074.spec \
  SPECFILES/moore/M075.spec \
  SPECFILES/moore/M076.spec \
  SPECFILES/moore/M078.spec \
  SPECFILES/moore/M081.spec \
  SPECFILES/moore/M082.spec \
  SPECFILES/moore/M083.spec \
  SPECFILES/moore/M084.spec \
  SPECFILES/moore/M085.spec \
  SPECFILES/moore/M086.spec \
  SPECFILES/moore/M087.spec \
  SPECFILES/moore/M088.spec \
  SPECFILES/moore/M092.spec \
  SPECFILES/moore/M093.spec \
  SPECFILES/moore/M094.spec \
  SPECFILES/moore/M095.spec \
  SPECFILES/moore/M096.spec \
  SPECFILES/moore/M098.spec \
  SPECFILES/moore/M101.spec \
  SPECFILES/moore/M102.spec \
  SPECFILES/moore/M103.spec \
  SPECFILES/moore/M104.spec \
  SPECFILES/moore/M105.spec \
  SPECFILES/moore/M106.spec \
  SPECFILES/moore/M108.spec \
  SPECFILES/moore/M109.spec \
  SPECFILES/moore/M110.spec \
  SPECFILES/moore/M111.spec \
  SPECFILES/moore/M112.spec \
  SPECFILES/moore/M113.spec \
  SPECFILES/moore/M114.spec \
  SPECFILES/moore/M115.spec \
  SPECFILES/moore/M116.spec \
  SPECFILES/moore/M117.spec \
  SPECFILES/moore/M118.spec \
  SPECFILES/moore/M119.spec \
  SPECFILES/moore/M120.spec \
  SPECFILES/moore/M121.spec \
  SPECFILES/moore/M122.spec \
  SPECFILES/moore/M123.spec \
  SPECFILES/moore/M126.spec \
  SPECFILES/moore/M127.spec \
  SPECFILES/moore/M128.spec \
  SPECFILES/moore/M129.spec \
  SPECFILES/moore/M131.spec \
  SPECFILES/moore/M133.spec \
  SPECFILES/moore/M134.spec \
  SPECFILES/moore/M135.spec \
  SPECFILES/moore/M136.spec \
  SPECFILES/moore/M137.spec \
  SPECFILES/moore/M138.spec \
  SPECFILES/moore/M139.spec \
  SPECFILES/moore/M140.spec \
  SPECFILES/moore/M141.spec \
  SPECFILES/moore/M142.spec \
  SPECFILES/moore/M143.spec \
  SPECFILES/moore/M144.spec \
  SPECFILES/moore/M145.spec \
  SPECFILES/moore/M152.spec \
  SPECFILES/moore/M157.spec \
  SPECFILES/moore/M159.spec \
  SPECFILES/moore/M160.spec \
  SPECFILES/moore/M169.spec \
  SPECFILES/moore/M170.spec \
  SPECFILES/moore/M172.spec \
  SPECFILES/moore/M173.spec \
  SPECFILES/moore/M175.spec \
  SPECFILES/moore/M176.spec \
  SPECFILES/moore/M180.spec \
  $email

#  Resubmit ourself for tomorrow.

exit
