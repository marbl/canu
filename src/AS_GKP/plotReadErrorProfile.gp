#! /usr/bin/env gnuplot

set multiplot
	
	set title "Read Error Profile"
	set lmargin 5
	set rmargin 5
	set tmargin 5
	set bmargin 5
	
	set key left top
	
	set xtics nomirror
	set ytics nomirror
	set xlabel "Base"
	set ylabel "Relative Errors (%)"
	set y2label "Average Read Position (%)"
	
	plot `echo \"$READ_ERROR_PROFILE\"` \
		using 1:($2 * 100) \
		title "Base vs. Relative Errors" \
		axis x1y1 \
		with lines
	
	set key right top
	
	set noxtic
	set noytic
	set y2tics nomirror
	
	#plot `echo \"$READ_ERROR_PROFILE\"` \
		using 1:($3 * 100) \
		title "Base vs. Average Read Position" \
		smooth bezier \
		axis x1y2 \
		with lines \
		linestyle 3
	
	plot `echo \"$READ_ERROR_PROFILE\"` \
		using 1:($4 * 100) \
		title "Base vs. % Mismatch Error" \
		smooth bezier \
		axis x1y2 \
		with lines \
		linestyle 3
	
	#plot `echo \"$READ_ERROR_PROFILE\"` \
		using 1:($5 * 100) \
		title "Base vs. % Insertion Error" \
		smooth bezier \
		axis x1y2 \
		with lines \
		linestyle 5
	
	#plot `echo \"$READ_ERROR_PROFILE\"` \
		using 1:($6 * 100) \
		title "Base vs. % Deletion Error" \
		smooth bezier \
		axis x1y2 \
		with lines \
		linestyle 6
	
unset multiplot
