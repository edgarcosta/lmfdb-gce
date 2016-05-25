# output as png image
set terminal pngcairo size 1000,800

# save file to "benchmark.png"
set output "ab_temp_storage_ModFormGL2.png"

# This sets the aspect ratio of the graph
#set size 1, 1

# The graph title
set title "Apache Benchmark"
set label "www-central2.lmfdb.xyz/L/ModularForm/GL2/Q/holomorfic/5/6/4/a/0" at graph 0.05, 0.95

# Where to place the legend/key
#set key left top

# Draw gridlines oriented on the y axis
set grid y

# Specify that the x-series data is time data
set xdata time

# Specify the *input* format of the time data
set timefmt "%s"

# Specify the *output* format for the x-axis tick labels

set format x "%H:%M"

# Label the x-axis
set xlabel "Time, h:min"

# Label the y-axis
set ylabel "Response time, ms"

# Tell gnuplot to use tabs as the delimiter instead of spaces (default)
set datafile separator '\t'

# Plot the data
plot \
"n1hcpu4_mmap__L_ModularForm_GL2_Q_holomorphic_5_6_4_a_0_5_10000.tsv" using 2:5 with points title 'MMAP', \
"n1hcpu4_wt__L_ModularForm_GL2_Q_holomorphic_5_6_4_a_0_5_10000.tsv" using 2:5 with points title 'WT', \
"n1hcpu4_wtzl__L_ModularForm_GL2_Q_holomorphic_5_6_4_a_0_5_10000.tsv" using 2:5 with points title 'WT zlib'
