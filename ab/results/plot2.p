# http://stackoverflow.com/a/20360148
#
#Several things:
#
#    If you want to output a jpg file, you must use set terminal jpeg first. But in any case I would suggest you to use the pngcairo terminal, if you need a bitmap image.
#
#    The tsv uses tabs as column separator. By default gnuplot uses any white space character as separator, in which case the fifth column is always 2013. So use set datafile separator '\t'.
#
#    In order to have some binning, you must use smooth frequency with an appropriate binning function, which bins your x-values. As y-values I use 1, so that smooth frequency just counts up.
#
#    Possibly you must skip the first line of your data file with every ::1.
#
#    In your case I would use boxes plotting style:

# output as png image
set terminal png
set output 'benchmark2.png'
set datafile separator '\t'
set style fill solid border
set boxwidth 8 absolute

# The graph title
set title "Apache Benchmark"

# Draw gridlines oriented on the y axis
set grid y

# Specify the *input* format of the time data
set timefmt "%s"

# Specify the *output* format for the x-axis tick labels
#set format x "%s"

# Label the x-axis
set xlabel "Response time, ms"

# Label the y-axis
set ylabel "Requests"




set yrange [0:*]
bin(x) = 10*floor(x/10.0)
#plot 'all.tsv' using (bin($5)):(1) every ::1 smooth frequency with boxes title 'ttime'
plot 'all.tsv' using (bin($5)):(1) smooth frequency with boxes notitle

