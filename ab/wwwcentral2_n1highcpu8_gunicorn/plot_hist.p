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
set terminal pngcairo size 1000,800
#set output 'ab_hist_storage_ElCurveBrowse.png'
set output 'ab_hist_storage_G2CurveQstats.png'
set datafile separator '\t'

# The graph title
set title "Apache Benchmark"
set label "www-central2.lmfdb.xyz/Genus2Curve/Q/stats" at graph 0.05, 0.95

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

#set style fill transparent solid noborder
set style fill transparent solid 0.3 noborder

set yrange [0:*]
set xrange [100:200]
bin(x,bwidth) = bwidth*floor(x/bwidth)
bwidth=1
plot \
"n1hcpu4_mmap__Genus2Curve_Q_stats5_10000.tsv" using (bin($5,bwidth)):(1) every ::1 smooth frequency with boxes title 'MMAP', \
"n1hcpu4_wt__Genus2Curve_Q_stats5_10000.tsv" using (bin($5,bwidth)):(1) every ::1 smooth frequency with boxes title 'WT', \
"n1hcpu4_wtzl__Genus2Curve_Q_stats5_10000.tsv" using (bin($5,bwidth)):(1) every ::1 smooth frequency with boxes title 'WT zlib'
#"n1hcpu4_mmap__EllipticCurve_browse5_10000.tsv" using (bin($5,bwidth)):(1) every ::1 smooth frequency with boxes title 'MMAP', \
#"n1hcpu4_wtzl__EllipticCurve_browse5_10000.tsv" using (bin($5,bwidth)):(1) every ::1 smooth frequency with boxes title 'WT zlib', \
#"n1hcpu4_wt__EllipticCurve_browse5_10000.tsv" using (bin($5,bwidth)):(1) every ::1 smooth frequency with boxes title 'WT'
