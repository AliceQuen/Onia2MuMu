#!/bin/bash
rm -f report.out
echo '-------------------------** Crab Report **-------------------------' > report.out
for i in crab_crab3_noPriFit*
do
	sed -i -e '$a\'"${i}" report.out
	crab status $i | sed -n -e '/running/p' -e '/jobs failed/p' -e '/finished/p' > tmp_report.out
	cat tmp_report.out | while read rows
	do
		sed -i -e '$a\'"${rows}" report.out	
	done 
	rm -f tmp_report.out
done

echo '-------------------------** Overall Summary **-------------------------' >> report.out

# Calculate and append overall finished percentage
grep 'finished' report.out | sed 's/.*(//; s/).*//' | awk -F/ '
{
    sum_finished += $1;
    sum_total += $2;
}
END {
    if (sum_total > 0) {
        percentage = (sum_finished / sum_total) * 100;
        printf "Overall Finished: %.2f%% (%d/%d)\n", percentage, sum_finished, sum_total;
    } else {
        print "No finished jobs found to calculate percentage.";
    }
}' >> report.out
