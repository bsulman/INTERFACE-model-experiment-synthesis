#!/usr/local/bin/perl

# File: lis2csv.pl
#
# For each .out file in the current directory, create a .csv (comma delimited) file
# The *.out files are not modified
#
# Usage: perl lis2csv.pl
#
# Melannie Hartman
# May 7, 2015
#

@outfiles = glob("*.lis");
foreach $outfile (@outfiles)
{
    print $outfile, "\n";
    $outfile =~ /(\S+)lis$/;
    $csvfile = $1 . "lis.csv";	# $1 is every character before the "lis" at the end of the filename
    print $csvfile, "\n";
    open(OUTFILE,"$outfile") || die "ERROR: Unable to open $outfile\n";
    open(CSVFILE,">$csvfile") || die "ERROR: Unable to open $csvfile\n";
    $lineCnt = 0;
    while ($outline = <OUTFILE>)
    {
        chomp($outline);
        $lineCnt++;
        $outline =~ s/^\s+//;      	# Remove whitespace at the beginning of the line
        if ($lineCnt == 1)
        {
            $outline =~ s/,/\./g;      # Replace all the commas with periods in the hdr line
        }

        $outline =~ s/\s+/,/g;		# Replace remaining whitespace delimiters with a comma
        print CSVFILE "$outline\n";
    }
    close(OUTFILE);
    close(CSVFILE);
}

