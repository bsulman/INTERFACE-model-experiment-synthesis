#!/usr/local/bin/perl

# File: out2csv.pl
#
# For each .out file in the current directory, create a .csv (comma delimited) file
# The *.out files are not modified
#
# Usage: perl out2csv.pl
#
# Melannie Hartman
# May 7, 2015
#

@outfiles = glob("*.out");
foreach $outfile (@outfiles)
{
    print $outfile, "\n";
    $outfile =~ /(\S+)out$/;
    $csvfile = $1 . "out.csv";	# $1 is every character before the "out" at the end of the filename
    print $csvfile, "\n";
    open(OUTFILE,"$outfile") || die "ERROR: Unable to open $outfile\n";
    open(CSVFILE,">$csvfile") || die "ERROR: Unable to open $csvfile\n";
    while ($outline = <OUTFILE>)
    {
        chomp($outline);
        $outline =~ s/^\s+//;      	# Remove whitespace at the beginning of the line
        $outline =~ s/\s+/,/g;		# Replace remaining whitespace delimiters with a comma
        print CSVFILE "$outline\n";
    }
    close(OUTFILE);
    close(CSVFILE);
}

