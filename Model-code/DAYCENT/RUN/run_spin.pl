#!/usr/local/bin/perl
#
#############################################################################
#  FILE:   run.pl
#
#  AUTHOR: Melannie Hartman
#          November 13, 2011
#          April 22, 2016
#          May 4, 2016
#          March 20, 2017
#
#  PURPOSE:
#
#  Run DayCent spinup for each litter type and site.
#
#  Example command lines with litter = "ACSAf" and site = "ARC"
#    ./DailyDayCent_srad_SOM -s experiment1.ACSAf -n experiment2.ACSAf.ARC -c ARC.cdi.txt
#    mv decomp.csv ACSAf.ARC.decomp.csv
#
#############################################################################

$DayCent = "./DailyDayCent_srad_SOM";

@litters = ("HQ500", "LQ500");
@soils = ("high_protect", "med_protect", "low_protect");


foreach $soil (@soils)
{
    foreach $litter (@litters)
    {
        $schfile = "spin_${litter}";
        $binfile = "spin_${litter}.${soil}";
        unlink("${binfile}.bin");
        $cdifile = "experiment_20C.cdi.txt";

        $command_string = sprintf "cp soils_%s.in soils.in", $soil;
        printf "$command_string\n";
        system(${command_string});

        $command_string = sprintf "%s  -s %s  -n %s -c %s", $DayCent, $schfile, $binfile, $cdifile;
        printf "$command_string\n";
        system(${command_string});

        $outputfile = "spin_${litter}.${soil}.decomp.csv";
        $command_string = sprintf "mv decomp.csv %s", $outputfile;
        printf "$command_string\n";
        system(${command_string});

        $outputfile = "spin_${litter}.${soil}.diagnostics.txt";
        $command_string = sprintf "mv diagnostics.txt %s", $outputfile;
        printf "$command_string\n";
        system(${command_string});

        $command_string = sprintf "%s  %s  %s  %s", $List100, $binfile, $binfile, "outvars-lidet.txt";
        printf "$command_string\n";
        system(${command_string});
    }
}



