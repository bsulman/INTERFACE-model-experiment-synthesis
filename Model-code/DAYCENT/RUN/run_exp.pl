#!/usr/local/bin/perl
#
#############################################################################
#  FILE:   run.pl
#
#  AUTHOR: Melannie Hartman
#          November 13, 2011
#          April 22, 2016
#          March 20, 2017
#
#  PURPOSE:
#
#  Run DayCent for each litter type and site.
#
#  Example command lines with litter = "ACSAf" and site = "ARC"
#    ./DailyDayCent_srad_SOM -s experiment1.ACSAf -n experiment2.ACSAf.ARC -c ARC.cdi.txt
#    mv decomp.csv ACSAf.ARC.decomp.csv
#
#############################################################################

$DayCent = "./DailyDayCent_srad_SOM";
$List100 = "./DailyDayCent_srad_list100";

@litters = ("HQ500", "LQ500");
@soils = ("high_protect", "med_protect", "low_protect");
@dtemps = ("20C", "22C", "25C");

# Increase temperature from 20C used in spinup to 22C and 25C 
#   Use experiment_22C.cdi.txt and experiment_25C.cdi.txt

foreach $soil (@soils)
{
    foreach $litter (@litters)
    {
        foreach $dtemp (@dtemps)
        {
            $schfile = "spin_${litter}";
            $ebinfile = "spin_${litter}.${soil}";
    
            $command_string = sprintf "cp soils_%s.in soils.in", $soil;
            printf "$command_string\n";
            system(${command_string});
    
            #Extend and change CDI file to simulate warming of +5C
    
            $cdifile = "experiment_${dtemp}.cdi.txt";
            $schfile = "exp_${litter}";
            $binfile2 = "exp_${litter}.${soil}.${dtemp}";
            unlink("${binfile2}.bin");
            unlink("${binfile2}.lis");
    
            $command_string = sprintf "%s  -s %s  -n %s -e %s -c %s", $DayCent, $schfile, $binfile2, $ebinfile, $cdifile;
            printf "$command_string\n";
            system(${command_string});
    
            $command_string = sprintf "%s  %s  %s  %s", $List100, $binfile2, $binfile2, "outvars.txt";
            printf "$command_string\n";
            system(${command_string});
    
            $outputfile = "exp_${litter}.${soil}.${dtemp}.decomp.csv";
            $command_string = sprintf "mv decomp.csv %s", $outputfile;
            printf "$command_string\n";
            system(${command_string});
        }
    }
}

# Change litter input amounts and quality, no temperature increase (experiment_20C.cdi.txt)
@litters = ("HQX2", "LQX2", "HQ650", "LQ650","HQL30", "LQL30", "LQ0", "HQ0");
@soils = ("high_protect", "med_protect", "low_protect");
@dtemps = ("20C", "22C");

foreach $soil (@soils)
{
    foreach $litter (@litters)
    {

        foreach $dtemp (@dtemps)
        {
            $schfile = "spin_${litter}";
            # Determine the correct spinup binary file to extend from
            if (index($litter, "HQ") != -1) 
            {
                $ebinfile = "spin_HQ500.${soil}";
            }
            else
            {
                $ebinfile = "spin_LQ500.${soil}";
            }
    
            $command_string = sprintf "cp soils_%s.in soils.in", $soil;
            printf "$command_string\n";
            system(${command_string});
    
            #Extend and change CDI file to simulate warming of +5C
    
            $cdifile = "experiment_${dtemp}.cdi.txt";
            $schfile = "exp_${litter}";
            $binfile2 = "exp_${litter}.${soil}.${dtemp}";
            unlink("${binfile2}.bin");
            unlink("${binfile2}.lis");
    
            $command_string = sprintf "%s  -s %s  -n %s -e %s -c %s", $DayCent, $schfile, $binfile2, $ebinfile, $cdifile;
            printf "$command_string\n";
            system(${command_string});
    
            $command_string = sprintf "%s  %s  %s  %s", $List100, $binfile2, $binfile2, "outvars.txt";
            printf "$command_string\n";
            system(${command_string});
    
            $outputfile = "exp_${litter}.${soil}.${dtemp}.decomp.csv";
            $command_string = sprintf "mv decomp.csv %s", $outputfile;
            printf "$command_string\n";
            system(${command_string});
    
            $outputfile = "exp_${litter}.${soil}.${dtemp}.cbalance.csv";
            $command_string = sprintf "mv cbalance.csv %s", $outputfile;
            printf "$command_string\n";
            system(${command_string});
    
            $outputfile = "exp_${litter}.${soil}.${dtemp}.nbalance.csv";
            $command_string = sprintf "mv nbalance.csv %s", $outputfile;
            printf "$command_string\n";
            system(${command_string});

        }

    }
}




