#!/usr/bin/env perl

###Usage###
# perl $PathToScript/Scriptname.pl -h

###Author###
#Margo Diricks (mdiricks@fz-borstel.de), with snippets from Christian Utpatel

###Dependecies###
use strict; # disables certain Perl expressions that could behave unexpectedly or are difficult to debug, turning them into errors
use warnings;
use Getopt::Long;
use File::Basename;
use Time::HiRes qw( time ); #Used to calculate elapsed time

###Version###
my $version = "haemophilus_v1.0.0";

###Get time###
my $datestring = localtime();
my $begin_time = time();

###initialize command-line parameter.###
my $opt_version         =      "";
my $species             =      "none";
my $help                =      "";

###get command-line parameter.
GetOptions('version|v'       =>    \$opt_version,
           'species|s:s'     =>    \$species,
           'help|h'          =>     \$help);
print "Species of interest: $species\n";

###print help message if specified or version if specified.
if($help eq '1') { help($version);     exit 0; } #if someone evokes the help function i.e. perl $PathToScript/Scriptname.pl -h
if($opt_version eq '1') { version($version);     exit 0; } # if someone evokes the version version i.e. perl $PathToScript/Scriptname.pl -v


###get files from command-line.
my @kraken_files = @ARGV; #kraken-report files
die "<ERROR>\t",$datestring,"\tNo Kraken-report files specified. Use --help for usage instructions.\n" if ! @kraken_files;
my $PATH_output = dirname($kraken_files[0]);
print "PATH_output is $PATH_output\n";

###Create info file###
open(info_s,">","$PATH_output/info_kraken_summary.txt") or die "\n<ERROR>\t",$datestring,"\tUnable to create kraken_info_summary.txt\n";
print info_s "###Info on kraken_summary.tsv file###\n";
print info_s "Start:\t $datestring\n";
print info_s "Used files located at:\t $PATH_output\n";
print info_s "Making a summary file with $species as species of interest\n";
close (info_s);

###run on all input kraken-report files.
my $krakenfile;

open(OUT,">","$PATH_output/kraken_summary.tsv") or die "\n<ERROR>\t",$datestring,"\tUnable to create kraken_summary.tsv\n";
  my $strain_header = "SampleID\tUnclassified_perc\tUnclassified_reads\tTotal_reads\tHuman_perc\tBacteria_perc\tHaemophilus_perc\tSpecies\tSpecies_perc\tSpecies_other_max\tSpecies_other_max_perc\tSpecies_Haemo_max\tSpecies_Haemo_max_perc\tSpecies_nonHaemo_max\tSpecies_nonHaemo_max_perc\tGenus_max\tGenus_max_perc\tGenus_max_nonHaemo\tGenus_max_nonHaemo_perc\n";
  print OUT $strain_header;
foreach my $file (@kraken_files) {
  next unless (-f "$file");
  next unless ($file =~ /^(.+).report/);
  $krakenfile = $1;
  my $sampleID = basename($krakenfile);
  print "\n<INFO>\t",$datestring,"\tProcessing $sampleID\n";
    
  open(Fin, "<", $file) or die "\n<ERROR>\t",$datestring,"\tUnable to open $file\n";

  my $u_perc = 0; #percentage of reads unclassified
  my $u_reads = 0; #amount of reads unclassified
  my $root_reads = 0; # amount of reads at root (i.e. classified reads)
  my $Total_reads = 0; # amount of total reads (classified + unclassified); 
  my $hum_perc = 0; #percentage human reads
  my $bac_perc = 0; #percentage bacterial reads
  my $haemophilus_perc = 0; #Percentage Haemophilus reads
  my $species_perc = 0; #percentage reads of species of interest
  my $species_max; #name of species with maximum amount of reads and percentage reads
  my $species_max_perc = 0; #perc of most prevalent species other than selected one
  my $species_max_nonHaemo = 0; #name of species with maximum amount of reads and percentage reads that does not belong to genus Haemophilus
  my $species_max_nonHaemo_perc = 0; #perc of most prevalent species that does not belong to genus Haemophilus
  my $species_max_Haemo = 0; #name of species with maximum amount of reads and percentage reads that does not belong to genus Haemophilus
  my $species_max_Haemo_perc = 0; #perc of most prevalent species that does not belong to genus Haemophilus
  my $genus_max = 0; #name of genus with maximum amount of reads
  my $genus_max_perc = 0; #perc of most prevalent genus
  my $genus_max_nonHaemo = 0; #name of genus with 2nd maximum amount of reads
  my $genus_max_nonHaemo_perc = 0; #perc of 2nd most prevalent genus
	
	while (my $line = <Fin>){
	$line          =~  s/\015?\012?$//; #Deal with breaklines
	my @line = split("\t",$line);
	my $perc = $line[0];
	my $read_clade = $line[1];
	my $read_taxon = $line[2];
	my $rank = $line[3];
	my $taxid = $line[4];
	my $name = $line[5];
	$name =~ s/^\s+//; #To remove spaces at the beginning
	
		if ($name eq "unclassified") {
		$u_reads = $read_clade;
		$u_perc = $perc;
		}
		if ($name eq "root") {
		$root_reads = $read_clade;
		}
		$Total_reads = add($u_reads , $root_reads);
		if ($rank eq "G" and $name eq "Haemophilus") {
		$haemophilus_perc = $perc;
		}		
		if ($rank eq "G" and $perc >= $genus_max_perc and $name ne "Homo") {
		$genus_max_perc = $perc;
		$genus_max = $name;
		}
		if ($rank eq "G" and $perc >= $genus_max_nonHaemo_perc and $name !~ /Haemophilus/){
		$genus_max_nonHaemo_perc = $perc;
		$genus_max_nonHaemo = $name;
		} 	
		if ($rank eq "S" and $perc >= $species_max_perc and $name ne "Homo sapiens"and $name ne $species) {
		$species_max_perc = $perc;
		$species_max=$name;
		}
		if ($rank eq "S" and $perc >= $species_max_Haemo_perc and $name ne "Homo sapiens"and $name =~ /Haemophilus/){
		$species_max_Haemo_perc = $perc;
		$species_max_Haemo=$name;
		} 		
		if ($rank eq "S" and $perc >= $species_max_nonHaemo_perc and $name ne "Homo sapiens"and $name !~ /Haemophilus/){
		$species_max_nonHaemo_perc = $perc;
		$species_max_nonHaemo=$name;
		} 		
		if ($rank eq "S" and $name eq $species){
		$species_perc = $perc;
		}
		if ($rank eq "S" and $name=~/^Homo\ssapiens$/){
		$hum_perc = $perc;
		}
		if ($rank eq "D" and $name=~/^Bacteria$/){
		$bac_perc = $perc;
		}
	
}
print OUT "$sampleID\t$u_perc\t$u_reads\t$Total_reads\t$hum_perc\t$bac_perc\t$haemophilus_perc\t$species\t$species_perc\t$species_max\t$species_max_perc\t$species_max_Haemo\t$species_max_Haemo_perc\t$species_max_nonHaemo\t$species_max_nonHaemo_perc\t$genus_max\t$genus_max_perc\t$genus_max_nonHaemo\t$genus_max_nonHaemo_perc\n";
close (Fin);
}
close (OUT);

###Finish info file###
my $end_time = time();
my $elapsed_time = sprintf("%.2f\n", $end_time - $begin_time);

open(info_e,">>","$PATH_output/info_kraken_summary.txt") or die "\n<ERROR>\t",$datestring,"\tUnable to create kraken_info_summary.txt\n";
print info_e "END:\t $datestring\n";
print info_e "Elapsed time (s):\t $elapsed_time";
close (info_e);

exit(0);

sub help { # print a help message.
   my $version = shift;
   print
   "

   kraken_parse_results.pl $version - For help: ask Margo Diricks (mdiricks\@fz-borstel.de)
   
   [USAGE]: [--OPTION PARAMETER] <.report file>
   
   Available OPTIONS:
   -s [--species]       Species of interest (mandatory)
   -h [--help]          This help message
   
   -v [--version]       Print version
   
   ";
   print "\n";
}


sub version { # print the version
   print "$version\n"
}

sub add
{
    my ($x,$y) = @_;
    my $res = $x + $y ;
    return $res ;   
}