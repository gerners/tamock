#!/usr/bin/env perl

#Extract all reads with taxa assigned to a reference genome (from taxa_w_refgenome.tsv)

use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path remove_tree);
use Getopt::Long;
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use IPC::Cmd qw(can_run run);
use Time::HiRes qw(gettimeofday tv_interval);

my $t0 = [gettimeofday];

my ($outdir,$centr_tab,$tax2brep,$modegz,$verbose);

GetOptions(	"outdir=s" => \$outdir,
			"taxa-refgenomes|x" => $tax2brep,
			"centrifuge-result=s" => \$centr_tab,
			"gzip" => \$modegz,
			"verbose" => \$verbose,			
			"help" => \&print_help);

#absolute path of project dir
my $wd = dirname($0);
$wd = File::Spec->rel2abs($wd);

if (! -r $centr_tab) {
	warn "ERROR: No centrifuge tabbed input given, please provide via -c\n"; print_help();
}

if (! -r $tax2brep) {
	warn "ERROR: No table containing taxa with reference genomes given, please provide via -x\n"; print_help();
}

die "No output directory given\n" unless $outdir;
if (! -d $outdir) {
	mkd("$outdir");
}

##########################################################################
#get all bacterial taxa from temporary file classification_out/taxa_w_refgenome.tsv (mock_profile_kreport.pl)
my %taxa_repl;
my $T2BR = r_file("$tax2brep");

while (my $line = <$T2BR>) {
	#skip header
	next if $. == 1;
	my @larr = split("\t",$line);
	$taxa_repl{$larr[0]} = $larr[1];
}
close $T2BR;

##########################################################################
#get all read IDs classified as bacterial from tabbed centrifuge output
warn "Loading read classifications ...\n";
my %bacreads;
my $CT = r_file($centr_tab);
while (my $line= <$CT> ) {
	chomp $line;
	my ($readid,$taxid,$numMatches) = (split("\t",$line))[0,2,7];
	if (exists $taxa_repl{$taxid}) {
		$bacreads{$readid} = $numMatches;
	}
}
close $CT;
my $nr_bacreads = keys %bacreads;
warn "'$nr_bacreads' reads will be replaced with simulated sequences\n";

##########################################################################
#read in FASTQ files and separate in bacterial and non bacterial
my @fastqarr;
my $printbac;
warn "Splitting Fastq files ...\n";
foreach my $fastqfile (@ARGV) {
	my $FQIN = r_file($fastqfile);
	my ($fastq_noref,$fastq_ref,$fpath);
	
	#adapt filepath for gzip/non-gzip output
	
	if ($modegz && $fastqfile !~ /\.gz$/ ) {
		$fpath = "${fastqfile}.gz";
	} elsif ($modegz && $fastqfile =~ /\.gz$/) {
		$fpath = $fastqfile;
	} elsif (! $modegz && $fastqfile !~ /\.gz$/ ) {
		$fpath = $fastqfile;
	} elsif (! $modegz && $fastqfile =~ /\.gz$/) {
		$fpath = substr($fastqfile,0,-3);
	} else {
		die "Unexpected filename of inputfile '$fpath', expected <name.suffix(es)>\n";
	}
	
	$fpath = basename($fpath);
	if ($fpath =~ /(\w+\.)(.+)/) {
		$fastq_noref = "$outdir/${1}noref.${2}";
		$fastq_ref = "$outdir/${1}ref.${2}";
	} else {
		die "Unexpected filename of inputfile '$fpath', expected <name.suffix(es)>\n";
	}
	my $FQNR = w_file($fastq_noref);
	my $FQR = w_file($fastq_ref);
	while (my $line = <$FQIN>) {
		#first line is ID of read
		if (($. == 1) || ($.-1)%4 == 0) {
			my $readid = (split(/\s/,$line))[0];
			$readid =~ s/^\@//;
			$readid =~ s#\/\d$##;
			$printbac = 1 if (exists $bacreads{$readid});
			$fastqarr[0] = $line;
		} elsif (($. == 2) || ($.-2)%4 == 0 ) {
			$fastqarr[1] = $line;
		} elsif (($. == 3) || ($.-3)%4 == 0) {
			$fastqarr[2] = $line;
		} elsif (($. == 4) || ($.-4)%4 == 0) {
			$fastqarr[3] = $line;
			
			if ($printbac) {
				print $FQR @fastqarr;
			} else {
				print $FQNR @fastqarr;
			}
			@fastqarr = ();
			$printbac = undef;
		}
	}
	close $FQIN;
	close $FQNR;
	close $FQR;
}

print "$0 took ",runtime(tv_interval($t0)), " to run\n";

##########################################################################
#subroutines
##########################################################################
sub r_file
{
	my $filepath = shift;
	my $FH;
	if ($filepath =~ /\.gz$/i) {
	 	$FH = IO::Uncompress::Gunzip->new("$filepath") or die "ERROR: couldn't open file '$filepath' : $GunzipError\n";
	} else {
		open $FH, "$filepath" or die "ERROR: couldn't open file '$filepath' : $!";
	}
	return $FH;
}
##########################################################################
sub w_file
{
	my $filepath = shift;
	my $FH;
	
	if ($filepath =~ /\.gz$/i) {
	 	$FH = IO::Compress::Gzip->new("$filepath") or die "ERROR: couldn't write to file '$filepath' : $GzipError\n";
	} else {
		open $FH, ">$filepath" or die "ERROR: couldn't write to file '$filepath' : $!";
	}
	return $FH;
}
##########################################################################
sub mkd
{
	my ($dir,$verbose) = @_,
	my $err;
	make_path("$dir",{verbose => $verbose, error => \$err});
	if (@$err) {
		for my $diag (@$err) {
			my ($file, $message) = %$diag;
			if ($file eq '') {
				print "General error: $message\n";
			} else {
				print "Problem creating $file: $message\n";
			}
		}
	}
}
##########################################################################
sub runtime
{
	#get runtime in seconds from e.g. "(time - $^T)"
	my $time = shift @_;
	my $rtime;
	#check if script ran for more than one minute
	if ($time > 59 ) {
		#or more than one hour
		if ($time > 3599 ) {
			#or more than one day
			if ($time > 86399) {
				$rtime = int($time / 86400) . "d " . int(($time % 86400) / 3600) . "h " . int((($time % 86400) % 3600) / 60) . "m " . (((time % 86400) % 3600) % 60) . "s";
				return $rtime;
			}
			$rtime = int($time / 3600) . "h " . int(($time % 3600) / 60) . "m " . (($time % 3600) % 60) . "s";
			return $rtime;	
		}
		$rtime = int($time / 60) . "m " . ($time % 60) . "s";
		return $rtime;
	}
	$rtime = $time . "s";
	return $rtime;
}
##########################################################################
sub print_help
{
	print STDERR <<EOD;

Usage: $0 [Parameters] <fastqfile1.fq> <fastfile2.fq> ...
	
	Extract classified bacterial sequences using centrifuge index
	
	==== Parameters ====
	
	-c/--centrifuge-result		Centrifuge classification outout in tab format
	-x/--taxa-refgenomes		Table containing all taxa with reference genomes assigned (taxa_w_refgenome.tsv)
	-o/--outdir         		Output directory
	
	-g/--gzip           		Gzip all output sequence files (off)
	
	-v/--verbose        		Verbose mode, print detailed information to screen
	
	
	-h/--help   		    	Prints this helpmessage
	
EOD
	exit;
}


