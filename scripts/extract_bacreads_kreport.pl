#!/usr/bin/env perl

#Extract all reads classified as bacterial/remove all reads classified as bacterial using kraken report
#either from Kraken or Centrifuge

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

my ($outdir,$centr_tab,$centr_ind,$modegz,$verbose);

GetOptions(	"outdir=s" => \$outdir,
			"centrifuge-result=s" => \$centr_tab,
			"index|x=s" => \$centr_ind,
			"gzip" => \$modegz,
			"verbose" => \$verbose,			
			"help" => \&print_help);

#absolute path of project dir
my $wd = dirname($0);
$wd = File::Spec->rel2abs($wd);

if (! -r $centr_tab) {
	warn "ERROR: No centrifuge tabbed input given, please provide via -c\n"; print_help();
}

die "No output directory given\n" unless $outdir;
if (! -d $outdir) {
	mkd("$outdir");
}


my $centi_path;
if ( -x "${wd}/../centrifuge/bin/centrifuge-inspect") {
	$centi_path = "${wd}/../centrifuge/bin/centrifuge-inspect";
} else {
	$centi_path = can_run('centrifuge-inspect')
		or die "'centrifuge-inspect' is not installed, please run tamock.pl --install-deps\n";
	warn "WARNING: Using centrifuge-inspect found at '$centi_path', version might be incompatible, consider running tamock.pl --install-deps\n";
}


##########################################################################
#get all bacterial taxa from index
#load taxonomy tree from centrifuge
#code snippet from https://github.com/infphilo/centrifuge/blob/master/centrifuge-kreport
warn "Loading nodes file ...\n";


my %taxonomy;

open my $IND, "-|", "$centi_path --taxonomy-tree $centr_ind" 
	or die "ERROR: can't open centrifuge index file: $!\n";

while (<$IND>) {
	chomp;
	my @fields = split /\t\|\t/;
	my ($node_id, $parent_id, $rank) = @fields[0,1,2];
	if ($node_id == 1) {
		$parent_id = 0;
	}
	
	#create array for childlist unless already set
	$taxonomy{$parent_id}{childlist} ||= [];

	push @{ $taxonomy{$parent_id}{childlist} }, $node_id;
	$taxonomy{$node_id}{rank} = $rank;
}
close $IND;

#create look up table for all taxa under bacteria
my %bacspecies;
get_childtaxa(2);

##########################################################################
#get all read IDs classified as bacterial from tabbed centrifuge output
warn "Loading read classifications ...\n";
my %bacreads;
my $CT = r_file($centr_tab);
while (my $line= <$CT> ) {
	chomp $line;
	my ($readid,$taxid,$numMatches) = (split("\t",$line))[0,2,7];
	if (exists $bacspecies{$taxid}) {
		$bacreads{$readid} = $numMatches;
	}
}
close $CT;
my $nr_bacreads = keys %bacreads;
warn "'$nr_bacreads' reads classified as bacterial detected\n";

##########################################################################
#read in FASTQ files and separate in bacterial and non bacterial
my @fastqarr;
my $printbac;
warn "Splitting Fastq files ...\n";
foreach my $fastqfile (@ARGV) {
	my $FQIN = r_file($fastqfile);
	my ($fastq_nobac,$fastq_bac,$fpath);
	
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
		$fastq_nobac = "$outdir/${1}nobac.${2}";
		$fastq_bac = "$outdir/${1}bac.${2}";
	} else {
		die "Unexpected filename of inputfile '$fpath', expected <name.suffix(es)>\n";
	}
	my $FQNB = w_file($fastq_nobac);
	my $FQB = w_file($fastq_bac);
	while (my $line = <$FQIN>) {
		#first line is ID of read
		if (($. == 1) || ($.-1)%4 == 0) {
			my $readid = (split(/\s/,$line))[0];
			$readid =~ s/^\@//;
			$printbac = 1 if (exists $bacreads{$readid});
			$fastqarr[0] = $line;
		} elsif (($. == 2) || ($.-2)%4 == 0 ) {
			$fastqarr[1] = $line;
		} elsif (($. == 3) || ($.-3)%4 == 0) {
			$fastqarr[2] = $line;
		} elsif (($. == 4) || ($.-4)%4 == 0) {
			$fastqarr[3] = $line;
			
			if ($printbac) {
				print $FQB @fastqarr;
			} else {
				print $FQNB @fastqarr;
			}
			@fastqarr = ();
			$printbac = undef;
		}
	}
	close $FQIN;
	close $FQNB;
	close $FQB;
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
sub get_childtaxa
{
	my $node = shift;
	my $refchild = $taxonomy{$node}{childlist};
	if ($refchild) {
		for my $child (@$refchild) {
			get_childtaxa($child);
			$bacspecies{$child} = 1;
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
	-x/--index	            	Centrifuge index 
	-o/--outdir         		Output directory
	
	-g/--gzip           		Gzip all output sequence files (off)
	
	-v/--verbose        		Verbose mode, print detailed information to screen
	
	
	-h/--help   		    	Prints this helpmessage
	
EOD
	exit;
}


