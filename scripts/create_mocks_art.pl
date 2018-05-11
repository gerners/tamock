#!/usr/bin/env perl

#create mock fastq files from the fullprofile.tsv file resulting form extract_bacterial_reads_from_kraken_report.pl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Copy;
use File::Path qw(make_path remove_tree);
use IO::Compress::Gzip qw(gzip $GzipError) ;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use IPC::Cmd qw(can_run run);
use Time::HiRes qw(gettimeofday tv_interval);

my $t0 = [gettimeofday];

#absolute path of project dir
my $wd = dirname($0);
$wd = File::Spec->rel2abs($wd);

#dependencies
my $art_path;
if ( -x "${wd}/../art/bin/art_illumina") {
	$art_path = "${wd}/../art/bin/art_illumina";
} else {
	$art_path = can_run('art_illumina') 
		or die "'art_illumina' is not installed, please run tamock.pl --install-deps\n";
		warn "WARNING: Using art_illumina found at '$art_path', version might be incompatible, consider running tamock.pl --install-deps\n";
}


my ($profile,$outdir,$refseq_folder,$model,$rlength,$mfl,$sfl,$modegz,$semode);
GetOptions(	"outdir=s" => \$outdir,
			"profile=s" => \$profile,
			"refgenomes|R=s" => \$refseq_folder,
			"length=i" =>\$rlength,
			"mean-fragment-length=i" => \$mfl,
			"sd-fragment-length=i" => \$sfl,
			"illumina-model|M=s" => \$model,
			"gzip" => \$modegz,
			"single-end" => \$semode,
			"help" => \&print_help);

#TODO keep alignemt files option

if (! -d $outdir) {
	mkd("$outdir");
}

my $PR = r_file($profile);
my $ref_seqs_dir = "${outdir}/tmp/mockreads";
if (! -d $ref_seqs_dir) {
	make_path("$ref_seqs_dir") or die "Couldn't create directory '$ref_seqs_dir': $!\n";
}

my $SUM = w_file("${outdir}/art_sampled_read_report.tsv");
print $SUM "Reference Name\tReference Seq ID\tReads/Refgen\tReads/Refseq\tLength/Refgen\tLength/Refseq\tRefgen_file\tRefseq_file\n";

while (my $line = <$PR>) {
	#skip header
	next if $.==1;
	
	my ($abund, $reffile,$refgenlgth) = (split("\t",$line))[0,3,4]; #/^GCF_\d+\.\d_\w+_genomic.fna/
	#check input
	if ($abund !~ /^\d+$/ || $reffile !~ /^GCF_.+_genomic\.fna.gz$/ || $refgenlgth !~ /^\d+$/) {
		die "Bacterial profile file does not have correct format. Required format:\n",
		"<Abundance in int(used)>\t<taxid(unused)>\t<name(unused)>\t<ref file name(used)>\t<ref genome length(used)>\n",
		"Reffile name is expected as in NCBI e.g. 'GCF_001027285.1_ASM102728v1_genomic.fna.gz'\n";
		
	}
	my $refname = substr($reffile,0,-4);
	
	
	#distribute number of reads per genome on respective chromosome/contigs of ref genome
	#split ref genomes if not present
	#print summary of split ref genomes and sampled reads
	 
	simulate_nreads_art($abund,"${refseq_folder}/${reffile}",$refname,$refgenlgth,$outdir,$ref_seqs_dir,$SUM);
}
close $SUM;


#move result file(s) to outdir from tmp
if ($semode) {
	my $filename = "mocked_bac.fq";
	$filename = "${filename}.gz" if $modegz;
	move("${outdir}/tmp/${filename}","${outdir}/${filename}") 
		or die "Couldn't move file '${outdir}/tmp/${filename}' to '${outdir}/${filename}' : $!";
} else {
	for my $i (1,2) {
		my $filename = "mocked_bac_${i}.fq";
		$filename = "${filename}.gz" if $modegz;
		move("${outdir}/tmp/${filename}","${outdir}/${filename}") 
			or die "Couldn't move file '${outdir}/tmp/${filename}' to '${outdir}/${filename}' : $!";
	}
}


#delete tmp dir
rmrf("${outdir}/tmp");


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
sub a_file
{
	my $filepath = shift;
	my $FH;
	
	if ($filepath =~ /\.gz$/i) {
	 	$FH = IO::Compress::Gzip->new("$filepath", Append => 1) or die "ERROR: couldn't append to file '$filepath' : $GzipError\n";
	} else {
		open $FH, ">>$filepath" or die "ERROR: couldn't append to file '$filepath' : $!";
	}
	return $FH;
}
##########################################################################
sub rmrf
{
	my ($goner,$verbose) = @_;
	if (-d $goner) {
		my ($err);
		remove_tree("$goner", {verbose => $verbose, error => \$err});
		if (@$err) {
			for my $diag (@$err) {
				my ($file, $message) = %$diag;
				if ($file eq '') {
					print "Error: $message\n";
				} else {
					print "Problem unlinking $file: $message\n";
				}
			}
		}
	} else {
		unlink $goner or warn "Could not unlink $goner: $!";
	}
}
##########################################################################
sub simulate_nreads_art
{
	my ($abund,$refpath,$refname,$reflgth,$outdir,$ref_seqs_dir,$SUM) = @_;
	my $RF = r_file($refpath);
	my $total_length;
	my %refseqs;
	my %refseqs_sort;
	{
		#change input record separator to read in one sequence "per line" in case of multiline fasta
		local $/=">";
		while(<$RF>) {
			chomp;
			next unless /\w/;
			s/>$//gs;
			my @chunk = split /\n/;
			my $header = shift @chunk;
			my $head1 = (split " ",$header)[0];
			my $seqlen = length join("", @chunk);
			#die "Double seqlength found at $head1 and existing $refseqs{$seqlen}{head}\n" if exists ($refseqs{$seqlen}{head}) ;
			$refseqs{$head1}{seqlen} = $seqlen;
			$refseqs_sort{$seqlen}{$head1} = 1;
			$refseqs{$head1}{seq} = [@chunk];
			$total_length += $seqlen;
		}
		local $/="\n";
		warn "Discrepancy of reference genome '$refname' detected. \nProfile length is '$reflgth' and total length in file is $total_length\n" unless ($reflgth == $total_length);
	}
	close $RF;
	
	my $sum_nreads;
	foreach my $head1 (keys %refseqs) {
		my $nreads = $abund * ($refseqs{$head1}{seqlen}/$total_length);
		#round nreads
		$refseqs{$head1}{nreads} = sprintf "%.0f", $nreads;
		$sum_nreads += $refseqs{$head1}{nreads};
	}
	my $diff = $abund - $sum_nreads;
	
	#correct rounding errors
	while ($diff != 0) {
		foreach my $seqlen (sort {$b <=> $a} keys %refseqs_sort) {
			foreach my $head1 (sort keys %{$refseqs_sort{$seqlen}}) {
				if ($diff > 0) {
					$refseqs{$head1}{nreads}++;
					$diff--;
				} elsif ($diff < 0) {
					$refseqs{$head1}{nreads}--;
					$diff++;
				} else {
					last;
				}
			}
		}
	}
	
	#split reference genome files if not existent, print report and create mockfiles
	foreach my $seqlen (sort {$b <=> $a} keys %refseqs_sort) {
		foreach my $head1 (sort keys %{$refseqs_sort{$seqlen}}) {
			if (! -f "$ref_seqs_dir/${head1}.fna") {
				my $FH = w_file("$ref_seqs_dir/${head1}.fna");
				print $FH ">${head1}\n",join("\n", @{$refseqs{$head1}{seq}}),"\n";
				close $FH;
			}
			print $SUM "$refname\t$head1\t$abund\t$refseqs{$head1}{nreads}\t$total_length\t$seqlen\t$refpath\t$ref_seqs_dir/${head1}.fna\n";
			next unless ($refseqs{$head1}{nreads});
			my $art_cmd = "$art_path -nf 0 -ss ${model} --rcount $refseqs{$head1}{nreads} -l $rlength";
			$art_cmd .= " --paired -m $mfl -s $sfl " unless $semode;
			$art_cmd .= " -i $ref_seqs_dir/${head1}.fna -o ${outdir}/tmp/mockreads/mock.${head1}_";
			chop $art_cmd if $semode; #remove trailing underscore in case of semode
			
			runcmd("$art_cmd",0);
			
			#concat to general base outfile
			if ($semode) {
				my $MR =r_file("${outdir}/tmp/mockreads/mock.${head1}.fq");
				my $outfile = "${outdir}/tmp/mocked_bac.fq";
				$outfile = "${outfile}.gz" if $modegz;
				my $RO = a_file("$outfile");
				while (my $line = <$MR>) {
					print $RO $line;
				}
				close $MR;
				close $RO;
				rmrf("${outdir}/tmp/mockreads/mock.${head1}.fq");
				rmrf("${outdir}/tmp/mockreads/mock.${head1}.aln");
				
			} else {
				my $MR1 =r_file("${outdir}/tmp/mockreads/mock.${head1}_1.fq");
				my $outfile1 = "${outdir}/tmp/mocked_bac_1.fq";
				$outfile1 = "${outfile1}.gz" if $modegz;
				my $RO1 = a_file("$outfile1");
				while (my $line = <$MR1>) {
					print $RO1 $line;
				}
				
				my $MR2 =r_file("${outdir}/tmp/mockreads/mock.${head1}_2.fq");
				my $outfile2 = "${outdir}/tmp/mocked_bac_2.fq";
				$outfile2 = "${outfile2}.gz" if $modegz;
				my $RO2 = a_file("$outfile2");
				while (my $line = <$MR2>) {
					print $RO2 $line;
				}
				for my $fh ($MR1,$RO1,$MR2,$RO2) {
					close $fh or die "Couldn't close filehandle\n";
				}
				for my $i (1,2) {
					rmrf("${outdir}/tmp/mockreads/mock.${head1}_${i}.fq");
					rmrf("${outdir}/tmp/mockreads/mock.${head1}_${i}.aln");
				}
			}
		}
	}
}
##########################################################################
sub runcmd
{
	my ($cmd,$verbose) = @_;
	my $buffer;
	my $success = scalar run (command => "$cmd",verbose => $verbose, buffer => \$buffer);
	if (! $success) {
		die $buffer,"$cmd failed\n" unless $verbose;
		die "'$cmd' failed\n";
	}
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
	
	Create mocked reads for the bacterial fraction of a classified sample 
	
	==== Parameters ====
	
	-p/--profile			Fullprofile table from mock_genomes_from_kraken_report.pl
	-R/--refgenomes			Directory to safe reference genomes in gzipped fasta format
	-x/--index   			Centrifuge index
	-o/--outdir      		Output directory
	
	-g/--gzip       		Gzip all output sequence files (off)
	--single-end     		Create single end mock sequence files instead of paired end (off)
	
	-l/--length	    		Read length of simulated reads (125)
	--mean-fragement-length		Mean size of fragments for paired-end simulations (200)
	--sd-fragment-length		Standard deviation of fragment size for paired-end simulations (19)
	-M/--illumina-model		Illumima error model for ART (HS25)
	                        	For available profiles check '$art_path -h'
	
	-v/--verbose			Verbose mode, print detailed information to screen
	
	
	-h/--help       		Prints this helpmessage
	
EOD
	exit;
}
