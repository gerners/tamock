#!/usr/bin/env perl

#select representative genomes from kraken_report of a sample classified by Kraken or Centrifuge
#and download/check corresponding genomes from NCBI Refseq

#if multiple reference genomes are found for one taxid, completion of assembly and then 
#date of release is tested for selection
#remaining cases are selected random by first occurence

#Rules
#	a)	Reads only classified at species level will be distributed to all strains with assigned reads of same species
#		using the same ratio as already assigned reads to respective strains
#	aa)	If lowest assignment of reads is on species level, reference strain according to NCBI assembly summary 
#		of said species will be sampled
#	b)	Reads of strains without a reference genome will be assigned to species level
#		and then reassigned to present strains or reference strains according to a) except no-reassign is set
#	c)	Species lvl assigned reads are rounded using the ratios of reads assigned to each strain

use strict;
use warnings;
use File::Basename;
use File::Copy;
use File::Path qw(make_path remove_tree);
use File::Fetch;
use Getopt::Long;
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Cwd;
use Time::HiRes qw(gettimeofday tv_interval);

my $t0 = [gettimeofday];

##########################################################################
#option handling
my ($kraken_report,$assembly_summary,$refseq_folder,$no_strain_ra,$debug);

#set defaults
my $outdir = getcwd;
my $verbose = 0;

GetOptions(	"kraken-report=s" => \$kraken_report,
			"assembly-summary=s" => \$assembly_summary,
			"refgenomes|R=s" => \$refseq_folder,
			"outdir=s" => \$outdir,
			"no-reassign" => \$no_strain_ra,
			"verbose+" => \$verbose,
			"debug+" => \$debug,
			"help" => \&print_help);
$verbose = 3;
#check mandatory options
if (! -f $kraken_report) {
	warn "ERROR: No Kraken input given/invalid file, please provide via -k\n"; print_help(1);
} elsif (! -f $assembly_summary) {
	warn "ERROR: No assembly summary input given, please provide via -a\n"; print_help(1);
} elsif (! $refseq_folder) {
	warn "ERROR: No refseq directory given, please provide via -r\n"; print_help();
}
if (! -d $refseq_folder) {
	warn "WARNING: Refseq directory not existing, creating directory '$refseq_folder'\n";
	mkd("$refseq_folder");
}

if (! -d $outdir) {
	mkd("$outdir");
}

##########################################################################
#read in kraken report file

#save level of line by counting leading whitespaces of Name field from Kraken report
my $sp_lvl = 0;
my $g_lvl = 0;
my ($prev_taxid,$sp_taxid,$prev_href);
my %species;

#flag to only read in bacterial species as only RefSeq for Bacteria are checked
my $bacflag = 0;

my $total_reads_nonspecieslvl;

my ($DB1,$DB2);
if ($debug) {
	open $DB1, ">$outdir/debug1" or die "ERROR: couldn't open file '$outdir/debug1' : $!";
	open $DB2, ">$outdir/debug2" or die "ERROR: couldn't open file '$outdir/debug2' : $!";
}

my $KR = r_file("$kraken_report");
while (my $line = <$KR>) {
	chomp $line;
	#skip header
	next if ($. == 1 && $line =~ /^Percentage/);
		
	my @larr = split("\t",$line);
	my $level = () = $larr[5] =~ /\G\s/g;
	
	my ($r_read, $r_ass, $rank, $taxid, $name) = @larr[1,2,3,4,5];
	$name =~ s/^\s+//;
	
	#save nr of unclassified, root and bacterial reads
	if ($. == 1 && $rank eq "U") {
		assign_species($taxid,$r_read,$r_ass,$rank,$.,$name);
		print $DB1 "$r_read\t$r_ass\t$rank\t$level\t$taxid\t$.\t$name\n" if $debug;
		
	} elsif ($. == 2 && $name eq "root") {
		assign_species($taxid,$r_read,$r_ass,$rank,$.,$name);
		print $DB1 "$r_read\t$r_ass\t$rank\t$level\t$taxid\t$.\t$name\n" if $debug;
		
	} elsif ($rank eq "D" && $name eq "Bacteria") {
		assign_species($taxid,$r_read,$r_ass,$rank,$.,$name);
		print $DB1 "$r_read\t$r_ass\t$rank\t$level\t$taxid\t$.\t$name\n" if $debug;
	}
	
	#check if bacterial line is reached or any other domain after that and skip before and afterwards
	if ($bacflag) {
		$bacflag = 0 if ($rank eq "D");
		warn "INFO: Last bacterial taxa found at linenr $.\n" if ($rank eq "D" && $verbose > 1);
	} else {
		$bacflag = 1 if ($rank eq "D" && $taxid == 2 && $name eq "Bacteria");
		warn "INFO: First bacterial taxa found at linenr" . ($. - 1) . "\n" if ($rank eq "D" && $taxid == 2 && $name eq "Bacteria" && $verbose > 1);
		next;
	}
	next unless ($bacflag == 1);
	
	if ($rank eq "G") {
		$g_lvl = $level;
		$sp_lvl = 0;
		print "$r_read\t$r_ass\t$rank\t$level\t$taxid\t$.\t$name\n" if $verbose > 2;
		$total_reads_nonspecieslvl += $r_ass if ($r_ass);
	
	} elsif ($rank =~ /^[PCOF]$/) {
		$sp_lvl = 0;
		print "$r_read\t$r_ass\t$rank\t$level\t$taxid\t$.\t$name\n" if $verbose > 2;
		$total_reads_nonspecieslvl += $r_ass if ($r_ass);
		
	} elsif ($level > $sp_lvl && $sp_lvl > 0) {
		die "Unexpected line below species line found at line nr '$.', level '$level', sp_level '$sp_lvl' and taxid $taxid/prev_taxid $prev_href->{taxid}\n" unless ($rank =~ /[-S]/);
		
		print $DB2 "$r_read\t$r_ass\t$rank\t$level\t$taxid\t$.\t$name\n" if $debug;
		
		if ($level == $prev_href->{level}) {
			
			#reset prev* to parent of previous since same level is maintained
			$prev_href = $species{$prev_href->{taxid}}{strainof};
			
			assign_strain($prev_href,$taxid,$r_read,$r_ass,$level,$.,$name);
			
			#set prev_href to current
			$prev_href = $prev_href->{strains}{$taxid};
			
		} elsif ($level > $prev_href->{level}) {
			
			assign_strain($prev_href,$taxid,$r_read,$r_ass,$level,$.,$name);
			
			#set prev* to current
			$prev_href = $prev_href->{strains}{$taxid};
			
		} elsif ($level < $prev_href->{level} ) {
			
			#set prev* to parent of parent, potentially going up more than one parent
			my $nr_parent = ($prev_href->{level} - $level)/2;
			#for one parent, jump up two, two parents up, three jumps up
			for (my $i = 0; $i <= $nr_parent; $i++) {
				$prev_href = $species{$prev_href->{taxid}}{strainof};
			}
			
			assign_strain($prev_href,$taxid,$r_read,$r_ass,$level,$.,$name);
			
			#set prev* to current level since current entry is not a child of previous
			$prev_href = $prev_href->{strains}{$taxid};
			
		} else {
			die "BUG: Error when determing taxonomic species/strain level at line $. with level '$level'/sp_level '$sp_lvl' and prev_level '$prev_href->{level}'\n";
		}
			
	#first level below genus will be considered top species level to encompass subgroups or species group is equal or higher to previous genus line
	} elsif ($level == ($g_lvl - 2) || $level >= $g_lvl && $rank eq "S") {
		$sp_lvl = $level;
		
		#save species with taxid as key
		assign_species($taxid,$r_read,$r_ass,$level,$.,$name);
		print $DB1 "$r_read\t$r_ass\t$rank\t$level\t$taxid\t$.\t$name\n" if $debug;
		
		$prev_href = $species{$taxid};
		
	} else {
		print $DB2 "$r_read\t$r_ass\t$rank\t$level\t$taxid\t$.\t$name\n" if $debug;;
		$total_reads_nonspecieslvl += $r_ass if ($r_ass);
	}
}
close $KR;
#########################################################
#create list with all toplevel species for later iterations
my @toplvl_species;
foreach my $taxid (sort {$a <=> $b} keys %species) {
	next if $taxid < 3;
	if (! $species{$taxid}{strainof}) {
		push @toplvl_species,$taxid;
	}
}

#################################################
#read in assembly summary from NCBI and download/check all refseq genome files for species/strains present in kraken output file 

my %genomes;
my %refstrains;

my $NCBI = r_file("$assembly_summary");

while (my $line = <$NCBI>) {
	next if ($line =~/^#/ );
	chomp $line;
	my @larr = split("\t",$line);
	my ($acc, $refcat, $taxid, $s_taxid, $name, $ass_lvl, $rel_date, $ftp) = @larr[0,4,5,6,7,11,14,19];
	$rel_date =~ s+\/++g;
	
	#check level of assembly and assign it a value between 1-4 for comparison
	if ($ass_lvl eq "Complete Genome") {
		$ass_lvl = 1;
	} elsif ($ass_lvl eq "Chromosome") {
		$ass_lvl = 2;
	} elsif ($ass_lvl eq "Scaffold") {
		$ass_lvl = 3;
	} elsif ($ass_lvl eq "Contig") {
		$ass_lvl = 4;
	} else {
		die "Unknown assembly level found at line '$.' and level '$ass_lvl'\n";
	}
	
	###################
	#check if another reference is already present for current ref genome, if so, test for better completeness 
	#via assembly level or most resent release date
	
	#SPECIES/STRAIN present in krakenstyle report
	if ($species{$taxid}) {
		
		#SPECIES/STRAIN & REFERENCE GENOME
		#always save reference genome
		if ($refcat eq "reference genome") { 
			
			assign_genome($taxid,$acc,$s_taxid,$name,$ass_lvl,$ftp,$rel_date,$refcat);
			$refstrains{$s_taxid} = $taxid;
			
		#SPECIES/STRAIN & NONREF GENOME
		} else {
			
			#if multiple reference genomes with the same taxid are present, rank in the order of
			#reference genome > representative genome > assembly status > date
			
			#strain/species has already a reference assigned -> check if current refseq genome should replace present entry
			if ($genomes{$taxid}) {
				
				next if ($genomes{$taxid}{category} eq "reference genome" || 
				($genomes{$taxid}{category} eq "representative genome" && $refcat eq "na") );
				
				if ($refcat eq "representative genome" && $genomes{$taxid}{category} eq "na") {
					assign_genome($taxid,$acc,$s_taxid,$name,$ass_lvl,$ftp,$rel_date,$refcat);
				} elsif ($ass_lvl < $genomes{$taxid}{assembly}) {
					assign_genome($taxid,$acc,$s_taxid,$name,$ass_lvl,$ftp,$rel_date,$refcat);
				} elsif ($ass_lvl == $genomes{$taxid}{assembly} && $rel_date > $genomes{$taxid}{release}) {
					assign_genome($taxid,$acc,$s_taxid,$name,$ass_lvl,$ftp,$rel_date,$refcat);
				} else {
					next if ($ass_lvl > $genomes{$taxid}{assembly} || $rel_date < $genomes{$taxid}{release});
					warn "INFO: Multiple reference genomes for same strain taxid and date/completeness found.\n",
					"\tReference genome with accession '$genomes{$taxid}{accession}' is kept while '$acc' at line '$.' is dropped\n" if ($verbose > 1);
					next;
				}
								
			#strain/species has no reference assigned -> assign current refseq genome as reference
			} else {
				assign_genome($taxid,$acc,$s_taxid,$name,$ass_lvl,$ftp,$rel_date,$refcat);
			}
			
			#check if current strain/species should be used as a reference genome on species level
			if ($refstrains{$s_taxid}) {
				
				next if ($genomes{$refstrains{$s_taxid}}{category} eq "reference genome" || 
				($genomes{$refstrains{$s_taxid}}{category} eq "representative genome" && $refcat eq "na") );
				
				if ($refcat eq "representative genome" && $genomes{$refstrains{$s_taxid}}{category} eq "na") {
					$refstrains{$s_taxid} = $taxid;
				} elsif ($ass_lvl < $genomes{$refstrains{$s_taxid}}{assembly}) {
					$refstrains{$s_taxid} = $taxid;
				} elsif ($ass_lvl == $genomes{$refstrains{$s_taxid}}{assembly} && $rel_date > $genomes{$refstrains{$s_taxid}}{release}) {
					$refstrains{$s_taxid} = $taxid;
				} else {
					next if ($ass_lvl > $genomes{$refstrains{$s_taxid}}{assembly} || $rel_date < $genomes{$refstrains{$s_taxid}}{release});
					warn "INFO: Multiple reference genomes for same species taxid and date/completeness found.\n",
					"\tReference genome with accession '$genomes{$refstrains{$s_taxid}}{accession}' is kept while '$acc' at line '$.' is dropped\n" if ($verbose > 1);
					next;
				}
				
			} else {
				$refstrains{$s_taxid} = $taxid;
			}
		}
	
	#STRAIN not present in Kraken report
	#save reference strain genomes even if strain is not present in data in case of reads only classified on species level 
	} elsif ($species{$s_taxid}) {
		
		assign_genome($taxid,$acc,$s_taxid,$name,$ass_lvl,$ftp,$rel_date,$refcat);
		
		#check if current strain/species should be used as a reference genome on species level
		if ($refstrains{$s_taxid}) {
			
			next if ($genomes{$refstrains{$s_taxid}}{category} eq "reference genome" || 
			($genomes{$refstrains{$s_taxid}}{category} eq "representative genome" && $refcat eq "na") );
			
			if ($refcat eq "representative genome" && $genomes{$refstrains{$s_taxid}}{category} eq "na") {
				$refstrains{$s_taxid} = $taxid;
			} elsif ($ass_lvl < $genomes{$refstrains{$s_taxid}}{assembly}) {
				$refstrains{$s_taxid} = $taxid;
			} elsif ($ass_lvl == $genomes{$refstrains{$s_taxid}}{assembly} && $rel_date > $genomes{$refstrains{$s_taxid}}{release}) {
				$refstrains{$s_taxid} = $taxid;
			} else {
				next if ($ass_lvl > $genomes{$refstrains{$s_taxid}}{assembly} || $rel_date < $genomes{$refstrains{$s_taxid}}{release});
				warn "INFO: Multiple reference genomes for same species taxid and date/completeness found.\n",
				"\tReference genome with accession '$genomes{$refstrains{$s_taxid}}{accession}' is kept while '$acc' at line '$.' is dropped\n" if ($verbose > 1);
				next;
			}
			
		} else {
			$refstrains{$s_taxid} = $taxid;
		}
		
		#assign detected reference strain for species
		if ($species{$s_taxid}{strainof}) {
			$prev_href = $species{$s_taxid}{strainof};
			$prev_href = $prev_href->{strains}{$s_taxid}
		} else {
			$prev_href = $species{$s_taxid};
		}
		assign_strain($prev_href,$taxid,0,0,undef,"NCBI-$.",$name);
	
	#STRAIN/SPECIES NOT PRESENT
	#species of current strain was not detected in Kraken output, therefore skip entry
	} else {
		next;
	}
}
close $NCBI;

#########################################################



#assign reference genomes to all species without assigned reference genomes, which have refgenomes from NCBI strains with zero reads assigned
my @missing_refs;
foreach my $taxid (keys %species) {
	my $href;
	#skip all entries without child strains or already a reference assigned
	if ($species{$taxid}{strainof}) {
		next if (! $species{$taxid}{strainof}->{strains}{$taxid}{strains} || $refstrains{$taxid});
		$href = $species{$taxid}{strainof}->{strains}{$taxid};
	} else {
		next if $refstrains{$taxid};
		$href = $species{$taxid};
	}
	
	#look for genome if no strains are present and directly assign refgenome
	if (scalar keys %{$href->{strains}} < 1) {
		if ($genomes{$taxid}) {
			$refstrains{$taxid} = $taxid;
			assign_strain($href,$taxid,0,0,undef,undef,$genomes{$taxid}{organism});
		}
	} else {
		push @missing_refs,$href;
	}
}

foreach my $href (@missing_refs) {
	my $sp_taxid = $href->{taxid};
	
	#all taxid's have to have strains present, otherwise they could not have ended up in @missing_refs
	foreach my $st_taxid ( keys %{$href->{strains}}) {
		
		#reference strain might have been set from a previous strain in this loop although actual reference strain coming around
		if ($refstrains{$sp_taxid}) {
			
			#skip strains/species without a reference genome present or where a reference genome is set
			next if (! exists $genomes{$st_taxid});
			next if ($genomes{$refstrains{$sp_taxid}}{category} eq "reference genome" || 
					($genomes{$refstrains{$sp_taxid}}{category} eq "representative genome" && $genomes{$st_taxid}{category} eq "na"));
			
			if ($genomes{$st_taxid}{category} eq "reference genome") { 
				$refstrains{$sp_taxid} = $st_taxid;
			} elsif ($genomes{$st_taxid}{category} eq "representative genome" && $genomes{$refstrains{$sp_taxid}}{category} eq "na" ) {
				$refstrains{$sp_taxid} = $st_taxid;
			} elsif ($genomes{$st_taxid}{assembly} < $genomes{$refstrains{$sp_taxid}}{assembly}) {
				$refstrains{$sp_taxid} = $st_taxid;
			} elsif ($genomes{$st_taxid}{assembly} == $genomes{$refstrains{$sp_taxid}}{assembly} && $genomes{$st_taxid}{release} > $genomes{$refstrains{$sp_taxid}}{release}) {
				$refstrains{$sp_taxid} = $st_taxid;
			} else {
				next if ($genomes{$st_taxid}{assembly} > $genomes{$refstrains{$sp_taxid}}{assembly} || $genomes{$st_taxid}{release} < $genomes{$refstrains{$sp_taxid}}{release});
				warn "INFO: Multiple reference genomes for same strain taxid (without reads assigned) and date/completeness found.\n",
				"\tReference genome with accession '$genomes{$refstrains{$sp_taxid}}{accession}' is kept while '$genomes{$refstrains{$st_taxid}}{accession}' ",
				"at line '$.' is dropped\n" if ($verbose > 1);
				next;
			}
			
		} else {
			if ($genomes{$st_taxid}) {
				$refstrains{$sp_taxid} = $st_taxid;
			}
		}
	}
}
undef @missing_refs;

#check if all species genomes linked to their respective reference genomes
foreach my $taxid (keys %refstrains) {
	if (! exists $genomes{$taxid}) {
		$genomes{$taxid} = $genomes{$refstrains{$taxid}};
	}
}

##########################################################################
#MODIFY READ COUNTS
#(reassign mode) add reads of strains with no refgenome to strain level with equal distribution to other strains with refgenome
#distribute all reads on species level equally to strains with reads assigned in the correct ratio of their reads assigned

#go through all species again and check corresponding refseq ref genome or warn if not present
#count all assigned reads at respective level for later calculations
my $total_reads;
#total unassigned reads
my $total_ua_reads;
#total reassigned reads, strain without refgenome added to species level
my $total_st2sp_reads;
my $total_st2sp_reassignments;
#total reassigned reads, species to strain(s)
my $total_sp2st_reads;
my $total_sp2st_reassignments;

#create hash for all taxids with a ref_genome assigned for later separation of FASTQ-files in extract_bacreads_kreport.pl
my %taxid_refgenome;

#Debug variable
my $sequence;

#reassign reads from strains without associated genomes to parent species/strain
foreach my $taxid (sort {$a <=> $b} @toplvl_species) {
	#print "\nsp: $taxid\t$species{$taxid}{root_ass}\t$species{$taxid}{name}\n";
	reassign_strain_reads($species{$taxid}{strains});
}

#subroutine for reassignment of all strains without a genome
sub reassign_strain_reads
{
	my $hstrains = shift;
	if ($hstrains) {
		for my $st_taxid (keys %{$hstrains}) {
			
			#start from lowest branch to push up all reads of strains without a refgenome
			reassign_strain_reads($hstrains->{$st_taxid}{strains});
			
			my $h_cur = $hstrains->{$st_taxid};
			
			#check if strain has reads assigned
			if ($h_cur->{root_ass}) {
				
				#STRAIN w/o refgenome
				#if mode reassign, reassign all reads of a strain without a reference to parent strain/species
				if (! ($genomes{$st_taxid} || $refstrains{$st_taxid}) ) {
					
					if ( ! $no_strain_ra) {
						
						#add reads to parent ID
						my $h_parent = $species{$st_taxid}{strainof};
						$h_parent->{root_ass} += $h_cur->{root_ass};
						
						warn "INFO:\t(st2sp) Reassigned '$h_cur->{root_ass}' reads from strain with no refgenome '$h_cur->{taxid}'/'$h_cur->{name}'",
						" to parent species '$h_parent->{taxid}'/'$h_parent->{name}' -> '$h_parent->{root_ass}' reads assigned to parent\n" if $verbose > 1;
						
						#remember taxid to filter out all classified sequences since they will be replaced by simulated sequences
						$taxid_refgenome{st2sp}{$h_cur->{taxid}} = $h_cur->{root_ass};
						
						#set assigned reads of strain to 0 as they are now assigned to parent species
						$total_st2sp_reads += $h_cur->{root_ass};
						$total_st2sp_reassignments++;
						$h_cur->{root_ass} = 0;
						$h_cur->{root_read} = 0;
						
					}
				} else {
					#remember taxid to filter out all classified sequences since they will be replaced by simulated sequences
					$taxid_refgenome{wref}{$st_taxid} = undef;
					print "wref1: $st_taxid\t-\n";
					
				}
#				$sequence++;
#				print "($sequence) st: $st_taxid\t$hstrains->{$st_taxid}{root_ass}\t$hstrains->{$st_taxid}{name}\n";
				
			} else {
#				$sequence++;
#				print "(X$sequence) st: $st_taxid\t$hstrains->{$st_taxid}{root_ass}\t$hstrains->{$st_taxid}{name}\n";
			}
			
		}
	}
}


#########################################################


#subroutine for reassignment of all strains without a genome
sub reassign_species_reads
{
	my $hspecies = shift;
	my $sp_taxid = $hspecies->{taxid};
	
	#check if reads are assigned to current level
	if ($hspecies->{root_ass}) {
		
		#check if there are any strains present where reads could be reassigned
		if ($hspecies->{strains}) {
			
			#sum up all reads assigned to the strains of current species
			my $cur_strainreads;
			my $nr_strains;
			
			for my $st_taxid (keys %{$hspecies->{strains}}) {
				next unless ($hspecies->{strains}{$st_taxid}{root_ass} > 0);				
				$cur_strainreads += $hspecies->{strains}{$st_taxid}{root_ass};
				$nr_strains++;
			}
			
			#distribute all reads at species level to all strains with reads assigned in the same ratio as they have reads assigned 
			#in regard to all reads rooted as species lvl
			if ($cur_strainreads) {
				my $roundcheck;
				foreach my $st_taxid (keys %{$hspecies->{strains}}) {
					
					#previously, all strains with no reference genome had their reads assigned to their parent,
					#so all strains with no reads assigned do not have a reference genome or represent a reference genome for the species 
					#which is not used as other strains with genomes present have reads assigned
					next unless ($hspecies->{strains}{$st_taxid}{root_ass} > 0);
					
					#add reads assigned to species level proportional to relative abundance of strains according to classification
					#in case of rounding errors, correct in next step
					#<reads assigned2species> += (<reads assigned2strain>/<all reads assigned to strains of current species>)*<reads_assigned2sspecies>
					my $reads2add = sprintf("%.0f", (($hspecies->{strains}{$st_taxid}{root_ass} / $cur_strainreads) * $hspecies->{root_ass}));
					$roundcheck += $reads2add;
					$hspecies->{strains}{$st_taxid}{root_ass} += $reads2add;
					
					warn "INFO:\t(sp2st) Reassigned '$reads2add' reads from species '$hspecies->{taxid}'/'$hspecies->{name}' to strain ",
						 "'$st_taxid'/'$hspecies->{strains}{$st_taxid}{name}' \t#from total parent '$hspecies->{root_ass} reads'\n" if ($reads2add && $verbose > 1);
				}
				
				my $rounding_diff = $roundcheck - $hspecies->{root_ass};
				#check if reads were gained/lost due to rounding
				if ($rounding_diff) {
					#correct rounding error by adding/removing a read starting from strain with highest number of reads assigned
					warn "INFO:\t(Correction) Proportional read assignment off by '" . ($roundcheck - $hspecies->{root_ass}) . 
						 "' read(s) for species: '$hspecies->{taxid}'/'$hspecies->{name} as $roundcheck/$hspecies->{root_ass} reads ",
						 "were distributed to strains from species level\n" if $verbose > 1;
					#sort strain taxid descending by assigned reads
					my @sorted_keys = sort { $hspecies->{strains}{$b}{root_ass} <=> $hspecies->{strains}{$a}{root_ass} } keys %{$hspecies->{strains}};
					
					foreach my $st_taxid (@sorted_keys) {
						if ($rounding_diff > 0) {
							$hspecies->{strains}{$st_taxid}{root_ass}--;
							$rounding_diff--;
							$roundcheck--;
							warn "INFO:\t(Correction) Correct by removing 1 read from strain '$st_taxid'/'$hspecies->{strains}{$st_taxid}{name}'\n" if $verbose > 1;
						} elsif ($rounding_diff < 0) {
							$hspecies->{strains}{$st_taxid}{root_ass}++;
							$rounding_diff++;
							$roundcheck++;
							warn "INFO:\t(Correction) Correct by adding 1 read to strain '$st_taxid'/'$hspecies->{strains}{$st_taxid}{name}'\n" if $verbose > 1;
						} else {
							last;
						}
					}
				}
				
				$taxid_refgenome{wref}{$sp_taxid} = 0;
				print "wref2: $sp_taxid\t0\n";
				
				#check if any strain was reassigned to species from which successfully sequences could be assigned to other strains with refgenomes
				foreach my $st_taxid ( keys %{$hspecies->{strains}}) {
					if (exists $taxid_refgenome{st2sp}{$st_taxid}) {
						$taxid_refgenome{wref}{$st_taxid} = 0;
						print "wref3: $st_taxid\t0\n";
					}
				}
				
				#save number of reassigned reads from species to strains
				$total_sp2st_reads += $roundcheck;
				$total_sp2st_reassignments++;
				
				$hspecies->{root_ass} = 0;
			
			#reads assigned to species, strain(s) present but no reads assigned to strain(s) -> assign to reference strain of species
			} elsif (! $cur_strainreads && $genomes{$sp_taxid}) {
				
				#remember taxid to filter out all classified sequences since they will be replaced by simulated sequences
				print "wref4s set $sp_taxid\tfrom $taxid_refgenome{wref}{$sp_taxid} add $hspecies->{root_ass}\n";
				$taxid_refgenome{wref}{$sp_taxid} = $hspecies->{root_ass} unless $taxid_refgenome{wref}{$sp_taxid};
				print "wref4s set $sp_taxid\t$hspecies->{root_ass}\n";
				
				#check if any strain was reassigned to species from which successfully sequences could be assigned to other strains with refgenomes
				foreach my $st_taxid ( keys %{$hspecies->{strains}}) {
					if (exists $taxid_refgenome{st2sp}{$st_taxid}) {
						$taxid_refgenome{wref}{$st_taxid} = 0;
						print "wref5: $st_taxid\t0\n";
					}
				}
	
			} else {
				
				#no refstrain for species present and no reads assigned to species -> unassigned species counts
				die "BUG: Species '$sp_taxid' should have a reference genome assigned, but has not despite a genome $genomes{$sp_taxid}{species_tax}/$genomes{$sp_taxid}{organism} being present for the species\n" 
					if ($genomes{$sp_taxid} || $refstrains{$sp_taxid});
				
				warn "INFO:\t(Unassigned): No reference genomes could be assigned for species '$sp_taxid' with $hspecies->{root_ass} reads assigned\n" if $verbose > 1;
				foreach my $st_taxid ( keys %{$hspecies->{strains}}) {
					if (exists $taxid_refgenome{st2sp}{$st_taxid}) {
						warn "INFO: (Unassigned): No reference genomes could be assigned for strain '$st_taxid'\n" if $verbose > 1;
					}
				}
			}
		
		#current species does not have a strain assigned	
		} else {
			if ($genomes{$sp_taxid} || $refstrains{$sp_taxid}) {
				$taxid_refgenome{wref}{$sp_taxid} = 0;
				print "wref6: $sp_taxid\t0\n";
			} else {
				warn "INFO: (Unassigned): No reference genomes could be assigned for species '$sp_taxid' with $hspecies->{root_ass} reads assigned\n" if $verbose > 1;
			}
		}
		
	}
	
	#reassign reads for subsequent strains if they have strains themself
	if ($hspecies->{strains}) {
		for my $st_taxid (keys %{$hspecies->{strains}}) {
			reassign_species_reads($hspecies->{strains}{$st_taxid});
		}
	}
}
print ">>>> before wref7\n";
foreach my $taxid (sort {$a <=> $b} keys %{$taxid_refgenome{wref}}) {
	print ">>>> $taxid\t$taxid_refgenome{wref}{$taxid}\n";
}

#reassign reads from strains without associated genomes to parent species/strain
foreach my $taxid (sort {$a <=> $b} @toplvl_species) {
	#print "\nsp: $taxid\t$species{$taxid}{root_ass}\t$species{$taxid}{name}\n";
	reassign_species_reads($species{$taxid});
}

if ($total_st2sp_reads) {
	warn "INFO: For '$total_st2sp_reassignments' strains with no refgenome reads were reassigned to species before distribution of species level reads ",
		 "to respective strains with reference genome\n";
	warn "\t#reads: $total_st2sp_reads/$species{2}{root_read} or ", sprintf("%.2f",(($total_st2sp_reads/$species{2}{root_read})*100)),"\% of all bacterial reads\n";
}
if ($total_sp2st_reads) {
	warn "INFO: For '$total_sp2st_reassignments' species, reads were reassigned to the respective strains\n";
	warn "\t#reads: $total_sp2st_reads/$species{2}{root_read} or ", sprintf("%.2f",(($total_sp2st_reads/$species{2}{root_read})*100)),"\%of all bacterial reads\n";
}

##########################################################################
#read in all currently present reference genome paths in provided directory
my %refgenomes;

my @reffiles = glob "'${refseq_folder}/*_genomic.fna*'";
foreach (@reffiles) {
	my $base = basename("$_");
	$refgenomes{$base} = 1;
}

##########################################################################
#write ART input files and all unassigned entries to file
my $GI = w_file("$outdir/genomeInfo.txt");
my $AF = w_file("$outdir/abundanceFile.txt");
my $IF = w_file("$outdir/fullprofile.tsv");

my $genome_count;
print $IF "Abundance\tNCBI TaxID\tName\tReference filename\tReference genome length\n";

my $UA = w_file("$outdir/norefgenome.tsv");
print $UA "#No reference genome found for the following species (and their strains)\nTaxid\tReads assigned\tName\n";

sub check_refgenomes {
	my $hspecies = shift;
	my $sp_taxid = $hspecies->{taxid};
	
	#check if reads are assigned to current level
	if ($hspecies->{root_ass}) {
		
		#check if reference genome is available or if reads were assigned to reference strain
		if (! $genomes{$sp_taxid}) {
			
			print $UA "$sp_taxid\t$hspecies->{root_ass}\t$hspecies->{name}\n";
			$total_ua_reads += $hspecies->{root_ass};
			
		} else {
			
			#remember taxid to filter out all classified sequences since they will be replaced by simulated sequences
			$taxid_refgenome{wref}{$sp_taxid} += $hspecies->{root_ass};
			print "wref7a add $hspecies->{root_ass} with result $sp_taxid\t$hspecies->{root_ass}\n" unless $taxid_refgenome{wref}{$sp_taxid};
			
			$total_reads += $hspecies->{root_ass};
			get_refgenome($sp_taxid);
			
			#write to ART input files
			my $base = basename("$genomes{$sp_taxid}{ftp}");
			$genome_count++;
			
			print $AF "${base}_genomic.fna\t$hspecies->{root_ass}\n";
			print $GI "${base}_genomic.fna\t$refgenomes{$sp_taxid}{genomelength}\t1\n";
			print $IF "$hspecies->{root_ass}\t$sp_taxid\t$genomes{$sp_taxid}{organism}\t${base}_genomic.fna.gz\t$refgenomes{$sp_taxid}{genomelength}\n";
			
		}
	}
	
	#reassign reads for subsequent strains if they have strains themself
	if ($hspecies->{strains}) {
		for my $st_taxid (keys %{$hspecies->{strains}}) {
			check_refgenomes($hspecies->{strains}{$st_taxid});
		}
	}
}
#load info for or download reference genomes for all entries with reads assigned
#count all successfully or unsucessfully assigned read counts
foreach my $taxid (sort {$a <=> $b} @toplvl_species) {
	check_refgenomes($species{$taxid});
}
close $UA;
close $AF;
close $GI;
close $IF;

print "Selected $genome_count genomes\n";

#write out all counts for all taxid's with assigned/reassigned refgenomes which will be replaced
my $TWR = w_file("$outdir/taxa_w_refgenome.tsv");
print $TWR "taxid\treads_after_reassignment\n";
foreach my $taxid (sort {$a <=> $b} keys %{$taxid_refgenome{wref}}) {
	print $TWR "$taxid\t$taxid_refgenome{wref}{$taxid}\n";
	print ">>>> $taxid\t$taxid_refgenome{wref}{$taxid}\n";
}
close $TWR;

print ">>>>>>>>>>";
print "Assigned bacterial reads at species level or lower with refgenome: $total_reads/$species{2}{root_read} or ", 
	sprintf("%.2f",(($total_reads/$species{2}{root_read})*100)),"\%\n";
	
print "Assigned bacterial reads at species level or lower with no refgenome: $total_ua_reads/$species{2}{root_read} or ", 
	sprintf("%.2f",(($total_ua_reads/$species{2}{root_read})*100)),"\%\n" if ($total_ua_reads);
	
print "Assigned bacterial reads higher than species level: $total_reads_nonspecieslvl/$species{2}{root_read} or ", 
	sprintf("%.2f",(($total_reads_nonspecieslvl/$species{2}{root_read})*100)),"\%\n" if ($total_reads_nonspecieslvl);
	
print "Number of bacterial reads of total reads: $species{2}{root_read}/",($species{0}{root_read} + $species{1}{root_read})," or ", 
	sprintf("%.2f",(($species{2}{root_read}/($species{0}{root_read} + $species{1}{root_read}))*100)),"\%\n";
	
if ((defined $total_ua_reads && defined $total_reads_nonspecieslvl ) && ($total_reads + $total_ua_reads + $total_reads_nonspecieslvl) != $species{2}{root_read}) {
	my $total_lost =  $species{2}{root_read} - $total_reads - $total_ua_reads - $total_reads_nonspecieslvl;
	print "Remaining $total_lost/$species{2}{root_read} or ", 
		sprintf("%.2f",(($total_lost/$species{2}{root_read})*100)),"\% are multimapping reads lost due to rounding\n";
}
print ">>>>>>>>>>";

print "$0 took ",runtime(tv_interval($t0)), " to run\n";


##########################################################################
#subroutines
##########################################################################
#check if refgenome is present, else download
sub get_refgenome {
	my ($taxid) = @_;
	my $base = basename("$genomes{$taxid}{ftp}");
	my $filename = "${base}_genomic.fna.gz";
	
	#download genome only if not present already
	if (! $refgenomes{$filename}) {
		fetchfile("$genomes{$taxid}{ftp}/${filename}","$refseq_folder");
		print "Downloaded '$filename'\n";		
	} 
	
	#save length of refgenome
	my $total_length;
	my $IN = r_file("$refseq_folder/$filename");
	local $/=">";
	while (<$IN>){
		chomp;
        next unless /\w/;
        s/>$//gs;
        my @chunk = split /\n/;
        my $header = shift @chunk;
        $total_length += length join "", @chunk;
	}
	close $IN;
	local $/="\n";
	$refgenomes{$taxid}{genomelength} = $total_length;
}
##########################################################################
sub assign_species {
	my ($taxid,$r_read,$r_ass,$level,$linenr,$name) = @_;
	$species{$taxid} = { 
		root_read => $r_read,
		root_ass => $r_ass,
		level => $level,
		line => $linenr,
		name => $name,
		taxid => $taxid
	};
}
##########################################################################
sub assign_strain {
	my ($prev_href,$taxid,$r_read,$r_ass,$level,$linenr,$name) = @_;
	die "$prev_href|$taxid|$r_read|$r_ass|$level|$linenr|$name" unless (ref($prev_href) eq "HASH");
	$prev_href->{strains}{$taxid} = { 
		root_read => $r_read,
		root_ass => $r_ass,
		level => $level,
		line => $linenr,
		name => $name,
		taxid => $taxid
	};
	$species{$taxid}{strainof} = $prev_href;
}
##########################################################################
sub assign_genome {
	my ($taxid,$acc,$s_taxid,$name,$ass_lvl,$ftp,$rel_date,$refcat) = @_;
	$genomes{$taxid} = {
		accession => $acc,
		species_tax => $s_taxid,
		organism => $name,
		assembly => $ass_lvl,
		ftp => $ftp,
		release => $rel_date,
		category => $refcat
	};
}
##########################################################################
sub fetchfile
{
	my ($uri,$targetdir) = @_;
	if (! $targetdir) {
		$targetdir = ".";
	}
	my $ff = File::Fetch->new(uri => "$uri");
	my $file = $ff->fetch( to => "$targetdir") or die $ff->error();
	return $file;
}
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
				$rtime = int($time / 86400) . "d " . int(($time % 86400) / 3600) . "h " . int((($time % 86400) % 3600) / 60) . "m " 
						 . (((time % 86400) % 3600) % 60) . "s";
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
	my $err = shift || 0;
	print STDERR <<EOD;

Usage: $0 [Parameters]
	
	Generate abundance profile with according references from RefSeq for classified NGS reads
	
	==== Parameters ====
	
	-k/--kraken-report		Kraken(style) report 
	-a/--assembly-summary		NCBI Assembly summary table
	
	-o/--outdir 			Output directory
	-R/--refgenomes			Directory to safe reference genomes in gzipped fasta format
	
	--no-reassign			No reassignment of read counts from strains without a reference to other 
	               			classified strains/reference genomes of same species (off)
	
	
	-v/--verbose			Print detailed information for each step, can be used multiple times
	
	
	-h/--help			Prints this helpmessage
	
EOD
	exit $err;
}
