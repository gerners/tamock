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
my $abs_species_threshold = 1;

GetOptions(	"kraken-report=s" => \$kraken_report,
			"assembly-summary=s" => \$assembly_summary,
			"refgenomes|R=s" => \$refseq_folder,
			"outdir=s" => \$outdir,
			"min-abund-species:i" => \$abs_species_threshold,
			"no-reassign" => \$no_strain_ra,
			"verbose+" => \$verbose,
			"debug+" => \$debug,
			"help" => \&print_help);

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
my $prev_level = 0;
my $prev_taxid;
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
	
	#save nr of unclassified reads
	if ($. == 1 && $rank eq "U") {
		assign_species($taxid,$r_read,$r_ass,$rank,$.,$name);
		print $DB1 "$r_read\t$r_ass\t$rank\t$level\t$taxid\t$.\t$name\n" if $debug;
	#save nr of all classified reads on root level
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
	}
	next unless ($bacflag == 1);
	
	#check if current line is lower level (= more leading whitespaces) than the last matched species level 
	#and therefore a strain, 
	#if not remove the value to skip all further lines until another species is matched
	if ($prev_level > 0 && $level > $prev_level) {
		die "Unexpected line below species line found at line nr '$.', level '$level', prevlevel '$prev_level'\n" unless ($rank =~ /[-S]/);
		
		#skip groups of strains with no reads assigned, as they are all assigned to the strains under 
		#the subspecies group in this case
		next unless $r_ass;
		
		assign_strain($prev_taxid,$taxid,$r_read,$r_ass,$rank,$.,$name);
		print $DB1 "$r_read\t$r_ass\t$rank\t$level\t$taxid\t$.\t$name\n" if $debug;
		$species{$taxid}{strainof} = $prev_taxid;
		next;
		
		
	} else {
		$prev_level = 0;
		$prev_taxid = undef;
	}
	
	#check if minimal number of reads is assigned at species level and save it if
	if ($rank eq "S" && $r_read >= $abs_species_threshold) {
		$prev_level = $level;
		$prev_taxid = $taxid;
		
		#save species with taxid as key
		assign_species($taxid,$r_read,$r_ass,$rank,$.,$name);
		print $DB1 "$r_read\t$r_ass\t$rank\t$level\t$taxid\t$.\t$name\n" if $debug;
	} else {
		$total_reads_nonspecieslvl += $r_ass if ($r_ass);
	}
}
close $KR;

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
	#skip strains which have no own strain ID but use the taxid of the species level while there is already a reference present
	next if ($s_taxid == $taxid && $refstrains{$s_taxid});
	
	#check if another reference is already present for current ref genome, if, test for better completeness 
	#via assembly level or most resent release date
	#first check, if genome is representative/reference genome i.e. not with "na", set as species lvl ref genome
	
	#STRAIN present in krakenstyle report
	if ($species{$taxid} ) {
		
		#STRAIN & REFERENCE/REPRESENTATIVE GENOME
		#always save reference/representative genomes
		if ($refcat ne "na") { 
			
			#reference > representative
			next if ($genomes{$taxid} && $genomes{$taxid}{category} eq "reference genome");
			
			assign_genome($taxid,$acc,$s_taxid,$name,$ass_lvl,$ftp,$rel_date,$refcat);
			$refstrains{$s_taxid} = $taxid;

		#STRAIN & NONREF
		} else {
			
			#if multiple reference genomes with the same taxid are present, rank in the order of
			#reference genome > assembly status > date
			
			#strain has already a reference assigned -> check if current refseq genome should replace present entry
			if ($genomes{$taxid}) {
				next if ($genomes{$taxid}{category} ne "na");		#refgenome for current taxid is reference genome -> do not replace
				if ($ass_lvl < $genomes{$taxid}{assembly}) {		#assembly level 1 => complete; 4 => contig
					assign_genome($taxid,$acc,$s_taxid,$name,$ass_lvl,$ftp,$rel_date,$refcat);
				} elsif ($ass_lvl == $genomes{$taxid}{assembly} && $rel_date > $genomes{$taxid}{release}) {
					assign_genome($taxid,$acc,$s_taxid,$name,$ass_lvl,$ftp,$rel_date,$refcat);
				} else {
					next if ($ass_lvl > $genomes{$taxid}{assembly} || $rel_date < $genomes{$taxid}{release});
					warn "INFO: Multiple reference genomes for same strain taxid and date/completeness found.\n\tReference genome with accession '$genomes{$taxid}{accession}' is kept while '$acc' at line '$.' is dropped\n" if ($verbose > 1);
					next;
				}
			#strain has no reference assigned -> assign current refseq genome as reference
			} else {
				assign_genome($taxid,$acc,$s_taxid,$name,$ass_lvl,$ftp,$rel_date,$refcat);
			}		
		}
	
	#STRAIN NOT PRESENT
	#SPECIES & REFERENCE STRAIN
	#save reference strain genomes even if strain is not present in data in case of reads only classified on species level 
	} elsif ($species{$s_taxid} && $refcat ne "na") {
		
		#reference genome > representative genome
		next if ($refstrains{$s_taxid} && $genomes{$refstrains{$s_taxid}}{category} eq "reference genome");
		
		assign_genome($taxid,$acc,$s_taxid,$name,$ass_lvl,$ftp,$rel_date,$refcat);
		
		#assign detected reference strain for species
		assign_strain($s_taxid,$taxid,0,0,"-",$species{$s_taxid}{line},$name);
		$species{$taxid}{strainof} = $s_taxid;
		$refstrains{$s_taxid} = $taxid;
	
	#STRAIN NOT PRESENT
	#OTHER STRAIN for species without any strain and no reference strain yet to have at least one strain as a reference genome available
	} else {
		if ($species{$s_taxid}) {
						
			next if ($refstrains{$s_taxid} && $genomes{$refstrains{$s_taxid}}{category} ne "na");
			
			#if there is already a reference/representative strain present for the species, skip
			if ($refstrains{$s_taxid}) {
				
				#some values might be unitialized
				 no warnings 'uninitialized';
				
				if ($ass_lvl < $genomes{$refstrains{$s_taxid}}{assembly}) {		#assembly level 1 => complete; 4 => contig
					assign_genome($taxid,$acc,$s_taxid,$name,$ass_lvl,$ftp,$rel_date,$refcat);
					assign_strain($s_taxid,$taxid,0,0,"-",$species{$s_taxid}{line},$name);
					$species{$taxid}{strainof} = $s_taxid;
					$refstrains{$s_taxid} = $taxid;
					
				} elsif ($ass_lvl == $genomes{$refstrains{$s_taxid}}{assembly} && $rel_date > $genomes{$taxid}{release}) {
					assign_genome($taxid,$acc,$s_taxid,$name,$ass_lvl,$ftp,$rel_date,$refcat);
					assign_strain($s_taxid,$taxid,0,0,"-",$species{$s_taxid}{line},$name);
					$species{$taxid}{strainof} = $s_taxid;
					$refstrains{$s_taxid} = $taxid;
					
				} else {
					next;
				}
				
			} else {
				assign_genome($taxid,$acc,$s_taxid,$name,$ass_lvl,$ftp,$rel_date,$refcat);
				assign_strain($s_taxid,$taxid,0,0,"-",$species{$s_taxid}{line},$name);
				$species{$taxid}{strainof} = $s_taxid;
				$refstrains{$s_taxid} = $taxid;
			}
		}
	}
}
close $NCBI;

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
#STRAINS
foreach my $taxid (sort {$a <=> $b} keys %species) {	

	#skip entry for unclassified (0) and root (1)
	next if $taxid < 3;

	#check if strain
	if ($species{$taxid}{strainof}) {
		
		#STRAIN w/o refgenome
		#if mode reassing, reassign all reads of a strain without a reference to species level from which it will be
		#added equally to all remaining strains
		if (! $genomes{$taxid} && $species{$species{$taxid}{strainof}}{strains}{$taxid}{root_ass}) {
			
			unless ($no_strain_ra) {
				#add reads to parent species ID
				$species{$species{$taxid}{strainof}}{root_ass} += $species{$species{$taxid}{strainof}}{strains}{$taxid}{root_ass};
				
				warn "INFO:\t(st2sp) Reassigned '$species{$species{$taxid}{strainof}}{strains}{$taxid}{root_ass}' reads from strain with no refgenome to '$taxid'/'$species{$species{$taxid}{strainof}}{strains}{$taxid}{name}' to species lvl '$species{$taxid}{strainof}'/'$species{$species{$taxid}{strainof}}{name}'\t#reads assigned now to species '$species{$species{$taxid}{strainof}}{root_ass}'\n" if $verbose;
				
				#set assigned reads of strain to 0 as they are now assigned to parent species
				$total_st2sp_reads += $species{$species{$taxid}{strainof}}{strains}{$taxid}{root_ass};
				$total_st2sp_reassignments++;
				$species{$species{$taxid}{strainof}}{strains}{$taxid}{root_ass} = 0;
			}
		} elsif ($genomes{$taxid}) {
			#set flag for genome found at current strain
			$species{$species{$taxid}{strainof}}{strains}{$taxid}{genome}++;
		}
	}
}

#SPECIES
foreach my $taxid (sort {$a <=> $b} keys %species) {
	
	#skip entry for unclassified (0) and root (1)
	next if $taxid < 3;		
	
	if ($species{$taxid}{root_ass} && $species{$taxid}{strains}) {
			
		#sum up all reads assigned to the strains of current species
		my $cur_strainreads;
		my $nr_strains;
		
		#check how many strains have actual reads assigned i.e. on which strains species reads will be assigned to
		#if no reads are assigned to any strains of current species, reads will be assigned to a single reference strain
		foreach my $strain ( keys %{$species{$taxid}{strains}}) {
			next unless ($species{$taxid}{strains}{$strain}{root_ass} > 0);				
			$cur_strainreads += $species{$taxid}{strains}{$strain}{root_ass};
			$nr_strains++;
		}
		
		#distribute all reads at species level to all strains with reads assigned in the same ratio as they have reads assigned in regard to all reads rooted as species lvl
		if ($cur_strainreads) {
			my $roundcheck;
			foreach my $strain (keys %{$species{$taxid}{strains}}) {
				
				#previously, all strains with no reference genome had their reads assigned to species level,
				#so all strains with no reads assigned do not have a reference genome or represent a reference genome for the species which is not needed as other strains with reference genome have reads assigned
				next unless ($species{$taxid}{strains}{$strain}{root_ass} > 0);
				
				#add reads assigned to species level proportional to relative abundance of strains according to classification
				#in case of rounding errors, correct in next step
				#<reads assigned2species> += (<reads assigned2strain>/<all reads assigned to strains of current species>)*<reads_assigned2sspecies>
				my $reads2add = sprintf("%.0f", (($species{$taxid}{strains}{$strain}{root_ass} / $cur_strainreads) * $species{$taxid}{root_ass}));
				$roundcheck += $reads2add;
				$species{$taxid}{strains}{$strain}{root_ass} += $reads2add;
				
				warn "INFO:\t(sp2st) Reassigned '$reads2add' reads from species '$taxid'/'$species{$taxid}{name}' to strain '$strain'/'$species{$taxid}{strains}{$strain}{name}' \t#reads species lvl '$species{$taxid}{root_ass}'\n" if ($reads2add && $verbose);
			
			#if reads are only present at species level but none are assigned to any strain, use a single reference strain for said species
			} 
			
			my $rounding_diff = $roundcheck - $species{$taxid}{root_ass};
			#check if reads were gained/lost due to rounding
			if ($rounding_diff) {
				#correct rounding error by adding/removing a read starting from strain with highest number of reads assigned
				warn "INFO:\t(Correction) Proportional read assignment off by '" . ($roundcheck - $species{$taxid}{root_ass}) . "' read(s) for species: '$taxid'/'$species{$taxid}{name} as $roundcheck/$species{$taxid}{root_ass} reads were distributed to strains from species level\n" if $verbose;
				#sort strain taxid descending by assigned reads
				my @sorted_keys = sort { $species{$taxid}{strains}{$b}{root_ass} <=> $species{$taxid}{strains}{$a}{root_ass} } keys %{$species{$taxid}{strains}};
				
				foreach my $strain (@sorted_keys) {
					if ($rounding_diff > 0) {
						$species{$taxid}{strains}{$strain}{root_ass}--;
						$rounding_diff--;
						$roundcheck--;
						warn "INFO:\t(Correction) Correct by removing 1 read from strain '$strain'/'$species{$taxid}{strains}{$strain}{name}'\n" if $verbose;
					} elsif ($rounding_diff < 0) {
						$species{$taxid}{strains}{$strain}{root_ass}++;
						$rounding_diff++;
						$roundcheck++;
						warn "INFO:\t(Correction) Correct by adding 1 read to strain '$strain'/'$species{$taxid}{strains}{$strain}{name}'\n" if $verbose;
					} else {
						last;
					}
				}
			}
			
			#save number of reassigned reads from species to strains
			$total_sp2st_reads += $roundcheck;
			$total_sp2st_reassignments++;
			
			$species{$taxid}{root_ass} = 0;
		
		#reads assigned to species, strain present but no reads assigned to strain
		} else {
			
			#check if a reference strain is found for a species without any reads classified to strains but reads assigned to species level
			if ($refstrains{$taxid}) {
				$species{$taxid}{strains}{$refstrains{$taxid}}{root_ass} += $species{$taxid}{root_ass};
				$species{$taxid}{root_ass} = 0;
			} else {
				#if no strains are found and no reference strain, check if a genome is assigned at species lvl
				#next if a reference genome on species lvl is found or no strains were found for the current species
				#warn "$species{$taxid}{root_ass}\n" unless $genomes{$taxid};
				next if ($genomes{$taxid} || ! $nr_strains);
				
				#safety check, else only one genome with the most complete and resent assembly should have been selected
				die "More than one strain ('$nr_strains'/'$cur_strainreads') with no reads assigned and which is not a reference strain found for species '$taxid' with '$species{$taxid}{root_ass}' reads\n" if ($nr_strains != 1);
				foreach my $strain ( keys %{$species{$taxid}{strains}}) {
					$species{$taxid}{strains}{$strain}{root_ass} += $species{$taxid}{root_ass};
					$species{$taxid}{root_ass} = 0;
				}
			}
		}
	}
}

if ($total_st2sp_reads) {
	warn "INFO: For '$total_st2sp_reassignments' strains with no refgenome reads were reassigned to species before distribution of species level reads to respective strains with reference genome\n";
	warn "\t#reads: $total_st2sp_reads/$species{2}{root_read} or ", sprintf("%.2f",(($total_st2sp_reads/$species{2}{root_read})*100)),"\% of all bacterial reads\n";
}
if ($total_sp2st_reads) {
	warn "INFO: For '$total_sp2st_reassignments' species, reads were reassigned to the respective strains\n";
	warn "\t#reads: $total_sp2st_reads/$species{2}{root_read} or ", sprintf("%.2f",(($total_sp2st_reads/$species{2}{root_read})*100)),"\%of all bacterial reads\n";
}


##########################################################################
#read in all currently present reference genomes in provided directory
my %refgenomes;

my @reffiles = glob "'${refseq_folder}/*_genomic.fna*'";
foreach (@reffiles) {
	my $base = basename("$_");
	$refgenomes{$base} = 1;
}

##########################################################################
#write all unassigned entries to file
my $UA = w_file("$outdir/norefgenome.tsv");
print $UA "#No reference genome found for the following species (and their strains)\nTaxid\tReads assigned\tName\n";

#open reference genome folder
foreach my $taxid (sort {$a <=> $b} keys %species) {
	next if $taxid < 3;
	
	#some saved reference strains with no reassigned reads might have unitialized values
	no warnings 'uninitialized';
	
	#check species first
	if (! $species{$taxid}{strainof}) {
		
		next unless $species{$taxid}{root_ass};
		
		#check if reference genome is available or if reads were assigned to reference strain
		if (! $genomes{$taxid}) {
			#skip if no reads were assigned to this species level
			next if ($species{$taxid}{root_ass} == 0);
			print $UA "$taxid\t$species{$taxid}{root_ass}\t$species{$taxid}{name}\n"; 
			$total_ua_reads += $species{$taxid}{root_ass} if ($species{$taxid}{root_ass} > 0);
		} else {
			$total_reads += $species{$taxid}{root_ass} if ($species{$taxid}{root_ass} > 0);
			check_refgenome($taxid);
		}
		
	#all remaining entries are strains
	} else {
		next if ($species{$species{$taxid}{strainof}}{strains}{$taxid}{root_ass} == 0);
		if (! $genomes{$taxid}) {
			die "BUG: Reads of strains with no reference genome not assigned to species level for taxid '$taxid'\n" unless $no_strain_ra;
			print $UA "Strain\t$taxid\t$species{$species{$taxid}{strainof}}{strains}{$taxid}{root_ass}\t$species{$species{$taxid}{strainof}}{strains}{$taxid}{name}\n"; 
			$total_ua_reads += $species{$species{$taxid}{strainof}}{strains}{$taxid}{root_ass};
		} else {
			$total_reads += $species{$species{$taxid}{strainof}}{strains}{$taxid}{root_ass};
			check_refgenome($taxid);
		}
	}
}
close $UA;

print "Assigned bacterial reads at species level or lower with refgenome: $total_reads/$species{2}{root_read} or ", 
	sprintf("%.2f",(($total_reads/$species{2}{root_read})*100)),"\%\n";
	
print "Assigned bacterial reads at species level or lower with no refgenome: $total_ua_reads/$species{2}{root_read} or ", 
	sprintf("%.2f",(($total_ua_reads/$species{2}{root_read})*100)),"\%\n" if ($total_ua_reads);
	
print "Assigned bacterial reads higher than species level: $total_reads_nonspecieslvl/$species{2}{root_read} or ", 
	sprintf("%.2f",(($total_reads_nonspecieslvl/$species{2}{root_read})*100)),"\%\n" if ($total_reads_nonspecieslvl);
	
print "Number of bacterial reads of total reads: $species{2}{root_read}/",($species{0}{root_read} + $species{1}{root_read})," or ", 
	sprintf("%.2f",(($species{2}{root_read}/($species{0}{root_read} + $species{1}{root_read}))*100)),"\%\n";
	
if (($total_reads + $total_ua_reads + $total_reads_nonspecieslvl) != $species{2}{root_read}) {
	my $total_lost =  $species{2}{root_read} - $total_reads - $total_ua_reads - $total_reads_nonspecieslvl;
	print "Remaining $total_lost/$species{2}{root_read} or ", 
		sprintf("%.2f",(($total_lost/$species{2}{root_read})*100)),"\% are multimapping reads lost due to rounding\n";
}



##########################################################################
#copy all respective reference genomes to outputdir, and write ART input files
my $GI = w_file("$outdir/genomeInfo.txt");
my $AF = w_file("$outdir/abundanceFile.txt");
my $IF = w_file("$outdir/fullprofile.tsv");

my $genome_count;
print $IF "Abundance\tNCBI TaxID\tName\tReference filename\tReference genome length\n";
foreach my $taxid (sort {$a <=> $b} keys %genomes) {
	
	#print only reference genomes for species/strains with reads assigned
	if (exists $species{$taxid}{strainof}) {
		next unless $species{$species{$taxid}{strainof}}{strains}{$taxid}{root_ass} > 0;
	} else {
		next unless $species{$taxid}{root_ass} > 0;
	}

	my $base = basename("$genomes{$taxid}{ftp}");
	
	#some saved reference strains/sepcies with no reassigned reads might have unitialized values
	no warnings 'uninitialized';
	my $abundance;
	if ($species{$taxid}{root_ass} > 0) {
		$abundance = $species{$taxid}{root_ass};
		$genome_count++;
	} elsif ($species{$species{$taxid}{strainof}}{strains}{$taxid}{root_ass} > 0) {
		$abundance = $species{$species{$taxid}{strainof}}{strains}{$taxid}{root_ass};
		$genome_count++;
	} else {
		warn "Couldn't find abundance for taxid '$taxid'\n";
	}
	print $AF "${base}_genomic.fna\t$abundance\n";
	print $GI "${base}_genomic.fna\t$refgenomes{$taxid}{genomelength}\t1\n";
	print $IF "$abundance\t$taxid\t$genomes{$taxid}{organism}\t${base}_genomic.fna.gz\t$refgenomes{$taxid}{genomelength}\n";
}
close $AF;
close $GI;
close $IF;

print "Selected $genome_count genomes\n";
print "$0 took ",runtime(tv_interval($t0)), " to run\n";
##########################################################################
#subroutines
##########################################################################
#check if refgenome is present, else download
sub check_refgenome {
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
	my ($taxid,$r_read,$r_ass,$rank,$linenr,$name) = @_;
	$species{$taxid} = { 
		root_read => $r_read,
		root_ass => $r_ass,
		rank => $rank,
		line => $linenr,
		name => $name
	};
}
##########################################################################
sub assign_strain {
	my ($prev_taxid,$taxid,$r_read,$r_ass,$rank,$linenr,$name) = @_;
	$species{$prev_taxid}{strains}{$taxid} = { 
		root_read => $r_read,
		root_ass => $r_ass,
		rank => $rank,
		line => $linenr,
		name => $name,
		genome => 0
	};
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
	my $err = shift || 0;
	print STDERR <<EOD;

Usage: $0 [Parameters]
	
	Generate abundance profile with according references from RefSeq for classified NGS reads
	
	==== Parameters ====
	
	-k/--kraken-report		Kraken(style) report 
	-a/--assembly-summary		NCBI Assembly summary table
	
	-o/--outdir 			Output directory
	-R/--refgenomes			Directory to safe reference genomes in gzipped fasta format
	
	
	--min-abund-species		Minimum number of classified reads for a species and all related strains
	                 		to be included into the mock community (1)
	--no-reassign			No reassignment of read counts from strains without a reference to other 
	               			classified strains/reference genomes of same species (off)
	
	
	-v/--verbose			Print detailed information for each step, can be used multiple times
	
	
	-h/--help			Prints this helpmessage
	
EOD
	exit $err;
}
