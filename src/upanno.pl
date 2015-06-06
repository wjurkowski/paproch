#!/usr/bin/perl -w
use strict;
use warnings;

if ($#ARGV != 0) {die "Program used with parameters [uniprot]\n";}

my @unip=open_file($ARGV[0]);
open(FRES, ">UniprotKBSeqF-aa.txt");
open(PAIR, ">UniprotKBSeqF-pair.txt");
open(RNGS, ">UniprotKBSeqF-range.txt");

my $licz=0;
my (@uniprotAC,$uniID);
my (%binding,%actsite,%metal,%site,%nonstd,%modres,%lipid,%carbohyd);
my (%disulfid,%crosslnk);
my (%transmem,%intramem,%repeat,%cabind,%dnabind,%znfing,%npbind,%motif);
my (%signal,%propep,%transit);
printf FRES "uniprotAC\tuniprotID\tBINDING\tACT_SITE\tMETAL\tSITE\tNON_STD\tMOD_RES\tLIPID\tCARBOHYD\n";
printf PAIR "uniprotAC\tuniprotID\tDISULFID\tCROSSLNK\n";
printf RNGS "uniprotAC\tuniprotID\tTRANSMEM\tINTRAMEM\tREPEAT\tCA_BIND\tDNA_BIND\tZN_FING\tNP_BIND\tMOTIF\tSIGNAL\tPROPEP\tTRANSIT\n";
for(my $i=0;$i<$#unip+1;$i++){	
	my @tab = split(/\s{2,}/,$unip[$i]);
	if($tab[0] eq "AC"){
	  $tab[1]=~s/\s+//;
	  my @new = split(/;/,$tab[1]);
	  push(@uniprotAC,@new);
	}
	elsif($tab[0] eq "ID"){
	  $tab[1]=~s/\s+//;
	  $uniID = $tab[1];
	}
	elsif($tab[0] eq "FT"){
	#Sites	
		if($tab[1] eq "BINDING"){
			if(exists $binding{$uniID}){$binding{$uniID}=$binding{$uniID}.",".$tab[2];}
			else{$binding{$uniID}=$tab[2];}
		}
		if($tab[1] eq "ACT_SITE"){
			if(exists $actsite{$uniID}){$actsite{$uniID}=$actsite{$uniID}.",".$tab[2];}
			else{$actsite{$uniID}=$tab[2];}
		}
		if($tab[1] eq "METAL"){
			if(exists $metal{$uniID}){$metal{$uniID}=$metal{$uniID}.",".$tab[2];}
			else{$metal{$uniID}=$tab[2];}
		}
		if($tab[1] eq "SITE"){
			if(exists $site{$uniID}){$site{$uniID}=$site{$uniID}.",".$tab[2];}
			else{$site{$uniID}=$tab[2];}
		}
	#Amino acids modifications	
		if($tab[1] eq "NON_STD"){
			if(exists $nonstd{$uniID}){$nonstd{$uniID}=$nonstd{$uniID}.",".$tab[2];}
			else{$nonstd{$uniID}=$tab[2];}
		}		
		if($tab[1] eq "MOD_RES"){
			if(exists $modres{$uniID}){$modres{$uniID}=$modres{$uniID}.",".$tab[2];}
			else{$modres{$uniID}=$tab[2];}
		}
		if($tab[1] eq "LIPID"){
			if(exists $lipid{$uniID}){$lipid{$uniID}=$lipid{$uniID}.",".$tab[2];}
			else{$lipid{$uniID}=$tab[2];}
		}
		if($tab[1] eq "CARBOHYD"){
			if(exists $carbohyd{$uniID}){$carbohyd{$uniID}=$carbohyd{$uniID}.",".$tab[2];}
			else{$carbohyd{$uniID}=$tab[2];}
		}
	#Amino acids modifications - covalent bonds	
		if($tab[1] eq "DISULFID"){
			if(exists $disulfid{$uniID}){$disulfid{$uniID}=$disulfid{$uniID}.",".$tab[2]."-".$tab[3];}
			else{$disulfid{$uniID}=$tab[2]."-".$tab[3];}
		}
		if($tab[1] eq "CROSSLNK"){
			if(exists $crosslnk{$uniID}){$crosslnk{$uniID}=$crosslnk{$uniID}.",".$tab[2]."-".$tab[3];}
			else{$crosslnk{$uniID}=$tab[2]."-".$tab[3];}
		}
	#Regions
		if($tab[1] eq "TRANSMEM"){
			if(exists $transmem{$uniID}){$transmem{$uniID}=$transmem{$uniID}.",".$tab[2]."-".$tab[3];}
			else{$transmem{$uniID}=$tab[2]."-".$tab[3];}
		}
		if($tab[1] eq "INTRAMEM"){
			if(exists $intramem{$uniID}){$intramem{$uniID}=$intramem{$uniID}.",".$tab[2]."-".$tab[3];}
			else{$intramem{$uniID}=$tab[2]."-".$tab[3];}
		}
		if($tab[1] eq "REPEAT"){
			if(exists $repeat{$uniID}){$repeat{$uniID}=$repeat{$uniID}.",".$tab[2]."-".$tab[3];}
			else{$repeat{$uniID}=$tab[2]."-".$tab[3];}
		}
		if($tab[1] eq "CA_BIND"){
			if(exists $cabind{$uniID}){$cabind{$uniID}=$cabind{$uniID}.",".$tab[2]."-".$tab[3];}
			else{$cabind{$uniID}=$tab[2]."-".$tab[3];}
		}
		if($tab[1] eq "DNA_BIND"){
			if(exists $dnabind{$uniID}){$dnabind{$uniID}=$dnabind{$uniID}.",".$tab[2]."-".$tab[3];}
			else{$dnabind{$uniID}=$tab[2]."-".$tab[3];}
		}
		if($tab[1] eq "ZN_FING"){
			if(exists $znfing{$uniID}){$znfing{$uniID}=$znfing{$uniID}.",".$tab[2]."-".$tab[3];}
			else{$znfing{$uniID}=$tab[2]."-".$tab[3];}
		}
		if($tab[1] eq "NP_BIND"){
			if(exists $npbind{$uniID}){$npbind{$uniID}=$npbind{$uniID}.",".$tab[2]."-".$tab[3];}
			else{$npbind{$uniID}=$tab[2]."-".$tab[3];}
		}
		if($tab[1] eq "MOTIF"){
			if(exists $motif{$uniID}){$motif{$uniID}=$motif{$uniID}.",".$tab[2]."-".$tab[3];}
			else{$motif{$uniID}=$tab[2]."-".$tab[3];}
		}
	#Molecule processing
		if($tab[1] eq "SIGNAL"){
			if(exists $signal{$uniID}){$signal{$uniID}=$signal{$uniID}.",".$tab[2]."-".$tab[3];}
			else{$signal{$uniID}=$tab[2]."-".$tab[3];}
		}
		if($tab[1] eq "PROPEP"){
			if(exists $propep{$uniID}){$propep{$uniID}=$propep{$uniID}.",".$tab[2]."-".$tab[3];}
			else{$propep{$uniID}=$tab[2]."-".$tab[3];}
		}
		if($tab[1] eq "TRANSIT"){
			if(exists $transit{$uniID}){$transit{$uniID}=$transit{$uniID}.",".$tab[2]."-".$tab[3];}
			else{$transit{$uniID}=$tab[2]."-".$tab[3];}
		}
	}
	elsif($tab[0] eq "//"){
		unless(exists $binding{$uniID}){$binding{$uniID}="NA";}
		unless(exists $actsite{$uniID}){$actsite{$uniID}="NA";}
		unless(exists $metal{$uniID}){$metal{$uniID}="NA";}
		unless(exists $site{$uniID}){$site{$uniID}="NA";}
		unless(exists $nonstd{$uniID}){$nonstd{$uniID}="NA";}
		unless(exists $modres{$uniID}){$modres{$uniID}="NA";}
		unless(exists $lipid{$uniID}){$lipid{$uniID}="NA";}
		unless(exists $carbohyd{$uniID}){$carbohyd{$uniID}="NA";}
		unless(exists $disulfid{$uniID}){$disulfid{$uniID}="NA";}
		unless(exists $crosslnk{$uniID}){$crosslnk{$uniID}="NA";}
		unless(exists $transmem{$uniID}){$transmem{$uniID}="NA";}
		unless(exists $intramem{$uniID}){$intramem{$uniID}="NA";}
		unless(exists $repeat{$uniID}){$repeat{$uniID}="NA";}    
		unless(exists $cabind{$uniID}){$cabind{$uniID}="NA";} 
		unless(exists $dnabind{$uniID}){$dnabind{$uniID}="NA";}
		unless(exists $znfing{$uniID}){$znfing{$uniID}="NA";}
		unless(exists $npbind{$uniID}){$npbind{$uniID}="NA";}
		unless(exists $motif{$uniID}){$motif{$uniID}="NA";}
		unless(exists $signal{$uniID}){$signal{$uniID}="NA";}
		unless(exists $propep{$uniID}){$propep{$uniID}="NA";}
		unless(exists $transit{$uniID}){$transit{$uniID}="NA";}
	  for(my $i=0;$i<=$#uniprotAC; $i++){
		my $uniAC=$uniprotAC[$i];
#print "kaka $uniID $uniAC \n";
		printf FRES "$uniID\t$uniAC\t$binding{$uniID}\t$actsite{$uniID}\t$metal{$uniID}\t$site{$uniID}\t$nonstd{$uniID}\t$modres{$uniID}\t$lipid{$uniID}\t$carbohyd{$uniID}\n";
		printf PAIR "$uniID\t$uniAC\t$disulfid{$uniID}\t$crosslnk{$uniID}\n";
		printf RNGS "$uniID\t$uniAC\t$transmem{$uniID}\t$intramem{$uniID}\t$repeat{$uniID}\t$cabind{$uniID}\t$dnabind{$uniID}\t$znfing{$uniID}\t$npbind{$uniID}\t$motif{$uniID}\t$signal{$uniID}\t$propep{$uniID}\t$transit{$uniID}\n";
	  }
	  undef @uniprotAC;
	}	
}


sub open_file{
        my ($file_name)=@_;
        open(INP1, "< $file_name") or die "Can not open an input file: $!";
        my @file1=<INP1>;
        close (INP1);
        chomp @file1;
        return @file1;
}

