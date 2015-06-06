#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

#Four files are needed: 1) list of variations 2-4) uniprotKB features (output of upanno.pl)

#control run options
my %opt;
getopts('hf:a:p:r:',\%opt);
die &usage() if (tutej_any($opt{f},$opt{a},$opt{p},$opt{r}));
&usage() if $opt{h};

#AA classification based on size and properties
my %size=(A => "t",G => "t",S => "t",T => "s",V => "s",C => "s",P => "s",N => "s",D => "s",I => "l", L => "l", M => "l", R => "l",H => "l",K => "l",E => "l",Q => "l",F => "l",Y => "l",W => "l");
my %charge=(H => "p",K => "p",R => "p",D => "n",E => "n",A => "z",G => "z",S => "z",T => "z",V => "z",C => "z",P => "z",N => "z",I => "z", L => "z", M => "z",Q => "l",F => "l",Y => "l",W => "l");
my %polarity=(C => "p",D => "p",E => "p",H => "p",K => "p",N => "p",P => "p",Q => "p",R => "p",S => "p",T => "p",A => "h", F => "h",I => "h",L => "h",V => "h",W => "h",M => "n",Y => "n",G => "n");

#alternative classification
#Class	dipole scale	V scale
#1	1	1	A,G,V
#2	1	2	I,L,F,P
#3	2	2	Y,M,T,S,C
#4	3	2	H,N,Q,W
#5	4	2	R,K
#6	5	2	D,E

#from	to	penalty comment
#1	1	0	similarly other "silent" mutations
#1	2	1
#1	3	2
#1	4	4
#1	5-6	5	
#2	3	2
#2	4	3
#2	5-6	4
#3	4	2
#3	5-6	3
#4	5-6	2
#5	6	3	

#read uniprot features
my ($uniAC,$uniID);
my (%binding,%actsite,%metal,%site,%nonstd,%modres,%lipid,%carbohyd);
my (%disulfid,%crosslnk);
my (%transmem,%intramem,%repeat,%cabind,%dnabind,%znfing,%npbind,%motif);
my (%signal,%propep,%transit);
if($opt{q}){
#AA
#uniprotAC;uniprotID;BINDING;ACT_SITE;LIPID;METAL;SITE;MOD_RES;CARBOHYD;SE_CYS
my @funcaa=open_file($opt{q});
for (my $i=1;$i<=$#funcaa;$i++){
  my $row=$funcaa[$i];
  my @prot=split(/\t/,$row);
  $uniID=$prot[0];
  $uniAC=$prot[1];
  $binding{$uniAC}=$prot[2];
  $actsite{$uniAC}=$prot[3];
  $metal{$uniAC}=$prot[4];
  $site{$uniAC}=$prot[5];
  $nonstd{$uniAC}=$prot[6];
  $modres{$uniAC}=$prot[7];
  $lipid{$uniAC}=$prot[8];
  $carbohyd{$uniAC}=$prot[9];
}
}
if($opt{p}){
#AA pairs, e.g. disulfides
#uniprotAC;uniprotID;DISULFID;THIOLEST;THIOETH
my @funcaapair=open_file($opt{p});
for (my $i=1;$i<=$#funcaapair;$i++){
  my $row=$funcaapair[$i];
  my @prot=split(/\t/,$row);
  $uniID=$prot[0];
  $uniAC=$prot[1];
  $disulfid{$uniAC}=$prot[2];
  $crosslnk{$uniAC}=$prot[3];
}
}
if($opt{r}){
#AA ranges e.g TM segments
#uniprotAC;uniprotID;TRANSMEM;SIGNAL;PROPEP
my @funcaarange=open_file($opt{r});
for (my $i=1;$i<=$#funcaarange;$i++){
  my $row=$funcaarange[$i];
  my @prot=split(/\t/,$row);
  $uniID=$prot[0];
  $uniAC=$prot[1];
  $transmem{$uniAC}=$prot[2];
  $intramem{$uniAC}=$prot[3];
  $repeat{$uniAC}=$prot[4];
  $cabind{$uniAC}=$prot[5];
  $dnabind{$uniAC}=$prot[6];
  $znfing{$uniAC}=$prot[7];
  $npbind{$uniAC}=$prot[8];
  $motif{$uniAC}=$prot[9];
  $signal{$uniAC}=$prot[10];
  $propep{$uniAC}=$prot[11];
  $transit{$uniAC}=$prot[12];
}
}

#read variants and compare with annotation
open(VAR, "< $opt{f}") or die "Can not open an input file: $!";
open (OUT, ">Annotated_prop.txt"); 
open (OUT1, ">Annotated_AA.txt"); 
#open (OUT2, ">Annotated_AApairs.txt"); 
#open (OUT3, ">Annotated_AAranges.txt"); 
printf OUT "uniprotAC\twtAA\tAApos\tmtAA\tviol_size\tviol_chrg\tviol_polr\tviol_start\n";
printf OUT1 "uniprotAC\twtAA\tAApos\tmtAA\tFUNC\tfuncAA\tfuncAApair\tfuncAArnge\tBINDING\tACT_SITE\tMETAL\tSITE\tNON_STD\tMOD_RES\tLIPID\tCARBOHYD\tDISULFID\tCROSSLNK\tTRANSMEM\tINTRAMEM\tREPEAT\tCA_BIND\tDNA_BIND\tZN_FING\tNP_BIND\tMOTIF\tSIGNAL\tPROPEP\tTRANSIT\n";

#my @variants=open_file($opt{f});
#for (my $i=1;$i<=$#variants;$i++){#GVkey, Symbol, Entrez,  Uniprot AC, Uniprot ID, Ref AA, Position, Var AA, disease
while(<VAR>){#for each data set
	my $var=$_;
	my ($func, $fa, $fp, $fr, $a1, $a2, $a3, $a4, $a5, $a6, $a7, $a8, $p1, $p2, $r1, $r2, $r3, $r4, $r5, $r6, $r7, $r8, $r9, $r10, $r11);
	$func = $fa = $fp = $fr = $a1 = $a2 = $a3 = $a4 = $a5 = $a6 = $a7 = $a8 = $p1 = $p2 = $r1 = $r2 = $r3 = $r4 = $r5 = $r6 = $r7 = $r8 = $r9 = $r10 = $r11 = 0;
	my ($sizeviol, $chargeviol, $polarviol, $startviol);	
	$sizeviol = $chargeviol = $polarviol = $startviol = 0;
	my @data=split("\t",$var);
	#my $pos=$aa[1]+$shift;
	my (@wt,@pos,@mt,@trans,@uniAC);
	$uniAC[0] = "null";
	if($opt{t} eq "v"){	
	  my (%info);
	  unless($var =~ /^#/){	
		##interpret info field
		my @inf=split(/;/,$data[7]);
		foreach my $field(@inf){
		  my @pair=split(/=/,$field);
		  $info{$pair[0]}=$pair[1];
		}
		@trans=split(/;/,$info{ET});
		my @vars=split(/;/,$info{VARAA});
		foreach my $ts(@trans){
			open(TMP, "> input") or die "Can not open an input file: $!";
			printf TMP "$ts\n";
			close(TMP);
			my $out = `R --no-save --args input < mapuniprot.R > output`;
			$out =~ s/[1] //;
			$out =~ s/"//g;
			push @uniAC, $out;
		}
		foreach my $v(@vars){
			$v=s/p.//;
			my @code=split(//,$v);
			my $ref=shift @code;
			my $alt=pop @code;
			my $aap=join("",@code);
			push @wt, $ref;
			push @mt, $alt;
			push @pos, $aap;
		}
	  }
	}
	elsif($opt{t} eq "b"){
		push @wt, $data[5];
		push @pos, $data[6];
		push @mt, $data[7];
		push @uniAC, $data[3];
	}
	
	#control types of amino acid changes
	for (my $i=1;$i<=$#wt;$i++){
	  if ($size{$wt[$i]} eq "t" and $size{$mt[$i]} eq "s"){$sizeviol = 1;} 
	  if ($size{$wt[$i]} eq "s" and $size{$mt[$i]} eq "t"){$sizeviol = -1;}
	  if ($size{$wt[$i]} eq "t" and $size{$mt[$i]} eq "l"){$sizeviol = 2;} 
	  if ($size{$wt[$i]} eq "l" and $size{$mt[$i]} eq "t"){$sizeviol = -2;}
	  if ($size{$wt[$i]} eq "s" and $size{$mt[$i]} eq "l"){$sizeviol = 1;} 
	  if ($size{$wt[$i]} eq "l" and $size{$mt[$i]} eq "s"){$sizeviol = -1;}
	  if ($charge{$wt[$i]} eq "p" and $charge{$mt[$i]} eq "n"){$chargeviol = -2;}
	  if ($charge{$wt[$i]} eq "n" and $charge{$mt[$i]} eq "p"){$chargeviol = 2;}
  	  if ($charge{$wt[$i]} eq "p" and $charge{$mt[$i]} eq "z"){$chargeviol = -1;}
	  if ($charge{$wt[$i]} eq "z" and $charge{$mt[$i]} eq "p"){$chargeviol = 1;}
	  if ($charge{$wt[$i]} eq "n" and $charge{$mt[$i]} eq "z"){$chargeviol = 1;}
	  if ($charge{$wt[$i]} eq "z" and $charge{$mt[$i]} eq "n"){$chargeviol = -1;}
	  if ($polarity{$wt[$i]} eq "p" and $polarity{$mt[$i]} eq "h"){$polarviol = 2;}
	  if ($polarity{$wt[$i]} eq "h" and $polarity{$mt[$i]} eq "p"){$polarviol = -2;}
	  if ($polarity{$wt[$i]} eq "p" and $polarity{$mt[$i]} eq "n"){$polarviol = 1;}
	  if ($polarity{$wt[$i]} eq "n" and $polarity{$mt[$i]} eq "p"){$polarviol = -1;}
	  if ($polarity{$wt[$i]} eq "h" and $polarity{$mt[$i]} eq "n"){$polarviol = -1;}
	  if ($polarity{$wt[$i]} eq "n" and $polarity{$mt[$i]} eq "h"){$polarviol = 1;}
	  #check translation start: Met in eucaryots; Ile in mitochondrion 
	  if($pos[$i] == 1){
		if($wt[$i] ne $mt[$i]){$startviol = 1;}
	  }
	  printf OUT "$trans[$i]\t$uniAC[$i]\t$wt[$i]\t$pos[$i]\t$mt[$i]\t$sizeviol\t$chargeviol\t$polarviol\t$startviol\n";
	}
	for (my $i=1;$i<=$#wt;$i++){	
	  #test single AA
	  unless($binding{$uniAC} eq "NA"){$a1=testAA($pos[$i],$binding{$uniAC});}
	  unless($actsite{$uniAC} eq "NA"){$a2=testAA($pos[$i],$actsite{$uniAC});}
	  unless($metal{$uniAC} eq "NA"){$a3=testAA($pos[$i],$metal{$uniAC});}
	  unless($site{$uniAC} eq "NA"){$a4=testAA($pos[$i],$site{$uniAC});}
	  unless($nonstd{$uniAC} eq "NA"){$a5=testAA($pos[$i],$nonstd{$uniAC});}
	  unless($modres{$uniAC} eq "NA"){$a6=testAA($pos[$i],$modres{$uniAC});}
	  unless($lipid{$uniAC} eq "NA"){$a7=testAA($pos[$i],$lipid{$uniAC});}
	  unless($carbohyd{$uniAC} eq "NA"){$a8=testAA($pos[$i],$carbohyd{$uniAC});}
	  if($a1 == 1 or $a2 == 1 or $a3 == 1 or $a4 == 1 or $a5 == 1 or $a6 == 1 or $a7 == 1 or $a8 == 1){$fa=1;}
	  #test AA pairs
	  unless($disulfid{$uniAC} eq "NA"){$p1=testAApair($pos[$i],$disulfid{$uniAC});}
	  unless($crosslnk{$uniAC} eq "NA"){$p2=testAApair($pos[$i],$crosslnk{$uniAC});}
	  if($p1 == 1 or $p2 == 1){$fp=1;}
	  #test AA ranges
	  unless($transmem{$uniAC} eq "NA"){$r1=testAArange($pos[$i],$transmem{$uniAC});}
	  unless($intramem{$uniAC} eq "NA"){$r2=testAArange($pos[$i],$intramem{$uniAC});}
	  unless($repeat{$uniAC} eq "NA"){$r3=testAArange($pos[$i],$repeat{$uniAC});}
	  unless($cabind{$uniAC} eq "NA"){$r4=testAArange($pos[$i],$cabind{$uniAC});}
	  unless($dnabind{$uniAC} eq "NA"){$r5=testAArange($pos[$i],$dnabind{$uniAC});}
	  unless($znfing{$uniAC} eq "NA"){$r6=testAArange($pos[$i],$znfing{$uniAC});}
	  unless($npbind{$uniAC} eq "NA"){$r7=testAArange($pos[$i],$npbind{$uniAC});}
	  unless($motif{$uniAC} eq "NA"){$r8=testAArange($pos[$i],$motif{$uniAC});}
	  unless($signal{$uniAC} eq "NA"){$r9=testAArange($pos[$i],$signal{$uniAC});}
	  unless($propep{$uniAC} eq "NA"){$r10=testAArange($pos[$i],$propep{$uniAC});}	
	  unless($transit{$uniAC} eq "NA"){$r11=testAArange($pos[$i],$transit{$uniAC});}	
	  if($r1 == 1 or $r2 == 1 or $r3 == 1 or $r4 == 1 or $r5 == 1 or $r6 == 1 or $r7 == 1 or $r8 == 1 or $r9 == 1 or $r10 == 1 or $r11 == 1){$fr=1;}
	  #general feature
	  if($fa == 1 or $fp == 1 or $fr ==1){$func=1;}
	  printf OUT1 "$uniAC[$i]\t$wt[$i]\t$pos[$i]\t$mt[$i]\t$func\t$fa\t$fp\t$fr\t$a1\t$a2\t$a3\t$a4\t$a5\t$a6\t$a7\t$a8\t$p1\t$p2\t$r1\t$r2\t$r3\t$r4\t$r5\t$r6\t$r7\t$r8\t$r9\t$r10\t$r11\n";
	}
}

sub testAA{
  my $pos = shift;
  my $string = shift;
#  my ($pos,$string)=@_;
  my $match = 0;
  my @values=split(/,/,$string);
  foreach my $val(@values){
    if($pos eq $val){$match=1;} 
  }
  return $match;
}

sub testAApair{
  my $pos = shift;
  my $string = shift;
  my $match = 0;
  my @values=split(/,/,$string);
  foreach my $val(@values){
    my @pair = split(/-/,$val);
    if($pos == $pair[0] or $pos == $pair[1]){$match=1;} 
  }
  return $match;
}

sub testAArange{
  my ($pos,$string)=@_;
  my $match = 0;
  my @values=split(/,/,$string);
  foreach my $val(@values){
    my @range = split(/-/,$val);
    if($pos >= $range[0] and $pos <= $range[1]){
    	my $dl=$pos-$range[0];
    	my $dr=$range[1]-$pos;
    	$match="1;".$dl.";".$dr;
    } 
  }
  return $match;
}

sub open_file{
        my ($file_name)=@_;
        open(INP1, "< $file_name") or die "Can not open an input file: $!";
        my @file1=<INP1>;
        close (INP1);
        chomp @file1;
        return @file1;
}

if ($#ARGV != 3) {die "Program used with parameters [list of variations] [uniprot funct. AA] [uniprot pair] [uniprot ranges]\n";}
sub usage(){
print STDERR << "EOF";
Usage: mapsnv.pl -f [vcf file] -a [annotation file] [arguments] 
 -h	: help message
 -q	: uniprot annotations for single aminoacids
 -p	: uniprot annotations for pairs
 -r	: uniprot annotations for ranges

EOF
exit;
}
sub tutej_any { ( grep $_, @_ ) < 3 }

