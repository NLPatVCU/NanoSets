#!/usr/bin/perl
=head1 NAME

nano-sets.pl - This Program returns a set of terms from a term queried

=head1 SYNOPSIS

This is a utility that takes as input either one term (DEFAULT)
or one CUI and the number of returns and returns a number of similar
concepts to that query.

=head1 USAGE

nano-sets.pl [OPTIONS] [TERM]

=head1 INPUT

=head2 [CUI|TERM]

The input is one term or one CUI associated to concepts in the UMLS.

=head1 OPTIONS

Optional command line arguements

=head2 General Options: 

=head3 --config FILE

This is the configuration file. There are a few configuration options
that can along with the UMLS::similarity packages config file options.

See umls-similarity documentation for more information.

=head3 --help

Displays a quick summary of program options and usage.

=head3 --version

Displays the version informations.

=head2 Input Options:

=head3 --infile FILE

A file containing a list of CUIS or terms in the following format:
Can be reversed: 

term1||cui1
term2||cui2
...
termN||cuiN

=head2 Debug Options:

=head3 --debug

Prints a log file containing progress and error information

=head 1 SYSTEM REQUIREMENTS

=over

=item * Perl (version 5.8.5 or better) - http://www.perl.org

=item * UMLS::Interface - http://search.cpan.org/dist/UMLS-Interface

=item * UMLS::Similarity - http://search.cpan.org/dist/UMLS-Similarity

=back

=head1 CONTACT US

    If you have any trouble installing and using UMLS-Similarity, 
    please contact us via the users mailing list :
    
      umls-similarity@yahoogroups.com
     
  You can join this group by going to:
    
      http://tech.groups.yahoo.com/group/umls-similarity/
     
  You may also contact us directly if you prefer :
    
      Bridget T. McInnes: bthomson at cs.umn.edu 
      


=head1 AUTHOR

 Robert S. VanDerzee
 
 Bridget T. McIness
 
=cut
###############################################################################

#                               THE CODE STARTS HERE
###############################################################################

#                           ================================
#                            COMMAND LINE OPTIONS AND USAGE
#                           ================================

                 
                                                                                
package UMLS::NanoSets;                                                         #NanoInformatics::NanoSets Package - Dr. McInnes -Robbie VanDerzee

#use strict;                                                                    #Strict-Rule -Disable when using Getopt::Long;
use warnings;                                                                   #Error-Warnings -Disable when experimenting
use Path::Class;                                                                #File-Pathing: Useful for file management
use autodie;                                                                    #Auto-or die for file editing
use List::Util qw( min max );                                                   #For calculating the max or minimum elements
use Getopt::Long;                                                               #Option-Config Handling - Useful for end-user and debug (Does not use strict)
use UMLS::Interface; 
    

#VERSION
our $VERSION = 0.01;                                                            #Version-Key::Modifier | Change if functionality is added. Document all changes with digital signature

#IMPLEMENTATION
my $getCUI          = UMLS::Interface->new();                                   #Import used to getAllCuis for the list of terms.
my $numb            = 0;
my $newterm         = "";
my $logdir          = "";

#   reference the Getoption cpan page
eval(GetOptions( "version", "help", "config=s", "term=s", "num=i", "log=s")) 
or die ("Please check the above mentioned option(s).\n");

#   if help is defined, print out help
if( defined $opt_help ) {
    $opt_help = 1;
    &showHelp;
    exit;
}

# if version is requested, show version
if( defined $opt_version ) {
    $opt_version = 1;
    &showVersion();
    exit;
}

if( defined $opt_config ) {
    $opt_config = 1;
    system "emacs config";
    exit;
    }
    
if( defined $opt_log ) {
    $logdir = $opt_log;
    chomp $logdir;
    
    if(-d $logdir) {
        } else {
            mkdir $logdir;
        }
    chdir ("$logdir");
}
    
print "--num NUMBER\n";
if(defined $opt_num) { 
    print "  You have defined --num as $opt_num\n\n";
    $numb = $opt_num;
}
else { 
    print "  You did not define the --num option\n\n";
    $numb = 6;
}

print "--term TERM\n";
if(defined $opt_term) { 
    print "  You have defined --term as $opt_term\n\n";
    $newterm = $opt_term;
}
else { 
    print "  You did not define the --string option\n\n";
    $newterm         = "Death";
}

##############################################################################                 

#  function to output help messages for this program                                         

##############################################################################                

sub showHelp() {

    print "This is a utility that provides an example of how to set\n";
    print "and use the GetOptions module for the programs written in\n"; 
    print "our lab -- this includes the help and showVersion information\n\n";

    print "Usage: NanoSets.pl [OPTIONS] \n\n";

    print "OPTIONS:\n\n";

    print "--num NUMBER           This option takes the max number of returned terms \n\n";

    print "--string STRING        This option takes the new term as a CUI or term \n\n"; 
    
    print "--log LOG              This option places all the data into the log director\n\n";

    print "--version                Prints the version number\n\n";

    print "--help                   Prints this help message.\n\n";
}

#  function to output the version number                                                   

##############################################################################              

sub showVersion {
    print '$Version: NanoSets.pl, 2016/08/010 12:50 btmcinnes';
    print "\nCopyright (c) 2015- Bridget McInnes\n"; 
}

#-----------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------

#TENTATIVE EXAMPLES->[USE THIS FOR AN ENTIRE PROGRAM RUN], OTHERWISE, IMPLEMENT YOUR OWN USE

        my @currentArray    = ();                                               #The current centroid array being used within a loop.
        my @scores          = ();                                               #An array containing every score within the umlsmatrix
        my $i               = 0;                                                #An incrementing variable with multiple purposes
        my $similarityscore = 0;                                                #A return value for saving the cosine_Similarity(@A, @B) return value
                                                                                #A term to create the score vector
        
        my @similarityscoresmetrics = ();
        my $SMQfile =
          retrieve_SMQS( "SMQ_Spreadsheet_flat_19_0.csv",
          "Test" );                                                             #File containing line by line $SMQHeader $Relationship $SMQ_Term
        my $newCUIfile = to_CUIs ($SMQfile);
        my $newScoreFile = score_NewTerm( $newterm, $newCUIfile, "termsvdata" );#A file matrix containing the score<>$newterm(CUI)<>$compared(CUI)
        my $umlsarray    = create_UMLSarray($newScoreFile);                     #New vector containing scores between it and every term in the database selected
        my $clusters     = get_Clusters("temp.clustering.3");                   #Returns the clusters foreach vector data
        my ( $scores, $length ) = get_Scores("temp");                           #Create a score array of every score in the matrix (WARNING::This may take long)                        
        my ( $centroids, $centroid ) =                                          #Create each centroid. This returns the centroid hashtable (sorted) (each key is a cluster, in order) and the newly created centroid file
          create_CentroidFile( $scores, $clusters,
          "index", "centroid.out", $length );  
        my %hashset = %$centroids;                                              #A copy of the centroids hashtable
        my @config = ();
           
        
#EXAMPLES --END
#-----------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------

my $errorhandler = "";
my $bestclustermeasure; 


foreach my $arrayRef ( sort keys %hashset )                                     #For every key in the hashset, that keys value is an array Reference
                                                                                #FUTURE VERSIONS:: REFERENCE THE HASHTABLE TO CONSERVE MEMORY
    {
    @currentArray = @${ $hashset{$arrayRef} };                                  #The current array being worked on is dereferenced as such
    $similarityscore = cosine_Similarity( $umlsarray, \@currentArray );         #Retrieve the cosine similarity between the new terms vector and the current array
    #print STDERR "\n Cluster # | $i | TERM | $newterm |"                       #For debugging or other (TENTATIVE) use.
    #  . " Cosine Similarity | $similarityscore\n";
    $similarityscoresmetrics[$i] = $similarityscore;                            #Conserves the score, and the index that score exists is its cluster identifier
    $i++;                                                                       #Incrememented to add to the @Scores array
    }
    for (my $y = 0; $y < scalar @similarityscoresmetrics; $y++)
    {
    if(max(@similarityscoresmetrics)==$similarityscoresmetrics[$y])
        {
        $bestclustermeasure = $y;
        last;
        }
    }
    
    returned_Terms ($newterm, "centroid.out", $bestclustermeasure, $numb);
    
=head3 score_NewTerm

Description:
Function to format a file to umls-similarity format

Input:

NewTerm CUIs NameOfOutPut

Output:

As file: NameOfOutPut.out
List by Line
NewTerm<>CUI1
NewTerm<>CUI2
    ...<>...
NewTerm<>CUIN

Example:

use UMLS::NanoSets;
my $umls = UMLS::NanoSets->new();
my $term = "Back Injury";
my $outfile;

$outfile = $umls->score_NewTerm($term, "CUIFILE", "OUTPUT");

=cut

sub score_NewTerm {
    my ( $newterm, $testCUIsfile, $outfile ) = @_;                              #Param 1, $newTerm, Param 2, $CUIs file, Param 3, Name of output file
    my $output;                                                                 #To write - file. Exists as $outfile name.out
    my $outtemp = "";                                                           #TEMPORARY -> Used to make $outfile a file object
    chomp( $newterm, $testCUIsfile, $outfile );                                 #Chomp any and all unnessecary spaces on all variable strings
    $outtemp = file($outfile);                                                  #Make a file object to being write. If it doesn't exists it will be created
    $output  = $outtemp->openw();                                               #Open and begin writting
    
    open( my $read, '<', $testCUIsfile ) or die "Failed to open $testCUIsfile"; 

    while ( my $line = <$read> )                                              
        {
        chomp $line;                                                           
        $output->print("$newterm<>$line\n");                                    #Print to $output, the newterm and the rest of the line (Which is a CUI)
        }                                                                    
        close $read;                                                            

    system                                                                      #A system call to get the UMLS-Similarity scores for each SMQ. 
"umls-similarity.pl --infile $outfile".
" --config /home/robbie/CONFIG > $outfile.out";
    return "$outfile.out";                                                      #Return this newly created file
}                                                                               #score_NewTerm -- end

=head3 to_CUIs

Description:
Function which converts terms to CUIS and removes them if no CUI is found

Input:

CUIFile

Output:

As file: CUIFile.out

List by Line
CUI1
CUI2
...
CUI3

Example:

use UMLS::NanoSets;
my $umls = UMLS::NanoSets->new();
my $outfile;

$outfile = $umls->to_CUIs("CUIFILE");

=cut

sub to_CUIs {                                                                   #A method that takes input of a file of terms and converts it into CUIS
    my ($CUIfile) = @_;                                                         #Parameter 1 is the file name
      chomp $CUIfile;                                                        
    my $CUIout = $CUIfile.".out";                                               #The output file name is equal to the file name plus ".out"                                                               
    
    my $CUIoutO = file($CUIout);                                              
    my $CUIoutW = $CUIoutO->openw();                                          
    
    open ( my $read, '<', $CUIfile) or die "Failed to open $CUIfile";          
    
    while ( my $line = <$read> )                                              
        {
            chomp $line;                                                      
            my $Cuis = $getCUI->getAllConcepts("$line");                        #Get CUI for the term on $line
            my $CUI = @{$Cuis}[0];                                              #The first CUI return is the CUI for that term
                if (length $CUI)                                                #If the return was initialized,
                {
                    $CUIoutW->print("$CUI\n");                                  #Print that CUI to file
                } else {
                    #print STDERR "'$line' was not found as a CUI\n";              #Print the line
                }
        }
        close $read;                                                            #Close the file
        print STDERR "\n";
    return $CUIout;                                                             #Return the newly written file as a String
    }                                                                           #to_CUIs --end

=head3 create_UMLSarray

Description:

Function to create a vector of scores between the new term and scores

Input:

UMLSFile

Output:

Array
1.000 -1.000 ... .1444
Example:

use UMLS::NanoSets;
my $umls = UMLS::NanoSets->new();
my $outarrayREF = ();

$outarrayREF = $umls->create_UMLSarray("UMLSFILE");

=cut

sub create_UMLSarray {                                                          #Create a new terms vector between that term<>Every cui in the database
    my ($umlsfile) = $_[0];                                                     #Umls file is parameter 0
    my $i          = 0;                                                         #Incrementor used for getting every 3rd element of the string
    my @split      = ();                                                        #Array used for splitting the line
    my @array      = ();                                                        #Array used for inputing the new set of information
    my $arrayIndex = 0;                                                         #Index value for moving along @array

    open( my $read, '<', "$umlsfile" ) or die "Failed to open $umlsfile";   
    
    while ( my $line = <$read> ) 
        {
        chomp $line;
        @split = split "<>", $line;

        foreach my $element (@split) 
        {
            if ( $i % 3 == 0 )                                                  #Every 3rd element in this file is the score 
            {
                $array[$arrayIndex] = $element;
                $arrayIndex++;
            }
            $i++;
        }
    }
    close $read;    

    return \@array;                                                             #Return that array as reference
}                                                                               #create_UMLSarray --end

=head3 retrieve_SMQS

Description:
A function that retrieves the SMQ concepts from the "SMQ_Spreadsheet_flat_XX_X.csv" release

Input:

SMQFILE OUTPUT

Output:

As file: OUTPUT

List by Line
CONCEPT1
CONCEPT2
...
CONCEPTN

Example:

use UMLS::NanoSets;
my $umls = UMLS::NanoSets->new();
my $outfile;

$outfile = $umls->retrieve_SMQS("SMQ_Spreadsheet_flat_19_0.csv", "SMQCONCEPTS");

=cut

sub retrieve_SMQS {                                                             #Creates a new file containing all the SMQ concepts from the SMQ file Spreadsheet
    my ( $SMQfile, $SMQfileW ) = @_;                                            #Para2 1 is the SMQ file, Param 2 is the new File 
    chomp( $SMQfile, $SMQfileW );
    my $Form = "";
    my $SMQHeader = "";
    my $file      = $SMQfileW;
    $SMQfileW = file($SMQfileW);
    $SMQfileW = $SMQfileW->openw();                                             #File opened to write

    open( my $read, '<', $SMQfile ) or die "Failed to open $file";

    while ( my $line = <$read> ) {
        chomp $line;
        if ( $line =~ /(SMQ)/ ) {                                               #Lines containing (SMQ)
            $line =~ s/,//g;                                                    #Globally remove ALL commas
            $SMQHeader = $line;
        }
        elsif ( $line =~ /\,+(Narrow|Broad)\,[A-Z]\,(.*?)\,[0-9]+/ ) {          #Array of index 1 any, Narrow or Broad, A-Z, All extra, Number
            $Form = $2;
            chomp $Form;
            $Form=~s/\"//g;                                                     #"
            $SMQfileW->print("$Form\n");                                        #Print information from file
        }
    }
    close $read;

    return $file;                                                               #Return new file
}                                                                               #retrieve_SMQS --end

=head3 create_Cluster

Description:
A system call to create the clustering file

Input:

UMLSMATRIXFILE NUMOFCLUSTERS METHODTYPE

Output:

As file: OUTPUT

List by Line
1
0
5
...
1

Example:

use UMLS::NanoSets;
my $umls = UMLS::NanoSets->new();
my $outfile;

$outfile = $umls->create_Cluster("UMLSMATRIX.FILE", "261", "vcluster");

=cut

sub create_Cluster {                                                            #Run vcluster and return that file created -- CLUTO FILE MUST BE BASHED
    my ( $umlsmatrix, $clusters, $method );
    system "$method" . "cluster $umlsmatrix $clusters";                         #Call for cluto from bash indexed vcluster.pl
    my $file = "$umlsmatrix.clustering.$clusters";                              #The file cluto will create

    return $file;                                                               #Create File
}                                                                               #create_Cluster --exit

=head3 get_Clusters

Description:
A function that reads the create_Clusters file and stores it

Input:

CLUTOFILE

Output:

ArrayREF

Example:

use UMLS::NanoSets;
my $umls = UMLS::NanoSets->new();
my $outREF;

$outREF = $umls->create_Cluster("CLUTOFILE");

=cut

sub get_Clusters {                                                              #Read the cluster file to get information about where terms are located
    my ($file) = $_[0];                                                         #The file is the first param
    my @array  = ();        
    my $i      = 0;
    open( my $read, '<', $file ) or die "Failed to open $file";

    while ( my $line = <$read> ) {
        chomp $line;
        $array[$i] = $line;                                                     #Adds that number to the array
        $i++;
    }
    close $read;                                                
    
    return \@array;                                                             #Returns that reference
}                                                                               #get_Clusters --end 

=head3 get_Scores

Description:
A system call to create the clustering file

Input:

UMLSFILE

Output:

As file: OUTPUT

ArrayREF LENGTHOFARRAY

Example:

use UMLS::NanoSets;
my $umls = UMLS::NanoSets->new();
my $outfile;
my $length;

($outfile, $length) = $umls->get_Scores("UMLSMATRIX.FILE");

=cut

sub get_Scores {                                                                #Get the scores from the umls file
    my ($umlsScoreFile) = $_[0];                                                #The first parameter is the umls score file
    chomp $umlsScoreFile;
    my @array  = ();
    my $length = 0;
    my $t      = 0;
    my $o      = 0;
    my @split  = ();

    open( my $read, '<', $umlsScoreFile )
      or die "Failed to open $umlsScoreFile";

    my $line = <$read>;
    @split = split " ", $line;
    $length = $split[0];

    while ( my $line = <$read> ) {
        chomp $line;
        @split = split " ", $line;
        for ( my $t = 0 ; $t < $length ; $t++ ) {
            $array[$o] = $split[$t];
            $o++;
        }
    }
    close $read;

    return ( \@array, $length );                                                #return the newly filled array as a reference and the length of that array
}

=head3 create_Centroid

Description:

A fucntion that creates a file containing important analysis information

Input:

SCOREFILE CLUSTERIDS INDEXFILE CENTROIDFILENAME LENGTHPERARRAY

Output:

CENTROIDS HASHTABLE As file: OUTPUT

List by Line
CUI INDEX CLUSTER
SCORES
CUI1 INDEX1 CLUSTER1
SCORES1
...
CUIN INDEXN CLUSTERN
SCORESN
END
CLUSTER ID
CENTROIDSCORESFORID
CLUSTER1 ID1
CENTROIDSCORESFORID1
...
CLUSTERIDN
CENTROIDSCORESFORIDN

Example:

use UMLS::NanoSets;
my $umls = UMLS::NanoSets->new();
my $outfile;

my ( $centroids, $centroid ) =                                        
          create_CentroidFile( $scores, $clusters,
          "index", "centroid.out", $length );

=cut

sub create_CentroidFile {                                                       #Create the centroid file. This contains a large amount of information, and calculates each clusters centroid
    
    my ( $scoresRef, $clusterNumRef, $indexFile, $centroidFile, $length ) = @_; #scoresRef is the array ref of scores, cluster ref is the clusters array, index is the index numbers, centroid file is where to write, length is the length of the array
    my $i              = 0;
    my $o              = 0;
    my @array          = @{$scoresRef};                                         #Array Derefencing 
    my @cluster        = @{$clusterNumRef};
    my %centroidarrays = ();
    my @count          = ();
    my $printVal       = "";
    chomp( $indexFile, $centroidFile, $length );

    my $centroidFileW = file($centroidFile);
    $centroidFileW = $centroidFileW->openw();

    open( my $read, '<', $indexFile );

    while ( my $line = <$read> ) {
        chomp $line;
        $centroidFileW->print("$line $cluster[$i] \n");                         #Print the CUI, Index, and ID

        for ( my $u = 0 ; $u < $length ; $u++ ) {
            $centroidFileW->print("$array[$o] ");                               #Print the scores until the length of the array is met. 
            $o++;
        }
        $centroidFileW->print("\n");
        $i++;
    }
    $centroidFileW->print("END");
    $o = 0;

    foreach my $element (@cluster) {                                            #For each cluster ID
        for ( my $u = 0 ; $u < $length ; $u++ ) {
            @${ $centroidarrays{$element} }[$u] += $array[$o];                  #Sum the elements for that ID
            $o++;                                                               #Next element in the array
        }
        $count[$element]++;
    }
    foreach my $key ( sort keys %centroidarrays ) {                             #For each key in the centroid arrays
        $centroidFileW->print("\nCluster # | $key | SIZE | $length\n");

        for ( my $j = 0 ; $j < $length ; $j++ ) {
            @${ $centroidarrays{$key} }[$j] =
              @${ $centroidarrays{$key} }[$j] / $count[$key];                   #Average the array values
            $printVal = @${ $centroidarrays{$key} }[$j];
            $printVal = pack( "A7", $printVal );
            $printVal =~ s/\s/0/g;
            $centroidFileW->print("$printVal ");
        }
    }
    close $read;

    return ( \%centroidarrays, $centroidFile );                                 #Return a reference to that hashtable and the file that was created. 
}

=head3 cosine_Similarity

Description:
Calculate the consine similarity between two vectors

Input:

VectorA VectorB

Output:

SCORE

Example:

use UMLS::NanoSets;
my $umls = UMLS::NanoSets->new();
my $score;

$score = $umls->get_Scores(\@A, \@B);

=cut

sub cosine_Similarity {                                                         #Utilize the euclidean dot product and magnitude to calculate the similarity between two vectors.
    my ( $a_ref, $b_ref ) = @_;
    my @arrayA = @{$a_ref};
    my @arrayB = @{$b_ref};

    my $c             = 0;
    my $cosine        = 0;
    my $aMagnitude    = 0;
    my $bMagnitude    = 0;
    my $EuProd        = 0;
    my $MagnitudeProd = 0;
    my $lowerBound    = 0;
    my $upperBound    = 0;
    $lowerBound = scalar @arrayA;
    $upperBound = scalar @arrayB;

    for ( my $c = 0 ; $c < max( $lowerBound, $upperBound ) ; $c++ ) {           #For the longer array
        if ( !$arrayA[$c] || $arrayA[$c] == -1 ) {
            $arrayA[$c] = 0;                                                    #null or -1 values reset
        }
        if ( !$arrayB[$c] || $arrayB[$c] == -1 ) {
            $arrayB[$c] = 0;
        }
        $EuProd     += $arrayA[$c] * $arrayB[$c];                               #product the values together
        $aMagnitude += $arrayA[$c]**2;                                          #Square the values and sum
        $bMagnitude += $arrayB[$c]**2;                                          #Square the values and sum
    }

    $aMagnitude = sqrt($aMagnitude);                                            #Sqrt to get magnitudes
    $bMagnitude = sqrt($bMagnitude);                                            #Sqrt to get magnitudes

    $MagnitudeProd = $aMagnitude * $bMagnitude;
    $cosine        = $EuProd / $MagnitudeProd;                                  #Divide to calculate
        
    return $cosine;                                                             #Return the decimal value
}

    
=head3 returned_Terms

Description:
Return terms given scores, the newterm, the cluster centroids, and how many desired returns. The final output method

Input:

NEWTERM CENTROIDFILE CLUSTERIDS NUMOFRETURNS

Output:

RANK | TERM
RANK | TERM1
RANK | TERM2
...
RANK | TERMN

Example:

use UMLS::NanoSets;
my $umls = UMLS::NanoSets->new();
my $finished;

$finished = returned_Terms ($newterm, "centroid.out", $scores[0], "100");
print $finished;

=cut

sub returned_Terms {                                                            #Return the terms based on the cluster chosen and the scores given
    my ( $newTerm, $centroidFile, $clusterID, $numOfReturns ) = @_;             #Param 1 new term, Param 2 Clusters, Param 3 Cluster ID, Param 4 num of returns
    my $terms = 0;
    my $i     = 0;
    chomp $newTerm;
    my @scores = ();
    my @terms2  = ();
    my @split = ();
    my $split2 = "";
    my $file = file("ClusterDataTerms");
    my $out = $file->openw();
    
    print STDERR "\n\n";
    
    open (my $read, '<', $centroidFile) or die "Failed to open $centroidFile"; 
    
    while ( my $line = <$read> ){
        chomp $line;
        
        if($line =~ "END")
            {
                last;
            }
            
        if($i % 2 == 0){
            
            @split = split " ", $line;
            if($split[2] == $clusterID)
                {
                $out->print("$newTerm<>$split[0]\n");
                }
            }
        $i++;    
        }
    close $read;    
        
    system "umls-similarity.pl --infile Data1".
    " --config /home/robbie/CONFIG > out";    
    
    open (my $read2, '<', "out") or die "Failed to open umls scores file, 'out'";
    $i = 0;
    while (my $line = <$read2>){
        @split = split "<>", $line;
        chomp $split[0];
        chomp $split[2];
        $scores[$i] = $split[0];
        $terms2[$i] = $split[2];
        $i++;
        }
    $i = 0;
    $terms = 0;
    
    print "Given the search query, $newTerm: \n\n";
    my %valueHash;
    
    for (my $s = 0; $s<scalar @terms2; $s++){
        push @{$valueHash{$scores[$s]}}, $terms2[$s]; 
        }
    
    my $h = 0;
    
    while($numOfReturns >= $h && $h < scalar @scores){
    foreach my $score (sort {$b <=> $a} keys %valueHash){
            foreach my $value (@{$valueHash{$score}}){
                if ($numOfReturns > $h && $h <= scalar @scores){
                    print "Rank | $h | Similarity | ".
                    "$score | Suggested Term | $value\n";
                }
                $h++; 
                }
            }
        }
    return "New terms compiled";                                                #Return that the operation has been completed. 
}
