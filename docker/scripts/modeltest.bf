/* This batch implements a model test of a 12 rate NREV model versus a standard GTR and strandGTR model. the strand GTR model constrains rates so that AG := TC, not AG := GA as in the GTR model.

Code by Wayne Delport, inspired by Darren P. Martin

wdelport@mac.com
19 February 2009

*/

/* user defined functions */
function PopulateNucleotideModelMatrix ( ModelMatrixName& ) {

	ModelMatrixName = { 4, 4 };
	modelString = "";
	modelString = 64;

	hv = 0;
	for ( h = 0; h < 4; h = h + 1 ) {
		for ( v = 0; v < 4; v = v + 1 ) {
			if ( h != v ) {
				modelString = ("ModelMatrixName["+h+"]["+v+"] := " + NREVBiasTerms [ hv ] + "t;\n");
				/* fprintf ( stdout, modelString ); */
				ExecuteCommands ( modelString );
				hv = hv + 1;
			}
		}		
	}
	return 0;
}


/* end of user defined functions */

/* ---------------------------- Model stuff -------------------------------------------- */


/* set some constraints */

INTERACTIVE = 1;


LIKELIHOOD_FUNCTION_OUTPUT = 3;

/* end of constraints */

/* define a nucleotide bias correction matrix AG := 1 */

/* built like this incase we want to scale to codon models with nuc bias terms */

global alpha = .35;
alpha:>0.01;alpha:<100;
category c = (4, EQUAL, MEAN, 
				GammaDist(_x_,alpha,alpha), 
				CGammaDist(_x_,alpha,alpha), 
				0 , 
		  	    1e25,
		  	    CGammaDist(_x_,alpha+1,alpha)
		  	 );
	
NREVBiasTerms = { { "c*AC*", "", "c*AT*", "c*CA*", "c*CG*", "c*CT*", "c*GA*", "c*GC*", "c*GT*", "c*TA*", "c*TC*", "c*TG*" } };
/* without * to remove constraints: see malemaLovesGoats */
ratesArray = { { "AC", "", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG" } };

global AC = 1;
global AT = 1;
global CA = 1;
global CG = 1;
global CT = 1;
global GA = 1;
global GC = 1;
global GT = 1;
global TA = 1;
global TC = 1;
global TG = 1;


/*SetDialogPrompt ("");*/
DataSet 		ds = ReadDataFile("/srv/shiny-server/rpnrm/0.fasta");




DataSetFilter nucFilter = CreateFilter (ds,1);
HarvestFrequencies (nucFreq,nucFilter,1,1,1);


fprintf ( stdout, "\nGTR\n" );
fprintf ( stdout,   "----\n" );
/* standard GTR */
CA := AC;
GA := 1;
GC := CG;
TA := AT;
TC := CT;
TG := GT;

PopulateNucleotideModelMatrix ( "GTRMatrix" );
Model GTRModel = ( GTRMatrix, nucFreq, 1 );

fscanf ("/srv/shiny-server/rpnrm/0.nwk", "Tree", givenTree);

Tree T = givenTree;
LikelihoodFunction lf_gtr = ( nucFilter, givenTree );
Optimize ( res_gtr, lf_gtr );
/*fprintf ( stdout, lf_gtr, "\n" );
fprintf ( stdout, res_gtr, "\n" );*/

/*fprintf ( stdout, "\nRate Matrix\n\n" );*/
hv = 0;
for ( h = 0; h < 4; h = h + 1 ) {
	for ( v = 0; v < 4; v = v + 1 ) {
		if ( v != h ) {
			if ( hv == 1 ) {
				fprintf ( stdout, "1.000 " );
				hv = hv + 1;
			}
			else {
				string = "";
				string*4;
				string = ratesArray [ hv ];
				ExecuteCommands ( "fprintf ( stdout, Format (" + string + ", 0, 3 ), \" \" );" );
				hv = hv + 1;
			}
		}
		else {
			fprintf ( stdout, "*     " );
		}
	}
	fprintf ( stdout, "\n" );
}
Export ( modelstr, GTRModel );
/*fprintf ( stdout, modelstr, "\n" );*/

/* forward loop in the spirit of Sergei's noMoreBush */
for ( malemaLovesGoats = 0; malemaLovesGoats < Columns ( ratesArray ); malemaLovesGoats = malemaLovesGoats + 1 ) {
	string = "";
	string * 32;
	string = ("ClearConstraints(" + ratesArray[malemaLovesGoats] + ");\n");
	/* fprintf ( stdout, string ); */
	ExecuteCommands ( string );
}



fprintf ( stdout, "\nstGTR\n" );
fprintf ( stdout,   "----\n" );
/* strand GTR */
CA := GT;
GA := CT;
GC := CG;
TA := AT;
TC := 1;
TG := AC;

PopulateNucleotideModelMatrix ( "stGTRMatrix" );
Model stGTRModel = ( stGTRMatrix, nucFreq, 1 );
Tree T = givenTree;
LikelihoodFunction lf_stgtr = ( nucFilter, givenTree );
Optimize ( res_stgtr, lf_stgtr );

hv = 0;
for ( h = 0; h < 4; h = h + 1 ) {
	for ( v = 0; v < 4; v = v + 1 ) {
		if ( v != h ) {
			if ( hv == 1 ) {
				fprintf ( stdout, "1.000 " );
				hv = hv + 1;
			}
			else {
				string = "";
				string*4;
				string = ratesArray [ hv ];
				ExecuteCommands ( "fprintf ( stdout, Format (" + string + ", 0, 3 ), \" \" );" );
				hv = hv + 1;
			}
		}
		else {
			fprintf ( stdout, "*     " );
		}
	}
	fprintf ( stdout, "\n" );
}

Export ( modelstr, stGTRModel );
/*fprintf ( stdout, modelstr, "\n" );*/	

/* forward loop in the spirit of Sergei's noMoreBush */
for ( malemaLovesGoats = 0; malemaLovesGoats < Columns ( ratesArray ); malemaLovesGoats = malemaLovesGoats + 1 ) {
	string = "";
	string * 32;
	string = ("ClearConstraints(" + ratesArray[malemaLovesGoats] + ");\n");
	/* fprintf ( stdout, string ); */
	ExecuteCommands ( string );
}

/* 12 rate model 

fprintf ( stdout, "\nNREV\n" );
fprintf ( stdout,   "----\n" );*/



PopulateNucleotideModelMatrix ( "NREVMatrix" );
Model NREVModel = ( NREVMatrix, nucFreq, 1 );
Tree T = givenTree;
LikelihoodFunction lf_nrev = ( nucFilter, givenTree );
Optimize ( res_nrev, lf_nrev );
/*fprintf ( stdout, lf_nrev, "\n" );*/

/* fprintf ( stdout, res_nrev, "\n" ); */ 

/*fprintf ( stdout, "\nRate Matrix\n\n" );*/
hv = 0;
for ( h = 0; h < 4; h = h + 1 ) {
	for ( v = 0; v < 4; v = v + 1 ) {
		if ( v != h ) {
			if ( hv == 1 ) {
				fprintf ( stdout, "1.000 " );
				hv = hv + 1;
			}
			else {
				string = "";
				string*4;
				string = ratesArray [ hv ];
				ExecuteCommands ( "fprintf ( stdout, Format (" + string + ", 0, 3 ), \" \" );" );
				hv = hv + 1;
			}
		}
		else {
			fprintf ( stdout, "*     " );
		}
	}
	fprintf ( stdout, "\n" );
}

Export ( modelstr, NREVModel );

gtr=2*(res_gtr[1][1]) - 2*res_gtr[1][0];
nrev6= 2*(res_stgtr[1][1]) - 2*res_stgtr[1][0];
nrev12= 2*(res_nrev[1][1]) - 2*res_nrev[1][0];
/*df<-data.frame(gtr,nrev6,nrev12);*/
/*fprintf ( "aic.csv","GTR\t"  ,gtr, "\nNREV6\t"  ,nrev6, "\nNREV12\t"  ,nrev12);*/
fprintf ( "aic.csv" ,gtr, "\n"  ,nrev6, "\n"  ,nrev12);

/*fprintf ( "aic.csv" , "GTR\n"  , "NREV6\n"  ,"NREV12\n");
fprintf ( "aic.csv",gtr  ,nrev6 ,nrev12);*/
/*seq=LAST_FILE_PATH; 
fprintf("MODELTEST RESULTS",seq);

fprintf ( "MODELTEST RESULTS", "\nLikelihood ratio tests\n" );
fprintf ("MODELTEST RESULTS", "----------------------\n" );

fprintf ( "MODELTEST RESULTS", "GTR vs NREV: LRT = ", 2*(res_nrev[1][0]-res_gtr[1][0]), "\tP = ", 1-CChi2 ( 2*(res_nrev[1][0]-res_gtr[1][0]), res_nrev[1][2]-res_gtr[1][2] ), "\n" );
	
fprintf ( "MODELTEST RESULTS", "stGTR vs NREV: LRT = ", 2*(res_nrev[1][0]-res_stgtr[1][0]), "\tP = ", 1-CChi2 ( 2*(res_nrev[1][0]-res_stgtr[1][0]), res_nrev[1][2]-res_stgtr[1][2] ), "\n" );

fprintf ("MODELTEST RESULTS", "\nAIC\n" );
fprintf ( "MODELTEST RESULTS", "---\n" );
	
fprintf ( "MODELTEST RESULTS", "GTR (all parameters): ", 2*(res_gtr[1][1]) - 2*res_gtr[1][0], "\n" );
fprintf ( "MODELTEST RESULTS", "GTR (rates only): ", 2*(res_gtr[1][2]) - 2*res_gtr[1][0], "\n" );
fprintf ("MODELTEST RESULTS", "stGTR (all parameters): ", 2*(res_stgtr[1][1]) - 2*res_stgtr[1][0], "\n" );
fprintf ( "MODELTEST RESULTS", "stGTR (rates only): ", 2*(res_stgtr[1][2]) - 2*res_stgtr[1][0], "\n" );
fprintf ( "MODELTEST RESULTS", "NREV (all parameters): ", 2*(res_nrev[1][1]) - 2*res_nrev[1][0], "\n" );
fprintf ( "MODELTEST RESULTS", "NREV (rates only): ", 2*(res_nrev[1][2]) - 2*res_nrev[1][0], "\n" ); */












