/* This batch implements a model test of a 12 rate NREV model versus a standard GTR and strandGTR model. the strand GTR model constrains rates so that AG := TC, not AG := GA as in the GTR model.




*/

/* user defined functions */

LoadFunctionLibrary ("libv3/UtilityFunctions.bf");
LoadFunctionLibrary ("libv3/IOFunctions.bf");
LoadFunctionLibrary ("libv3/models/parameters.bf");
LoadFunctionLibrary ("libv3/all-terms.bf");
LoadFunctionLibrary ("libv3/tasks/trees.bf");

utility.SetEnvVariable ("LIKELIHOOD_FUNCTION_OUTPUT",1);
utility.SetEnvVariable ("ACCEPT_ROOTED_TREES",TRUE);
utility.SetEnvVariable ("AUTOMATICALLY_CONVERT_BRANCH_LENGTHS",TRUE);


global alpha = 0.35;

category c = (4, EQUAL, MEAN,
				GammaDist(_x_,alpha,alpha),
				CGammaDist(_x_,alpha,alpha),
				0 ,
		  	    1e5,
		  	    CGammaDist(_x_,alpha+1,alpha)
		  	 );

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
NREVBiasTerms = { { "AC*", "", "AT*", "CA*", "CG*", "CT*", "GA*", "GC*", "GT*", "TA*", "TC*", "TG*" } };
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

DataSet 		ds = ReadDataFile("seq.fasta");

DataSetFilter nucFilter = CreateFilter (ds,1);
HarvestFrequencies (nucFreq,nucFilter,1,1,1);

/* 12 rate model */

fprintf ( stdout, "\nNREV\n" );
fprintf ( stdout,   "----\n" );

PopulateNucleotideModelMatrix ( "NREVMatrix" );
Model NREVModel = ( NREVMatrix, nucFreq, 1 );

fscanf (PROMPT_FOR_FILE, "Tree", givenTree);


Tree T = givenTree;
LikelihoodFunction lf_nrev = ( nucFilter, givenTree );
Optimize ( res_nrev, lf_nrev );

function maximumlikelihood ( LIKELIHOODArray& ) {

	LIKELIHOODArray = { { res_nrev[1][0]} };
	for (likelihood = 0; likelihood < Columns ( LIKELIHOODArray ); likelihood = likelihood + 1 ){
                 if(LIKELIHOODArray[likelihood] is maximum){
                     fprintf ( bestree,T);
                  }


      }fprintf ( bestree,T);


}


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

/*The given tree with the highest likelihood*/
fprintf ( "data.csv",res_nrev[1][0],"\n");

fprintf ("data1.nwk",T,"\n");









/*fprintf ( da,LIKELIHOODArray);*/
