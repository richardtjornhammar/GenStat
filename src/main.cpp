//// C/C++ STUFF
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <math.h>
#include <cstdlib>
#include <cstring>
#include <random>

//	MMDB	STUFF
#include "mmdb/mmdb_manager.h"

//	GSL	STUFF
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>

//	XTC	STUFF
#include <time.h>
#include <float.h>
#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_trr.h"
#include "xdrfile/xdrfile_xtc.h"

void output_matrix( gsl_matrix * );

//	SIZE OF D IS 10E-11 to 10E-10 m2/s FOR BIOMOLECULES
//	DIFF EQ:: d/dt(phi) = D d^2/dt^2(phi)	FICKS SECOND
//
//	CONS	VALUE				UNIT
#define KBOL	1.38064889E-23		//	m2 kg s-2 K-1
#define RGAS	8.31446210		//	J K−1 mol−1
#define NAVO	6.022140857E23		//

// COMPILE:: g++ -std=c++11 main.cpp -lmmdb -lgsl -lblas -o rich_dyn
// WXTC	  :: g++ -std=c++11 src/*.cpp -lmmdb -lgsl -lblas -lxdrfile -o rich_dyn


std::pair< int , std::vector< std::string > >
argparser( std::pair < int, char ** > mp, std::vector<int>& v_set, 
	std::vector < std::pair < int, std::pair<std::string, std::string > > > opts_n_defaults ) {
	// ARGS OUTPUT ORDER CORRESPONDS TO OPTS_N_DEFAULTS FIRST VALUE
	std::vector< std::string > args;
	std::pair  < int, std::vector< std::string > > ret_args;
	//std::vector< int > v_set;

	for( int i=0 ; i<mp.first ; i++ ) {
		args.push_back			( mp.second[i] );
	}
	for( int i=0 ; i<opts_n_defaults.size() ; i++ ) {
		ret_args.second.push_back	( opts_n_defaults[i].second.second );
		//v_set.push_back(0);
	}
	// COMMAND INPUT
	int arg 	= 0 ;
	ret_args.first 	= 0 ;
	int a_size = mp.first==1?2:mp.first;
	int help=0,failed=0;
	while ( ++arg < a_size  ) {
		if(  args[arg] == "-h" ||  args[arg] =="--h" || args[arg] =="--help" || args[arg] =="-help" || mp.first==1) {
			help=1;
			break;
		}else{
			int ia=arg,iap=++arg;
			for(int i=0;i<opts_n_defaults.size();i++) {
				if ( args[ia] == opts_n_defaults[i].second.first ) {
					ret_args.second[i] = args[iap];
					v_set[i]=1;
				}
			}
		}
	}
	for(int i=0;i<v_set.size();i++) {
		if( opts_n_defaults[i].first > 0 && v_set[i] == 0 ) {
			ret_args.first = 1;
			failed=1;
			break;
		}
	}
	if ( help==1 || failed==1 ) {
		std::cout << "HELP:: VIABLE INPUT OPTIONS ARE:\nHELP:: [X] \t OPTION \t VALUE\n" << std::endl;
		for(int i=0;i<opts_n_defaults.size();i++)
			std::cout << "HELP:: [" << opts_n_defaults[i].first << "] \t " 
				<< opts_n_defaults[i].second.first << " \t " 
				<< ret_args.second[i] << std::endl;
		ret_args.first = 1;
		std::cout << std::endl;	
		if(failed==1) {
			std::cout << "ERROR::OBLIGATORY ARGUMENTS NOT SET ( MARKED [1] )\n" << std::endl;
		}
	}
	//v_set.swap( vis );
	//std::cout << " " << v_set.size() << " " << vis.size() << std::endl;
	return ret_args;
}
void
fatal(void){
	std::string	author("Richard Tjörnhammar (e-mail: richard.tjornhammar@gmail.com)");
	std::cout << "INFO:: PROGRAM FAILED" << std::endl;
	std::cout << "PLEASE CONTACT " << author << "\nWITH ANY QUESTIONS REGARDING THE FUNCTIONALITY" << std::endl;
	exit(1);
}

void
fatal( std::string errstr ){
	std::string	author("Richard Tjörnhammar (e-mail: richard.tjornhammar@gmail.com)");
	std::cout << "INFO:: PROGRAM FAILED" << std::endl;
	std::cout << "INFO:: " << errstr << std::endl;
	std::cout << "PLEASE CONTACT " << author << "\nWITH ANY QUESTIONS REGARDING THE FUNCTIONALITY" << std::endl;
	exit(1);
}

void
die( std::string errstr ){
	fatal(errstr);
}

void write_atoms_to_xtc( gsl_matrix *cell , CMMDBManager *mmdb , std::string fpname , std::string ftname, std::vector<int> vimnr, int debug) 
{
	std::string testfn(ftname);
	std::string fname(ftname);
	char *testncfstr = new char[testfn.length() + 1];
	std::strcpy(testncfstr, testfn.c_str());

	int bytes = fpname.length();
	char * cpstrp	= new char[bytes + 1];
	std::strcpy (cpstrp, fpname.c_str() );

	int selHnd	= mmdb->NewSelection();
	int nAtoms	= 0;
	PPCAtom cur_atoms;

//	HERE BUILD A NEW MMDB MANAGER 
//	AND OUTPUT FIRST FRAME TO PDB 

	CMMDBManager	mmdb_single; 
	PPCModel	model;
	int nModels;
	mmdb -> GetModelTable ( model , nModels );
	double a, b, c, alf, bet, gam, cell_vol;
	int ivol;
	if (mmdb -> isCrystInfo( ) ) 
		mmdb -> GetCell	( a , b , c , alf , bet , gam , cell_vol ,ivol );
	if (mmdb -> isSpaceGroup() )  {
		char *sg_str = mmdb -> GetSpaceGroup();
		int RC = mmdb_single.SetSpaceGroup ( sg_str );
	}

	mmdb_single . DeleteAllModels (  );
	mmdb_single . SetCell	( a , b , c , alf , bet , gam , ivol );
	mmdb_single . AddModel	( model[vimnr[0]] );
	selHnd	= mmdb_single . NewSelection( );
	mmdb_single . Select	( selHnd,STYPE_ATOM,0,"*",
                		  ANY_RES,"*",ANY_RES,"*",
				  "HOH", "*", "*" , "*", SKEY_NEW );
	mmdb_single . GetSelIndex ( selHnd, cur_atoms, nAtoms );
	for (int i=0;i<nAtoms;i++)
		delete cur_atoms[i];
	mmdb_single . FinishStructEdit();
	int rval	= mmdb_single.WritePDBASCII( cpstrp );	//	READ PDB
	if(rval!=0) { fatal(); }
//
//	GMX COORDINATES ARE IN NM NOT AA
//
	selHnd	= mmdb -> NewSelection( );
	XDRFILE *xd;
	xd = xdrfile_open( fname.c_str(), "w" );
	int step1=0,result;
	float time1=0.0,prec1=10000;
	float gmx_cell[DIM][DIM];

	for(int i=0;i<DIM;i++)
		for(int j=0;j<DIM;j++)
			gmx_cell[i][j]=gsl_matrix_get(cell,i,j)*0.1;
	if (NULL == xd)
		die("Opening xdrfile for writing");
	for( int k=1; k<vimnr.size(); k++ ) {
		if(vimnr[k]==0)
			continue;
		mmdb->Select (	selHnd  , STYPE_ATOM , vimnr[k] , "*",
				ANY_RES , "*" , ANY_RES,"*" ,
				"*", "*", "*" , "*", SKEY_NEW );
		mmdb->Select (	selHnd  , STYPE_ATOM , vimnr[k] , "*",
				ANY_RES , "*" , ANY_RES ,"*" ,
				"HOH", "*", "*" , "*", SKEY_CLR );
		mmdb->GetSelIndex ( selHnd, cur_atoms, nAtoms );
		rvec *x;
		x = (float (*)[3]) calloc(nAtoms, sizeof(*x) );
		for(int iatom=0;iatom<nAtoms;iatom++) { // fill the damn thing...
			x[iatom][0] = cur_atoms[iatom]->x*0.1;
			x[iatom][1] = cur_atoms[iatom]->y*0.1;
			x[iatom][2] = cur_atoms[iatom]->z*0.1;
		}
		result = write_xtc(xd,nAtoms,step1+k,time1+k,gmx_cell,x,prec1);
		free(x);
		if ( result != 0 )
			die("Writing xtc file");
	}
	xdrfile_close(xd);

	if( debug ) {
		int natoms2,step2;
		rvec *x2;
		float toler=1e-2,time2,prec2; // toler 1e-3 nm is 1e-2 AA
		float gmx_cell2[DIM][DIM];
		result = read_xtc_natoms( testncfstr , &natoms2 );
		if ( exdrOK != result )
			die("read_xtc_natoms");
		if ( natoms2 != nAtoms ) 
			die("Number of atoms incorrect when reading xtc");
		x2 = (float (*)[3]) calloc(natoms2, sizeof(*x2) );
		if (NULL == x2)
			die("Allocating memory for x2");
		xd = xdrfile_open(testfn.c_str(),"r");
		if (NULL == xd)
			die("Opening xdrfile for reading");
		int k = 0;
		do {
			result = read_xtc(xd,natoms2,&step2,&time2,gmx_cell2,x2,&prec2);
			if (exdrENDOFFILE != result)
			{
				if (exdrOK != result)
					die("read_xtc");
				if (natoms2 !=nAtoms)
					die("natoms2 != natoms1");
				if (step2-step1 != k)
					die("incorrect step");
				if (fabs(time2-time1-k) > toler)
					die("incorrect time");
				if (fabs(prec2-prec1) > toler)
					die("incorrect precision");
				for(int i=0; (i<DIM); i++)
					for(int j=0; (j<DIM); j++)
						if (fabs(gmx_cell2[i][j] - gmx_cell[i][j]) > toler)
							die("box incorrect");
			}
			if( debug>1 ) { // check coordinates in first frame
				mmdb->Select (	selHnd  , STYPE_ATOM , vimnr[k] , "*",
						ANY_RES , "*"	, ANY_RES,"*" ,
						"*", "*", "*" 	, "*", SKEY_NEW );	
				mmdb->Select (	selHnd  , STYPE_ATOM , vimnr[k] , "*",
						ANY_RES , "*" 	, ANY_RES ,"*" ,
						"HOH", "*", "*"	, "*", SKEY_CLR );
				mmdb->GetSelIndex ( selHnd, cur_atoms, nAtoms );
	
				for(int iatom=0;iatom<nAtoms;iatom++) {
					std::cout << "CRD::" << x2[iatom][0] << " " << x2[iatom][1] << " " << x2[iatom][2] ;
					float res	=	x2[iatom][0] - cur_atoms[iatom]->x*0.1 +
								x2[iatom][1] - cur_atoms[iatom]->y*0.1 +
								x2[iatom][2] - cur_atoms[iatom]->z*0.1 ;
					if(res > toler)
						std::cout << " WARNING " << std::endl;
					else
						std::cout << std::endl;
				}
			}
			k++;
			std::cout << "*";
		} while (result == exdrOK);
		std::cout << " \nALL IS GOOD " << std::endl;
	}
//	SHOULD ALSO CHECK THE COORDINATES
}

static void test_xtc()
{ // GENERATES SOME BOGUS FOR TESTING VISUALIZATION ...
	std::string testfn("test.xtc");
	char *testncfstr = new char[testfn.length() + 1];
	strcpy(testncfstr, testfn.c_str());

	XDRFILE *xd;
	int result,i,j,k,nframes=13;
	int natoms2,natoms1=173;
	int step2,step1=1993;
	float time2,time1=1097.23;
	matrix box2,box1;
	rvec *x2,*x1;
	float prec2,prec1=1000;
	float toler=1e-3;
	
	printf("INFO:: TESTING XTC FUNCTIONALITY \t");
	for(i=0; (i<DIM); i++)
		for(j=0; (j<DIM); j++)
			box1[i][j] = (i+1)*3.7 + (j+1);
	x1 = (float (*)[3]) calloc(natoms2, sizeof(*x1) );
	if (NULL == x1)
		die("Allocating memory for x1 in test_xtc");
	
	for(i=0; (i<natoms1); i++)
		for(j=0; (j<DIM); j++)
			x1[i][j] = (i+1)*3.7 + (j+1);
	xd = xdrfile_open(testfn.c_str(),"w");
	if (NULL == xd)
		die("Opening xdrfile for writing");
	for(k=0; (k<nframes); k++)
		{
			result = write_xtc(xd,natoms1,step1+k,time1+k,box1,x1,prec1);
			if (0 != result)
				die("Writing xtc file");
		}
	xdrfile_close(xd);
	
	result = read_xtc_natoms( testncfstr , &natoms2 );
	if (exdrOK != result)
		die("read_xtc_natoms");
	if (natoms2 != natoms1)
		die("Number of atoms incorrect when reading trr");
	x2 = (float (*)[3]) calloc(natoms2, sizeof(*x2) );
	if (NULL == x2)
		die("Allocating memory for x2");
		
	xd = xdrfile_open(testfn.c_str(),"r");
	if (NULL == xd)
		die("Opening xdrfile for reading");
	
	k = 0;
	do{
		result = read_xtc(xd,natoms2,&step2,&time2,box2,x2,&prec2);
		if (exdrENDOFFILE != result) {
			if (exdrOK != result)
				die("read_xtc");
			if (natoms2 != natoms1)
				die("natoms2 != natoms1");
			if (step2-step1 != k)
				die("incorrect step");
			if (fabs(time2-time1-k) > toler)
				die("incorrect time");
			if (fabs(prec2-prec1) > toler)
				die("incorrect precision");
			for(i=0; (i<DIM); i++)
				for(j=0; (j<DIM); j++)
					if (fabs(box2[i][j] - box1[i][j]) > toler)
						die("box incorrect");
			for(i=0; (i<natoms1); i++)
				for(j=0; (j<DIM); j++)
					if (fabs(x2[i][j] - x1[i][j]) > toler)
						die("x incorrect");
		}
		k++;
	} while (result == exdrOK);	
	xdrfile_close(xd);

	printf("[ PASSED ]\n");
}
void output_matrix( gsl_matrix *mat) {
	std::cout << "\n[========================] " << std::endl;
	for(int i=0;i<mat->size1;i++) {
		for(int j=0;j<mat->size2;j++) {
			std::cout << " " << gsl_matrix_get(mat,i,j);
		}
		std::cout << "\n";
	}
	std::cout << "[========================] " << std::endl;
}
	
void calc_cell_matrix( gsl_matrix *cell , CMMDBManager   *mmdb ) {
	if( !(cell->size1==DIM&&cell->size2==DIM) )
		fatal();
	realtype a, b, c, alf, bet, gam, cell_vol;
	realtype d2r=M_PI/180.0;
	int ivol;
	if(mmdb->isCrystInfo()) 
		mmdb -> GetCell( a , b , c , alf , bet , gam , cell_vol ,ivol );
	else
		fatal();
	std::cout << "GOT CELL: a, b, c, ang1, ang2, ang3 = " 
		<< a << ", " << b << ", " << c << ", " << alf << ", " << bet << ", " << gam << std::endl;
	alf*=d2r; bet*=d2r; gam*=d2r;

	gsl_matrix_set(cell,0,0,a);	gsl_matrix_set(cell,1,0, b*cos(gam) 	);	
	gsl_matrix_set(cell,0,1,0);	gsl_matrix_set(cell,1,1, b*sin(gam)	);
	gsl_matrix_set(cell,0,2,0);	gsl_matrix_set(cell,1,2, 0 		);

	double c21=c*(cos(alf)-cos(gam)*cos(bet))/sin(gam);
	gsl_matrix_set(cell, 2, 0, c*cos(bet) 	);
	gsl_matrix_set(cell, 2, 1, c21 );
	gsl_matrix_set(cell, 2, 2, sqrt( c*c*(1-cos(bet)*cos(bet)) - c21*c21 ) ); //check
}

int main ( int argc, char ** argv ) {
	int verbose = 0;

//	PREPARING ARGPARSER
	std::pair < int, char ** >	margs;
	margs.first = argc; margs.second = argv;

	std::vector < std::pair < int, std::pair<std::string, std::string > > >	opts_n_defaults;
	std::pair	< int, std::pair<std::string,std::string > > 		vipss;
//	WOULD BE BETTER WITH HASH INSTEAD SO WE CAN FORGET THE ORDER...
	vipss.first=1;	vipss.second.first	= "-ifile"; 
			vipss.second.second	= "i_default.pdb";
	opts_n_defaults.push_back(vipss);
	vipss.first=0;	vipss.second.first	= "-ofile"; 
			vipss.second.second	= "o_stationary.pdb";
	opts_n_defaults.push_back(vipss);
	vipss.first=0;	vipss.second.first	= "-nsteps";
			vipss.second.second	= "100";
	opts_n_defaults.push_back(vipss);
	vipss.first=0;	vipss.second.first	= "-verbose";
			vipss.second.second	= "0";
	opts_n_defaults.push_back(vipss);
	vipss.first=0;	vipss.second.first	= "-traj";
			vipss.second.second	= "traj.xtc";
	opts_n_defaults.push_back(vipss);
	vipss.first=0;	vipss.second.first	= "-frame";
			vipss.second.second	= "frame.pdb";
	opts_n_defaults.push_back(vipss);
	vipss.first=0;	vipss.second.first	= "-dt";
			vipss.second.second	= "1.0";
	opts_n_defaults.push_back(vipss);
	// test_xtc();
	if(argc<2)
		fatal("NOT ENOUGH ARGUMENTS");

//	EXECUTING ARGPARSE
	std::vector<int> v_set;
	for( int i=0 ; i < opts_n_defaults.size() ; i++ )
		v_set.push_back(0);
	std::pair< int , std::vector< std::string > > run_args = argparser( margs, v_set, opts_n_defaults );
	if(run_args.first) { // FAILED
		fatal();
	}

//	HERE WE SHOULD SET SOME INPUT PARAMETERS
	double dt	= 1.0E-12;				// TRAJECTORY TIMESTEP
	double gen_time	= atof(run_args.second[2].c_str())*dt;
	std::cout << "INFO:: WILL MIMIC SIMULATION TIME OF "<< gen_time*1E12 << " [ ps ] " << std::endl;
	verbose		= atoi(run_args.second[3].c_str());
	std::cout << "INFO:: SETTING VERBOSE LEVEL "<< verbose << std::endl;

//	THIS TOOL CURRENTLY ONLY FOR PDB
	CMMDBManager   mmdb;
	int itype, otype;
	itype = MMDB_FILE_PDB;

//	INITIALIZE MMDB
	InitMatType();
//
// 	SETTING SLACK FLAGS
//
	mmdb.SetFlag(	MMDBF_IgnoreDuplSeqNum | 
                    	MMDBF_IgnoreBlankLines |
                        MMDBF_IgnoreRemarks    |
                        MMDBF_IgnoreHash    |
                        MMDBF_IgnoreNonCoorPDBErrors |
                        MMDBF_IgnoreSegID |
                        MMDBF_IgnoreElement |
                        MMDBF_IgnoreCharge);
//
// 	READ COORDINATE FROM USING THE SUPPLIED NAME
//
	std::cout << "INFO:: WILL READ FILE " << run_args.second[0] << std::endl;
	int bytes = run_args.second[0].length();
	char * cpstr	= new char[bytes + 1];
	std::strcpy (cpstr, run_args.second[0].c_str() );
	int rval	= mmdb.ReadPDBASCII( cpstr );	//	READ PDB
	if(rval!=0) {
		fatal();
	}

//	START TO WORK WITH THE DATA (TABLES)
	PPCModel	model_T;
	PPCChain	chain_T1, chain_T2;
	PPCResidue	resid_T1, resid_T2;
	PPCAtom		atoms_T1, atoms_T2;
	int		nModels;
	int 		nChains1,nResidues1,nAtoms1;
	int 		nChains2,nResidues2,nAtoms2;
	int		imod,jmod,icha,ires,iat;

	mmdb.GetModelTable( model_T, nModels );
	std::cout << "INFO:: HAVE " << nModels << " MODELS" << std::endl;
	if(nModels <= 1 )
		fatal("NO PDB CONTAINS TO FEW MODELS");

	int nErr=0;
	
	gsl_matrix *M = gsl_matrix_calloc(nModels,nModels); // MSD
	gsl_matrix *W = gsl_matrix_calloc(nModels,nModels); // WEIGTHS
	gsl_matrix *P = gsl_matrix_calloc(nModels,nModels); // PROB
	gsl_matrix *S = gsl_matrix_calloc(nModels,nModels); // CUMULATIVE SUM

	gsl_matrix *cell = gsl_matrix_calloc(DIM,DIM);
	calc_cell_matrix( cell , &mmdb );
	// NOTE EACH MATRIX ROW IS A CELL VECTOR I IS ROW, J IS COLUMN. SEE OUTPUT
	output_matrix(cell);

	double gnrm	= 1.0/sqrt(2.0*M_PI);
	double l_jump	= 0.15E-9;	// JUMP DIFFUSION LENGTH OF PENTANE ( TYPICALLY 1-3 Å )

	double stationary_dist[nModels][2];
	for ( imod=0 ; imod<nModels ; imod++ ) {
		stationary_dist  [imod][0] = 0.0;
		stationary_dist  [imod][1] = 0.0;
		for ( jmod=imod ; jmod<nModels ; jmod++ ) {
			double X = 0.0, X2 = 0.0, Nc = 0.0;
			//	NOTE THAT TWO MODELS CAN HAVE DIFFERENT AMOUNT OF ATOMS
			//	ONLY USE THE ATOMS THAT ARE UNIQUELY DETERMINED IN BOTH
			//	SETS
			nChains1 = mmdb.GetNumberOfChains( imod+1 ); 
			nChains2 = mmdb.GetNumberOfChains( jmod+1 ); 
			if( nChains1 == nChains2 ) {
				// SAME PROTEIN IN DIFFERENT MODELS SO GO CHAIN BY CHAIN
				for ( icha = 0 ; icha < nChains1 ; icha++ ) { 
					nResidues1 = mmdb.GetNumberOfResidues( imod+1 , icha ); 
					nResidues2 = mmdb.GetNumberOfResidues( jmod+1 , icha ); 
					if ( nResidues1 == nResidues2 ) {
						for ( ires = 0 ; ires < nResidues1 ; ires++ ) { 
							mmdb.GetAtomTable    ( imod+1 ,icha ,ires , atoms_T1, nAtoms1 );
							mmdb.GetAtomTable    ( jmod+1 ,icha ,ires , atoms_T2, nAtoms2 );
							if( nAtoms1 == nAtoms2 ) { // MIGHT BE WATER
								// HERE WE CALCULATE
								for (iat = 0; iat<nAtoms1 ; iat++ ) {
									double dx  = atoms_T1[iat]->x - atoms_T2[iat]->x;
									double dy  = atoms_T1[iat]->y - atoms_T2[iat]->y;
									double dz  = atoms_T1[iat]->z - atoms_T2[iat]->z;
									double dr2 = dx*dx+dy*dy+dz*dz;
									double dr  = sqrt(dr2);
									X+=dr; X2+=dr2; Nc+=1.0;
								}
							}
						}
					} else {
						continue;
					}
				}
			} else {
				std::cout << "ERROR:: CHAIN MISSMATCH " << nChains1 << " AND " << nChains2 << std::endl;
				continue;
			}
			double Mij = X2/Nc-X*X/Nc/Nc;				// MSD
			double Wij = exp( -0.5*(Mij*1E-20)/l_jump/l_jump ); 	// WEIGHT ACCORDING TO JUMP DIFF PROB
			gsl_matrix_set( M, imod, jmod, Mij ); gsl_matrix_set( M, jmod, imod, Mij );
			gsl_matrix_set( W, imod, jmod, Wij ); gsl_matrix_set( W, jmod, imod, Wij );
			if( verbose==1 )
				std::cout << "INFO ( \t "<< imod << " , " << jmod << " ) = \t M = " 
					<< Mij << " \t W = " 
					<< Wij << std::endl;
		}
	}
	// CALCULATE MARKOV MATRIX FROM THE WEIGHTS
	gsl_vector *col = gsl_vector_calloc(nModels);
	gsl_vector *row = gsl_vector_calloc(nModels);
	gsl_vector *csr = gsl_vector_calloc(nModels);

	for ( imod=0 ; imod<nModels ; imod++ ) {
		gsl_matrix_get_row ( row , W , imod );
		double rsum = gsl_blas_dasum( row );
		double csum = 0.0;
		double rcur = 0.0;
		gsl_vector_scale( row , 1.0/rsum );
		for(int ir=0;ir<nModels;ir++) {
			rcur  = gsl_vector_get(row,ir);
			csum += rcur;
			gsl_vector_set(csr,ir,csum);
		}
		gsl_vector_scale( csr  ,  1.0/csum  );
		gsl_matrix_set_row ( S , imod , csr );
		gsl_matrix_set_row ( P , imod , row );
	}
	switch(verbose){
		case 3:
	// CHECK THE ROWS
			for ( imod=0 ; imod<nModels ; imod++ ) {
				gsl_matrix_get_row ( row , P , imod );
				gsl_matrix_get_row ( csr , S , imod );
				double csum = gsl_blas_dasum( csr );
				double rsum = gsl_blas_dasum( row );
				std::cout << "INFO " << rsum << " \t csr( " << csum << " ) "<< std::endl; 
			}
		case 2:
			for ( imod=0 ; imod<nModels ; imod++ ) {
				for ( jmod=0 ; jmod<nModels ; jmod++ ) {
					std::cout << gsl_matrix_get(P,imod,jmod) << " ";
				}
				std::cout << std::endl;
			}
		case 1:
			for ( imod=0 ; imod<nModels ; imod++ ) {
				for ( jmod=0 ; jmod<nModels ; jmod++ ) {
					std::cout << gsl_matrix_get(S,imod,jmod) << " ";
				}
				std::cout << std::endl;
			}
		default: break;
	}

	// GENERATE MARKOV TRAJECTORY 
	int step	= 0;
	double s_time	= 0.0;

	std::mt19937 generator(1871865);
	std::uniform_real_distribution<double> distribution(0.0,1.0);
	double dice_roll = distribution(generator),dModels=1.0;
	dModels=((double)nModels);			
	int I0 = floor( nModels*dice_roll ); 
	int I  = 0;
	std::cout << std::endl;
	std::cout << "===========================================" << std::endl;
	std::cout << "WILL START TO GENERATE THE TRANSITION ORDER" << std::endl;
	std::cout << "===========================================" << std::endl;
	double tot_cnt=0.0;
	std::vector<int> vimnr;
	while( s_time < gen_time ) {
		s_time+=dt;
		gsl_matrix_get_row( row , S , I0 );
		do {	I 		= floor(distribution(generator)*dModels);
			dice_roll 	= 	distribution(generator);
		} while( gsl_vector_get(row,I) < dice_roll );
		vimnr.push_back(I);
		if(verbose>=1)
			std::cout << " T: \t " << I0 << " \t -> \t " << I << " \t time "<< s_time << std::endl;
		I0=I;
		stationary_dist[I][0] += 1.0;
		tot_cnt+=1.0;
	}

	if( v_set[5] || v_set[4] ) {
		// NOTE THAT THE FIRST FRAME IS SENT TO PDB FOR TOPOLOGY
		write_atoms_to_xtc( cell , &mmdb , run_args.second[5] , run_args.second[4] , vimnr , 0 );
	}

	double X1 = 0.0;
	double X2 = 0.0;
	double tot_sum = 0.0;
	for( int i=0 ; i<nModels ; i++ ) {
		double val = stationary_dist[i][0]/tot_cnt;
		X1 += val;
		X2 += val*val;
		stationary_dist[i][1] = val;
		tot_sum+=stationary_dist[i][1];
	//	std::cout << " @ " << i << " " << val << std::endl;
	}
	double fModels = nModels;
	std::cout	<< "INFO:: IDEAL:: " << 1/fModels << " MEAN " << X1/fModels 
			<< " STD " << sqrt((X2-X1*X1/fModels)/fModels) << " MSTATS" << std::endl;
	std::cout	<< "INFO:: SUM:: " << tot_sum << std::endl;
//
	PPCAtom 	cur_atoms;
	int		nAtoms	= 0;
	int		nAtoms0	= 0;

	int selHnd0 = mmdb.NewSelection();
	mmdb.Select (	selHnd0    , STYPE_ATOM , 1 , "*",
			ANY_RES    , "*" , ANY_RES,"*" ,
			"*", "*"   , "*" , "*", SKEY_NEW );
	mmdb.Select (	selHnd0    , STYPE_ATOM , 1 , "*",
			ANY_RES    , "*" , ANY_RES ,"*" ,
			"HOH", "*" , "*" , "*", SKEY_CLR );
	mmdb.GetSelIndex ( selHnd0 , cur_atoms, nAtoms );
	std::cout << "INFO::ZEROTH " << nAtoms << std::endl;

	gsl_matrix *XM1 = gsl_matrix_calloc(nAtoms,3);
	gsl_matrix *XM2 = gsl_matrix_calloc(nAtoms,3);
	gsl_vector *rx1 = gsl_vector_calloc(3);
	gsl_vector *rx2 = gsl_vector_calloc(3);

	for( int imodel	= 1; imodel <= nModels ; imodel++ ) {

		int selHnd = mmdb.NewSelection();
		mmdb.Select (	selHnd  , STYPE_ATOM , imodel , "*",
				ANY_RES , "*" , ANY_RES,"*" ,
				"*", "*", "*" , "*", SKEY_NEW );
		mmdb.Select (	selHnd  , STYPE_ATOM , imodel , "*",
				ANY_RES , "*" , ANY_RES ,"*" ,
				"HOH", "*", "*" , "*", SKEY_CLR );
		mmdb.GetSelIndex ( selHnd,cur_atoms, nAtoms );

		double wi = stationary_dist[imodel][1];
		for(int iatom=0;iatom<nAtoms;iatom++) {
			gsl_matrix_get_row( rx1, XM1, iatom );
			gsl_matrix_get_row( rx2, XM2, iatom );
			double x1 = gsl_vector_get(rx1,0);
			double y1 = gsl_vector_get(rx1,1);
			double z1 = gsl_vector_get(rx1,2);
			x1+=cur_atoms[iatom]->x*wi;
			y1+=cur_atoms[iatom]->y*wi;
			z1+=cur_atoms[iatom]->z*wi;
			double x2 = gsl_vector_get(rx2,0);
			double y2 = gsl_vector_get(rx2,1);
			double z2 = gsl_vector_get(rx2,2);
			x2+=cur_atoms[iatom]->x*cur_atoms[iatom]->x*wi;
			y2+=cur_atoms[iatom]->y*cur_atoms[iatom]->y*wi;
			z2+=cur_atoms[iatom]->z*cur_atoms[iatom]->z*wi;
			gsl_vector_set(rx1,0,x1);
			gsl_vector_set(rx1,1,y1);
			gsl_vector_set(rx1,2,z1);
			gsl_vector_set(rx2,0,x2);
			gsl_vector_set(rx2,1,y2);
			gsl_vector_set(rx2,2,z2);
			
			gsl_matrix_set_row( XM1, iatom, rx1 );
			gsl_matrix_set_row( XM2, iatom, rx2 );
		}

		if(imodel==1)
			nAtoms0 = nAtoms;

		if(nAtoms0-nAtoms!=0)
			std::cout << "INFO::WARNING::BAD_ATOM_NUMBERS_IN_MODELS" << std::endl;
		else
			if(verbose)
				std::cout << "INFO::" << imodel << " " << nAtoms << std::endl;
	}

	if( v_set[1] ) {
	// FINALLY BUILD THE END MODEL
	// WRITTEN IN A REDUNDANT WAY
	// STAT OUTP IO
		PPCModel	model_tmp;
		int Nmods;
		mmdb.GetModelTable( model_tmp , Nmods );
		CModel	*model0 = new CModel;
		model0	-> Copy (  model_tmp[0]	);	// ZEROTH ONE USED FOR TMP STORAGE
		mmdb	. DeleteAllModels(	);	
		mmdb	. AddModel ( model0 	);

		int selHnd = mmdb.NewSelection();
		mmdb.Select (	selHnd  , STYPE_ATOM , 1 , "*",
				ANY_RES , "*" , ANY_RES,"*" ,	
				"*", "*", "*" , "*", SKEY_NEW );
		mmdb.Select (	selHnd  , STYPE_ATOM , 1 , "*",
				ANY_RES , "*" , ANY_RES ,"*" ,
				"HOH", "*", "*" , "*", SKEY_CLR );
		mmdb.GetSelIndex ( selHnd, cur_atoms, nAtoms );

		for(int iatom=0;iatom<nAtoms;iatom++) {
			double xm = gsl_matrix_get(XM1,iatom,0);
			double ym = gsl_matrix_get(XM1,iatom,1);
			double zm = gsl_matrix_get(XM1,iatom,2);
	
			double xs = gsl_matrix_get(XM2,iatom,0);
			double ys = gsl_matrix_get(XM2,iatom,1);
			double zs = gsl_matrix_get(XM2,iatom,2);

			cur_atoms[iatom]->x=xm;
			cur_atoms[iatom]->y=ym;
			cur_atoms[iatom]->z=zm;
//	looks odd
			double B = (xs-xm*xm)+(ys-ym*ym)+(zs-zm*zm);
//			double B = (xs-xm)*(xs-xm)+(ys-ym)*(ys-ym)+(zs-zm)*(zs-zm); // this is wrong... but I forgot why...

		//std::cout << "INFO::B:: \t "<< iatom << "   " << B << std::endl;
			cur_atoms[iatom]->tempFactor=B;
		}
		std::cout << "INFO:: WILL WRITE FILE " << run_args.second[1] << std::endl;
		int bytes = run_args.second[1].length();
		char * cpstr	= new char[bytes + 1];
		std::strcpy (cpstr, run_args.second[1].c_str() );
		int rval	= mmdb.WritePDBASCII( cpstr );	//	READ PDB
		if(rval!=0) {
			fatal();
		}
	}

	gsl_matrix_free( XM1 );
	gsl_matrix_free( XM2 );
	gsl_vector_free( rx1 );
	gsl_vector_free( rx2 );
	gsl_matrix_free(  M  );
	gsl_matrix_free(  W  );
	gsl_matrix_free(  P  );
	gsl_matrix_free(  S  );
	gsl_vector_free( csr );
	gsl_vector_free( row );
	gsl_vector_free( col );

	return 0;
}
