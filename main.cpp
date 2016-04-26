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

//	SIZE OF D IS 10E-11 to 10E-10 m2/s FOR BIOMOLECULES
//	DIFF EQ:: d/dt(phi) = D d^2/dt^2(phi)	FICKS SECOND
//
//	CONS	VALUE				UNIT
#define KBOL	1.38064889E-23		//	m2 kg s-2 K-1
#define RGAS	8.31446210		//	J K−1 mol−1
#define NAVO	6.022140857E23		//

// COMPILE:: g++ -std=c++11 main.cpp -lmmdb -lgsl -lblas -o rich_dyn

std::pair< int , std::vector< std::string > >
argparser( std::pair < int, char ** > mp, std::vector < std::pair < int, std::pair<std::string, std::string > > > opts_n_defaults ) {

	// ARGS OUTPUT ORDER CORRESPONDS TO OPTS_N_DEFAULTS FIRST VALUE
	std::vector< std::string > args;
	std::pair  < int, std::vector< std::string > > ret_args;
	std::vector< int > v_set;

	for( int i=0 ; i<mp.first ; i++ ) {
		args.push_back			( mp.second[i] );
	}
	for( int i=0 ; i<opts_n_defaults.size() ; i++ ) {
		ret_args.second.push_back	( opts_n_defaults[i].second.second );
		v_set.push_back(0);
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

	return ret_args;
}

void
fatal(void){
	std::string	author("Richard Tjörnhammar (e-mail: richard.tjornhammar@gmail.com)");
	std::cout << "INFO:: PROGRAM FAILED" << std::endl;
	std::cout << "PLEASE CONTACT " << author << "\nWITH ANY QUESTIONS REGARDING THE FUNCTIONALITY" << std::endl;
	exit(1);
}

int main ( int argc, char ** argv ) {
	int verbose = 0;

//	PREPARING ARGPARSER
	std::pair < int, char ** >	margs;
	margs.first = argc; margs.second = argv;

	std::vector < std::pair < int, std::pair<std::string, std::string > > >	opts_n_defaults;
	std::pair	< int, std::pair<std::string,std::string > > 		vipss;
	vipss.first=1;	vipss.second.first = "-ifile"; 
			vipss.second.second = "i_default.pdb";
	opts_n_defaults.push_back(vipss);
	vipss.first=0;	vipss.second.first = "-ofile"; 
			vipss.second.second = "o_default.pdb";
	opts_n_defaults.push_back(vipss);
	vipss.first=0;	vipss.second.first = "-nsteps";
			vipss.second.second = "100";
	opts_n_defaults.push_back(vipss);
	vipss.first=0;	vipss.second.first = "-verbose";
			vipss.second.second = "0";
	opts_n_defaults.push_back(vipss);

//	EXECUTING ARGPARSE
	std::pair<int, std::vector< std::string > > run_args = argparser( margs, opts_n_defaults );
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
	int nErr=0;
	
	gsl_matrix *M = gsl_matrix_calloc(nModels,nModels); // MSD
	gsl_matrix *W = gsl_matrix_calloc(nModels,nModels); // WEIGTHS
	gsl_matrix *P = gsl_matrix_calloc(nModels,nModels); // PROB
	gsl_matrix *S = gsl_matrix_calloc(nModels,nModels); // CUMULATIVE SUM

	double gnrm	= 1.0/sqrt(2.0*M_PI);
	double l_jump	= 0.15E-9;	// JUMP DIFFUSION LENGTH OF PENTANE ( TYPICALLY 1-3 Å )

	double stationary_dist[nModels][2];
	for ( imod=0 ; imod<nModels ; imod++ ) {
		stationary_dist  [imod][0] = 0.0;
		stationary_dist  [imod][1] = 0.0;
		for ( jmod=imod ; jmod<nModels ; jmod++ ) {
			double X = 0.0, X2 = 0.0,Nc = 0.0;
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
	double dice_roll = distribution(generator);			
	int I0 = floor( nModels*dice_roll ); 
	int I  = 0;
	std::cout << std::endl;
	std::cout << "===========================================" << std::endl;
	std::cout << "WILL START TO GENERATE THE TRANSITION ORDER" << std::endl;
	std::cout << "===========================================" << std::endl;
	double tot_cnt=0.0;
	while( s_time < gen_time ) {
		dice_roll = distribution(generator);	
		s_time+=dt;
		gsl_matrix_get_row( row , S , I0 );
		for(int is=0;is<nModels;is++){
			I=is;
			if( gsl_vector_get(row,is) >= dice_roll )
				break;
		}
		if(verbose>=1)
			std::cout << " T: \t " << I0 << " \t -> \t " << I << " \t time "<< s_time << std::endl;
		I0=I;
		stationary_dist[I][0]+=1.0;
		tot_cnt+=1.0;
	}
	double X1 = 0.0;
	double X2 = 0.0;
	double tot_sum=0.0;
	for( int i=0 ; i<nModels ; i++ ) {
		double val = stationary_dist[i][0]/tot_cnt;
		X1 += val;
		X2 += val*val;
		stationary_dist[i][1] = val;
		tot_sum+=stationary_dist[i][1];
		std::cout << " @ " << i << " " << val << std::endl;
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
		double B = (xs-xm*xm)+(ys-ym*ym)+(zs-zm*zm);
		
		std::cout << "INFO::B:: \t "<< iatom << "   " << B << std::endl;
		cur_atoms[iatom]->tempFactor=B;
	}

	if( run_args.second[1].length()>0 ) {
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
