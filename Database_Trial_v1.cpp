// FILE: model_prop_table.cpp
//
// PROGRAMMERS: Elisa Antolini, Jeremy Heyl
//
// SEDRIS C++ SDK Release 4.1.4 - July 8, 2013

/*
 *  This program show how to create a new transmittal. It contains only a <Transmittal Root>, 
    a model library with ten models.

    Each star model has a Feature Model and a Feature Model Instance.
    Each feature model instance has a property table with a single
    regular axis running from time zero to 20,000,000,000 years.  The
    table also has three Table Property Descriptions (one for mass,
    one for radius and one for temperature).  Finally under each Table
    Property Description there are ten Property Values that have been
    set to reasonable but fake values (with the correct units).

    The number of the current star model (i.e. star_0 through star_9) and 
    the number of the table description (i.e. mass, radius or temperature 0, 1 or 2) are
    stored in a static variable, so that the lower levels of the tree
    can access these data.

 *  This can be used as a workspace for building object and adding
 *  objects from other transmittals.
 *
 *  Created by elisa antolini on 08/07/13.
 *  Updated by jeremy heyl on 09/07/13.
 *  Updated by jeremy heyl on 12/07/13.
 *  Copyright 2013 elyan. All rights reserved.
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>

#include "seWorkspace.h"
#include "seTransmittal.h"
#include "seDRMAll.h"
#include "seHelperDataTable.h"
#include "seDRMEnvironmentRoot.h"
#include "seIterator.h"


//added for loadfile function
#define ALLOCBLOCK 512
#define MAXCOLUMNS 100
#define MAXLINESIZE 5192
const char * const LOADFILE_FS=" \t";

using namespace sedris;

const char * const SedrisVersionString = "4.1.x";
const char * const ToolName = "STF Database";
const char * const ToolVersionString = "4.1.4.0";
const char * const BaseToolString = "%s v%s\n    (compatible with SEDRIS C++ SDK %s)\n\n\n";

static const char *OUT_FNAME = "Database_Trial_v1.stf"; //Output file_name


double *data[5]; //double *data[NSTARMODELS][NCOLUMNS][Insert with Time] allocare la memoria with malloc 
//and fill with for (i=0;i<NSTAR;i++) data[i]=loadfile(…)
//nines = loadfile("",data[i]); data[i][5] = (double *) malloc(nlines*sizeof(double));

unsigned int col[]={1,2,34,35,40,41}; //Evolution Parameter of the STAR in function of time history.data
unsigned int col2[]={1,3,4,5,6,18,22,23,24,25}; //Evolution Parameter of the STAR in function of the header of the profile100.data
unsigned int col3[]={1,3,6,26,29,36,62,67,71,72,73,74,75}; //Evolution Parameter of the STAR in function of radius in profile100.data

static int  nlines;
static int nModels = 34; // Total numbers of models
//int mesa_lines;
double *Ages[1];
seObject star_model;
seObject *iterator = new seObject[10000];
double SM= 1.98855e30; //1 Solar Mass [Kg]
double SR= 6.960e8; //1 Solar Radius [m]
double SL= 3.827e26; //1 Solar Luminosity [Watt]
double SD= 150 ; //1 Solar density in the center [gr/cm3]
char model_name[2000],files_name[2000];
char All_Models[2000][40];

double minT;
int minrow;




//load catalogue

unsigned int
loadfilestr_fileptr_fs(FILE *in, unsigned int ncolumns, unsigned int columns[], char **data[], const char *fs) {
	unsigned int i, j, ialloc;
	char buffer[MAXLINESIZE];
	char **ap, *argv[MAXCOLUMNS], *inputstring;
	
	for (j=0;j<ncolumns;j++) {
		assert ((data[j]=(char **) malloc(sizeof(char *)*512))!=NULL);
	}
	ialloc=512;
	
	i=0;
	while (fgets(buffer,MAXLINESIZE-1,in)) {
		if (buffer[0]!='#'){
			inputstring=buffer;
            
			/* break line into up to MAXCOLUMNS columns */
			for (ap = argv; (*ap = strsep(&inputstring, fs)) != NULL;)
				if (**ap != '\0')
					if (++ap >= &argv[MAXCOLUMNS])
						break;
			/* were there any tokens? */
			if (ap>argv) {
				/* assign columns to the data arrays; missing values given nan */
				for (j=0;j<ncolumns;j++) {
					if (argv+columns[j]<=ap) {
						assert( (data[j][i]=(char *) malloc((strlen(argv[columns[j]-1])+1)*sizeof(char)))!=NULL);
						strcpy(data[j][i],argv[columns[j]-1]);
					} else {
						assert( (data[j][i]=(char *) malloc(sizeof(char)))!=NULL);
						strcpy(data[j][i],"");
					}
				}
				i++;
				/* do we need to allocate more memory? */
				if (i==ialloc) {
					for (j=0;j<ncolumns;j++) {
						assert ((data[j]=(char **) realloc((void *) data[j],sizeof(char *)*(i+512)))!=NULL);
					}
					ialloc=i+512;
				}
			}
		}
	}
	/* free excess memory */
	for (j=0;j<ncolumns;j++) {
		assert ((data[j]=(char **) realloc((void *) data[j],sizeof(char *)*i))!=NULL);
	}
	
	/* return the number of lines read */
	return i;
}

unsigned int
loadfilestr_fs(const char *filename, unsigned int ncolumns, unsigned int columns[], char **data[], const char *fs) {
	FILE *in;
	unsigned int retval;
	
	assert((in=fopen(filename,"r"))!=NULL);
	
	retval=loadfilestr_fileptr_fs(in, ncolumns, columns, data, fs);
	fclose(in);
	return retval;
}

unsigned int
loadfilestr_fileptr(FILE *in, unsigned int ncolumns, unsigned int columns[], char **data[]) {
	return loadfilestr_fileptr_fs(in,ncolumns,columns,data,LOADFILE_FS);
}

unsigned int
loadfilestr(const char *filename, unsigned int ncolumns, unsigned int columns[], char **data[]) {
	return loadfilestr_fs(filename,ncolumns,columns,data,LOADFILE_FS);
}


	
//loadfile code

unsigned int
loadfile_fileptr_fs(FILE *in, unsigned int ncolumns, unsigned int columns[], double *data[], const char *fs, int skipFirstRow) {
	unsigned int i, j, ialloc;
	char buffer[MAXLINESIZE];
	char **ap, *argv[MAXCOLUMNS], *inputstring;
	
	for (j=0;j<ncolumns;j++) {
		assert ((data[j]=(double *) malloc(sizeof(double)*512))!=NULL);
	}
	ialloc=512;
	
	i=0;
	while (fgets(buffer,MAXLINESIZE-1,in)) {
        
        
		if (buffer[0]!='#') {
			inputstring=buffer;
            
			/* break line into up to MAXCOLUMNS columns */
			for (ap = argv; (*ap = strsep(&inputstring, fs)) != NULL;)
				if (**ap != '\0')
					if (++ap >= &argv[MAXCOLUMNS])
						break;
			/* were there any tokens? */
			if (ap>argv) {
				/* assign columns to the data arrays; missing values given nan */
				for (j=0;j<ncolumns;j++) {
                    if (skipFirstRow == 0){
                        data[j][i]=(argv+columns[j]<=ap ? atof(argv[columns[j]-1]) : 0.0/0.0);
                    }
                    if(skipFirstRow == 1 and i > 1){
                        data[j][i]=(argv+columns[j]<=ap ? atof(argv[columns[j]-1]) : 0.0/0.0);
                        //printf("%f \n",data[j][i]);
                        
                    }
				}
				i++;
				/* do we need to allocate more memory? */
				if (i==ialloc) {
					for (j=0;j<ncolumns;j++) {
						assert ((data[j]=(double *) realloc((void *) data[j],sizeof(double)*(i+512)))!=NULL);
					}
					ialloc=i+512;
				}
			}
		}
        
        
	}
	/* free excess memory */
	for (j=0;j<ncolumns;j++) {
		assert ((data[j]=(double *) realloc((void *) data[j],sizeof(double)*i))!=NULL);
	}
	
	/* return the number of lines read */
	return i;
}
 
unsigned int
loadfile_fs(const char *filename, unsigned int ncolumns, unsigned int columns[], double *data[], const char *fs,int skipFirstRow) {
	FILE *in;
	unsigned int retval;
	
    
	assert((in=fopen(filename,"r"))!=NULL);
	
	retval=loadfile_fileptr_fs(in, ncolumns, columns, data, fs, skipFirstRow);
    
	fclose(in);
    
	return retval;
}

unsigned int
loadfile_fileptr(FILE *in, unsigned int ncolumns, unsigned int columns[], double *data[], int skipFirstRow) {
    

	return loadfile_fileptr_fs(in,ncolumns,columns,data,LOADFILE_FS,skipFirstRow);
}

unsigned int
loadfile(const char *filename, unsigned int ncolumns, unsigned int columns[], double *data[], int skipFirstRow) {
   
	return loadfile_fs(filename,ncolumns,columns,data,LOADFILE_FS,skipFirstRow);
}


//Create a Relation between STAR models and <Feature Model Instance>å

seObject create_relation(seObject &rel, int num)
{
	
	int j=num;
	seObject k = rel;
	iterator[j]=k;
	//k.print();
	return k;
}

/*
 *-----------------------------------------------------------------------------
 *
 * FUNCTION: add_next_level
 *
 * PARAMETERS:
 *
 *   sc -- storage collection to which new objects are being added
 *
 *   parent -- object whose components are being generated by this function
 *
 *   level -- the current level within Model Library's aggregation tree
 *
 * RETURNS:
 *
 *   nothing
 *
 *-----------------------------------------------------------------------------
 */


void add_next_level
(
 seTransmittal &transmittal,
 seObject &parent,
 int level,
 int loc_index,
 int n_lines,
 int n_lines2,
 int n_lines3,
 //double *fill_data[6],
 //unsigned int num_lines[3],
 double *fill_data[6],
 double *fill_data2[10],
 double *fill_data3[13],
 char Model_name[100]
 )
{
	/* these two variables store which star model and which table we are working on */
	/* so the lower levels can use this information to complete their property values */
	static int model_i, propdescription_i;
	int num_comp_no_link_this_level;
    
	seDRMPropertyTable ptObj;
	SE_DRM_Class level_type[] =
    {
		SE_CLS_DRM_MODEL,
		SE_CLS_DRM_FEATURE_MODEL,
		SE_CLS_DRM_FEATURE_MODEL_INSTANCE,
		SE_CLS_DRM_PROPERTY_TABLE,
		SE_CLS_DRM_TABLE_PROPERTY_DESCRIPTION,
		SE_CLS_DRM_PROPERTY_VALUE
        
    };
    
	
    //n_lines = num_lines[0];
    
    if ( level >= (sizeof(level_type)/sizeof(SE_DRM_Class)) )
        return;
    else if ( level == 0 )
        num_comp_no_link_this_level = 1; //number of models
    else if ( level == (sizeof(level_type)/sizeof(SE_DRM_Class) - 1) )
		num_comp_no_link_this_level = n_lines;
    else if  ( level_type[level] == SE_CLS_DRM_TABLE_PROPERTY_DESCRIPTION )
	/* one table each for mass, radius,temperature, Color Index (B_V) and Absolute Magnitude (V_M) */
		num_comp_no_link_this_level = 5;
	else if ( level_type[level] == SE_CLS_DRM_PROPERTY_VALUE)
		num_comp_no_link_this_level = n_lines; //number of property values = nlines in file
	else  num_comp_no_link_this_level = 1;
	
    for ( int i = 0; i < num_comp_no_link_this_level; i++ )
    {
        seObject new_obj;
		
        transmittal.createObject(new_obj, level_type[level]);
        
        if (level_type[level] == SE_CLS_DRM_TABLE_PROPERTY_DESCRIPTION){
            printf("%d\n",i);
        }
		
		
        if ( level_type[level] == SE_CLS_DRM_MODEL )
        {
			//char model_name[100], file_name[200];
			seDRMModel model(new_obj);
		
			
			// store the model number for the lower levels
			
			model_i = i;
			
			
			//sprintf(model_name, "STAR_%dMesa", i+1);
			model.set_name(Model_name);
			/*sprintf(file_name,"/Volumes/Macintosh HD/sedris_cpp_sdk_4.1.4/src/apps/final_database_v2/%s",model_name);
			printf("%s\n",file_name);
			n_lines=loadfile(file_name,4,col,fill_data);  //return the lines number of the raw in each column
			printf("%d\n",n_lines);*/
			
			
			/*for (int j=0;j<n_lines;j++) {
			 //printf("|%g|%g|%g|%g|\n",fill_data[0][j]*1e-03,fill_data[1][j]*SM*1e-03*1e-24,pow(10,fill_data[2][j])*SR*1e-09,fill_data[3][j]*1e-09);				
			 }*/
			
		}
     
	//Create a link between <Feature Model> (<Model Library>) and a <Feature Model Instance> (<Environment Root>)
		
		if ( level_type[level] == SE_CLS_DRM_FEATURE_MODEL )
		{
			
			star_model = new_obj;
			create_relation(star_model,model_i);
			//printf("%d\n",model_i);
					
		}
	
	
		
		if ( level_type[level] == SE_CLS_DRM_PROPERTY_TABLE ) 
		{
            //Add irregular Axis
			seDRMIrregularAxis irraxisObj;
			SE_Short_Integer_Positive vc;
            
			
			// Here I just made an array of sv
			// SE_Single_Value sv[10];
			
			// But you won't know how may values you have at compile time, so this is better
			//int nlines=10;  // you will know nlines from the loadfile command and this will probably be a global variable
			SE_Single_Value *sv=new SE_Single_Value[n_lines+1];
			
			// time axis first
			// here we set up a regular time axis
			// we probably will want to use seDRMIrregularAxis
			// which requires us to list the time points
			// here it might make sense to use the model_i variable because
			// each type of star might have a different set of time points.
			transmittal.createObject(irraxisObj);
			irraxisObj.set_axis_type(EAC_TIME_COORDINATE);
			irraxisObj.set_value_unit(EUC_YEAR);
			irraxisObj.set_value_scale(ESC_MEGA);//1e6
			irraxisObj.set_interpolation_type(SE_INTERPTYP_DIAGONALIZATION);
			irraxisObj.set_axis_value_count(n_lines);
			
			
			for (vc=0;vc<n_lines+1;vc++) {
				sv[vc].value.long_float_value = fill_data[0][vc]*1e-06; //fill up with time
				sv[vc].value_type = SE_SINGVALTYP_LONG_FLOAT;
				irraxisObj.set_axis_value_array(sv,vc);
			}
			
			//Get Lines from Mesa
			
			//getMesaLines(n_lines);
			
			/*for (vc=0;vc<n_lines;vc++) {
			getAge(fill_data[3][vc]*1e-09,vc,n_lines); //Get the Ages from MESA
			}*/
			
			
			new_obj.addComponent(irraxisObj);
            
            
            ///////// PROPERTY TABLE OF THE PROFILE ///////////
            
            seDRMPropertyTable ProfileTable;
            transmittal.createObject(ProfileTable);
            new_obj.addComponent(ProfileTable);
            
            
            
            SE_DRM_Class level_type2 = ProfileTable.getClassType();
            seDRMClassificationData ClassData;

            transmittal.createObject(ClassData);
            ClassData.set_tag(ECC_PROPERTY_SET);
                
            //ProfileTable.addComponent(ClassData);
            level_type2 = ClassData.getClassType();
            
        
            
            if(level_type2 == SE_CLS_DRM_CLASSIFICATION_DATA){
            
                
                
                seDRMPropertyValue pv2(ClassData);
                
                transmittal.createObject(pv2);
                
                SE_Property_Code elemType2;
                EDCS_Attribute_Value av2;
                
                
                elemType2.code_type = SE_PROPCODTYP_VARIABLE;
                elemType2.code.variable = SE_VARCOD_REF_VEC_V0;
                pv2.set_meaning(elemType2);
                
                // I pulled this working code from stf_test.cpp
                // I think that it is important that all of these things are set.
                av2.attribute_value_type = EDCS_AVT_REAL;
                //av.value.real_value.unit_scale = ESC_UNI;
                av2.value.real_value.numeric_value_type = EDCS_NVT_SINGLE_VALUE;
                

                av2.value.real_value.unit = EUC_KILOGRAM;
                av2.value.real_value.unit_scale = ESC_YOTTA;//e24
                av2.value.real_value.value.single_value = fill_data[1][1]*SM*1e-24;//YOTTA Kg
            
                pv2.set_value(av2);

                }
            
               ProfileTable.addComponent(ClassData);
            
                //TABLE PROPERTY DESCRIPTIONS OF THE PROFILES
                //seDRMTablePropertyDescription propDescObj2(new_obj);
                
                //ProfileTable.addComponent(propDescObj2);

                
            }
            
            
		//}
        
        
		
		if (level_type[level] ==  SE_CLS_DRM_TABLE_PROPERTY_DESCRIPTION )
		{
			
			propdescription_i=i;
			seDRMTablePropertyDescription propDescObj(new_obj);
			SE_Element_Type elemType;
			elemType.code_type = SE_ELEMTYPCOD_ATTRIBUTE;
            
			if (i==0) {
				elemType.code.attribute = EAC_MASS; //Solar Mass to be defined
			}
            
			if (i==1) {
				elemType.code.attribute = EAC_TEMPERATURE;
			}
            
			if (i==2) {
				elemType.code.attribute = EAC_POWER_GENERATION_RATE; // Solar Luminosity to be defined -> EAC_POWER Is not there
			}
            
			if (i==3) {
				elemType.code.attribute = EAC_AIR_DENSITY; //Density -> EAC_VOLUMIC_MASS is not there
			}
			if (i==4) {
				elemType.code.attribute = EAC_PRESSURE_ALTIMETER_SETTING;// Pressure -> EAC_PRESSURE is not there
			}
            
            //if (i==5) {
				//elemType.code.attribute = EAC_OBJECT_ASPECT;// To Be Defined as EAC_MAGNITUDE
			//}
				
			propDescObj.set_meaning(elemType);
			
			/* is the table for mass in units of mass (gram), radius (metres) or temperature (Kelvin)? */
			// store the type of table for the lower levels
		
			
			if (i==0) {
				propDescObj.set_value_unit(EUC_KILOGRAM);
				propDescObj.set_value_scale(ESC_YOTTA);//e24
			}
            
			if (i==1) {
				propDescObj.set_value_unit(EUC_KELVIN);
				propDescObj.set_value_scale(ESC_KILO);//e3
                
			}
            
			if (i==2){
				propDescObj.set_value_unit(EUC_WATT);
				propDescObj.set_value_scale(ESC_YOTTA);//e24
			}
			if (i==3){
		
				propDescObj.set_value_unit(EUC_GRAM_PER_CUBIC_CM); // Density
				propDescObj.set_value_scale(ESC_MEGA);//
			}
			
			if(i==4){
                propDescObj.set_value_unit(EUC_PASCAL);
				propDescObj.set_value_scale(ESC_UNI);//e9

			}
            
            //else{
			//	propDescObj.set_value_unit(EUC_UNITLESS);//Absolute Magnitude
			//	propDescObj.set_value_scale(ESC_UNI);//
			//}
			
			propDescObj.set_value_type(SE_DTDATAVALTYP_SINGLE_INTEGER);

		}
            
        
		
		if ( level_type[level] == SE_CLS_DRM_PROPERTY_VALUE) 
		{
			seDRMPropertyValue pv(new_obj);
			SE_Property_Code elemType;
			EDCS_Attribute_Value av;
			
			
			elemType.code_type = SE_PROPCODTYP_VARIABLE;
			elemType.code.variable = SE_VARCOD_REF_VEC_V0;
			pv.set_meaning(elemType);
			
			// I pulled this working code from stf_test.cpp
			// I think that it is important that all of these things are set.
			av.attribute_value_type = EDCS_AVT_REAL;
			//av.value.real_value.unit_scale = ESC_UNI;
			av.value.real_value.numeric_value_type = EDCS_NVT_SINGLE_VALUE;
			
			if (propdescription_i==0) {         // table of masses
				av.value.real_value.unit = EUC_KILOGRAM;
				av.value.real_value.unit_scale = ESC_YOTTA;//e24
				av.value.real_value.value.single_value = fill_data[1][i]*SM*1e-24;//YOTTA Kg
				
			}
            
			if (propdescription_i==1) {  // table of temperature
                av.value.real_value.unit = EUC_KELVIN;
				av.value.real_value.unit_scale = ESC_KILO;//e3
				av.value.real_value.value.single_value = pow(10,fill_data[2][i])*1e-03; //Kilo Kelvin
				
			}
            
			if (propdescription_i==2) {  // table of Luminosity
                av.value.real_value.unit = EUC_WATT;
				av.value.real_value.unit_scale = ESC_YOTTA;//e24
				av.value.real_value.value.single_value = fill_data[3][i]*SL*1e-24;// Yotta Watt
				
			}
			if (propdescription_i==3) {  // table of Density in the Center
				av.value.real_value.unit = EUC_GRAM_PER_CUBIC_CM;
				av.value.real_value.unit_scale = ESC_MEGA;
				av.value.real_value.value.single_value = pow(10,fill_data[4][i])*1e-06*SD; // Mega gr/cm3
			}
			if (propdescription_i==4) {  // table of Pressure in the Center
				av.value.real_value.unit = EUC_PASCAL;
				av.value.real_value.unit_scale = ESC_UNI;
				av.value.real_value.value.single_value = fill_data[5][i];
			}
            
            //else if (propdescription_i==5) {  // table of Magnitudes (M_V)
				//av.value.real_value.unit = EUC_UNITLESS;
				//av.value.real_value.unit_scale = ESC_UNI;
				//av.value.real_value.value.single_value = fill_data[5][i];
			//}*/
             
			pv.set_value(av);
		}
		
        parent.addComponent(new_obj);
		
        add_next_level(transmittal, new_obj, level+1, i, n_lines,n_lines2,n_lines3,fill_data,fill_data2,fill_data3, Model_name);
		
		if (level==0) {
			/*Free memory for Evolution Parameter*/ 
			free( (void *) data[0]);
			free( (void *) data[1]);
			free( (void *) data[2]);
			free( (void *) data[3]);
		}
    }
	
}



/*
 *-----------------------------------------------------------------------------
 *
 * FUNCTION: add_env_root
 *
 * PARAMETERS:
 *
 *   sc -- storage collection to which new objects are being added
 *
 *   parent -- object whose components are being generated by this function
 *
 *   level -- the current level within Model Library's aggregation tree
 *
 * RETURNS:
 *
 *   nothing
 *
 *-----------------------------------------------------------------------------
 */



void add_env_root
(
 seTransmittal &transmittal,
 seObject &parent,
 int level,
 int loc_index,
 int nlines_str,
 char **fill_data_str[8],
 int nlines_dbl,
 double *fill_data_dbl[8]
 )
{
	
	
	
	static int propervalue_i, properfeature_i;
	int num_comp_no_link_this_level;

	SE_DRM_Class level_type[] =
	{
		SE_CLS_DRM_ENVIRONMENT_ROOT,
		SE_CLS_DRM_FEATURE_MODEL_INSTANCE,
		SE_CLS_DRM_WORLD_TRANSFORMATION,
		SE_CLS_DRM_EI_3D_LOCATION
		
	};
		
	if ( level >= (sizeof(level_type)/sizeof(SE_DRM_Class)) )
        return;
    if ( level == 0 )
        num_comp_no_link_this_level = 1; //Number of <Environment Root>
	else if ( level == 1) 
		num_comp_no_link_this_level = nlines_dbl; //Numbers of <Feature Model Instance>
	else  num_comp_no_link_this_level = 1; //Number of all the other components
	
    for ( int i = 0; i < num_comp_no_link_this_level; i++ )
    {
        seObject new_obj;
		
        transmittal.createObject(new_obj, level_type[level]);
		
		if ( level_type[level] == SE_CLS_DRM_ENVIRONMENT_ROOT )
        {
		
			seDRMEnvironmentRoot my_envRoot(new_obj);
			SRM_SRF_Parameters_Info my_srf_parameters_info;

			
			//set up SRF (Spatial Reference frame)
			
			my_srf_parameters_info.srf_params_info_code = SRM_SRFPARAMINFCOD_TEMPLATE;
			my_srf_parameters_info.value.srf_template.template_code = SRM_SRFTCOD_EQUATORIAL_INERTIAL;
			my_srf_parameters_info.value.srf_template.orm_code = SRM_ORMCOD_EARTH_INERTIAL_J2000r0;
			
			my_envRoot.set_srf_context_info_srf_parameters_info(my_srf_parameters_info);
			
			my_envRoot.set_srf_context_info_angular_unit(EUC_DEGREE_ARC);
			my_envRoot.set_srf_context_info_linear_unit(EUC_PARSEC);
			my_envRoot.set_srf_context_info_linear_scale(ESC_UNI);
			
		}
		
		//Create 10  <Property Value> for each <Feature Model Instance> (8 String and 2 Real)
		
		if ( level_type[level] == SE_CLS_DRM_FEATURE_MODEL_INSTANCE)
		{ 
			
			properfeature_i= i; //defines which value fills the Property Value
			//printf("%d\n",i);
			seObject make_rel;
			make_rel = iterator[0]; //I loaded Just one model for now 
			make_rel.print();
			new_obj.addAssociate(make_rel);
			
		}
		
		
		if ( level_type[level] == SE_CLS_DRM_WORLD_TRANSFORMATION)
		{ 
				
			
			for(int j=0; j<10; j++)
			{
			 
				seObject new_pvobj;
				transmittal.createObject(new_pvobj, SE_CLS_DRM_PROPERTY_VALUE);
				seDRMPropertyValue pv(new_pvobj);
				propervalue_i= j; //defines two string and two double Property value
				
				SE_Property_Code elemType;
				EDCS_Attribute_Value av;
			
				int length = 0;
				
				//Define Star ID as String (column 1)
				
				if(propervalue_i == 0){
					
					length = sizeof fill_data_str[0][properfeature_i];
					
					elemType.code_type = SE_PROPCODTYP_ATTRIBUTE; //Defined as Attribute
					elemType.code.attribute = EAC_OBJECT_LABEL; //Contains a String
					pv.set_meaning(elemType);
					av.attribute_value_type = EDCS_AVT_STRING;//Contains a String
					av.value.string_value.length= length;
					av.value.string_value.characters= fill_data_str[0][properfeature_i]; 
					length = 0;
				}
				
				//Define HIP number as String (column 2)
				
				if(propervalue_i == 1){
					
					length = sizeof fill_data_str[1][properfeature_i];
					
					elemType.code_type = SE_PROPCODTYP_ATTRIBUTE; //Defined as Attribute
					elemType.code.attribute = EAC_OBJECT_LABEL; //Contains a String
					pv.set_meaning(elemType);
					av.attribute_value_type = EDCS_AVT_STRING;//Contains a String
					av.value.string_value.length= length;
					av.value.string_value.characters= fill_data_str[1][properfeature_i];
					length = 0;
				}
				
				//Define HD number as String (column 3)
				
				if(propervalue_i == 2){
					
					length = sizeof fill_data_str[2][properfeature_i];
					
					elemType.code_type = SE_PROPCODTYP_ATTRIBUTE; //Defined as Attribute
					elemType.code.attribute = EAC_OBJECT_LABEL; //Contains a String
					pv.set_meaning(elemType);
					av.attribute_value_type = EDCS_AVT_STRING;//Contains a String
					av.value.string_value.length= length;
					av.value.string_value.characters= fill_data_str[2][properfeature_i]; 
					length = 0;
					
				}
				
				//Define HR number as String (column 4)
				
				if(propervalue_i == 3){
					
					length = sizeof fill_data_str[3][properfeature_i];
					
					elemType.code_type = SE_PROPCODTYP_ATTRIBUTE; //Defined as Attribute
					elemType.code.attribute = EAC_OBJECT_LABEL; //Contains a String
					pv.set_meaning(elemType);
					av.attribute_value_type = EDCS_AVT_STRING;//Contains a String
					av.value.string_value.length= length;
					av.value.string_value.characters= fill_data_str[3][properfeature_i]; 
					length = 0;
					
				}
				
				//Define Gliese number as String (column 5)
				
				if(propervalue_i == 4){
					
					length = sizeof fill_data_str[4][properfeature_i];
					
					elemType.code_type = SE_PROPCODTYP_ATTRIBUTE; //Defined as Attribute
					elemType.code.attribute = EAC_OBJECT_LABEL; //Contains a String
					pv.set_meaning(elemType);
					av.attribute_value_type = EDCS_AVT_STRING;//Contains a String
					av.value.string_value.length= length;
					av.value.string_value.characters= fill_data_str[4][properfeature_i]; 
					length = 0;
					
				}
				
				//Define BayerFlamsteed Designation as String (column 6)
				
				if(propervalue_i == 5){
					
					length = sizeof fill_data_str[5][properfeature_i];
					
					elemType.code_type = SE_PROPCODTYP_ATTRIBUTE; //Defined as Attribute
					elemType.code.attribute = EAC_OBJECT_LABEL; //Contains a String
					pv.set_meaning(elemType);
					av.attribute_value_type = EDCS_AVT_STRING;//Contains a String
					av.value.string_value.length= length;
					av.value.string_value.characters= fill_data_str[5][properfeature_i]; 
					length = 0;
					
				}
			
				//Define ProperName String (column 7)
			
				if(propervalue_i == 6){
					
					length = sizeof fill_data_str[6][properfeature_i];
				
					elemType.code_type = SE_PROPCODTYP_ATTRIBUTE; //Defined as Attribute
					elemType.code.attribute = EAC_OBJECT_LABEL; //Contains a String
					pv.set_meaning(elemType);
					av.attribute_value_type = EDCS_AVT_STRING;//Contains a String
					av.value.string_value.length= length;
					av.value.string_value.characters= fill_data_str[6][properfeature_i]; 
					length = 0;
				}
			
			
				//Define SpectralType as String (column 16)
			
				else if(propervalue_i == 7){
				
					length = sizeof fill_data_str[7][properfeature_i];
				
					elemType.code_type = SE_PROPCODTYP_ATTRIBUTE; //Defined as Attribute
					elemType.code.attribute = EAC_CLASSIFICATION_NAME; //Contains a String
					pv.set_meaning(elemType);
					av.attribute_value_type = EDCS_AVT_STRING;//Contains a String
					av.value.string_value.length= length;
					av.value.string_value.characters= fill_data_str[7][properfeature_i];
					length = 0;
				}
			
			
				//Define AbsMag as Double (column 15)
			
				else if(propervalue_i == 8){
				
					elemType.code_type = SE_PROPCODTYP_ATTRIBUTE; //Defined as Attribute
					elemType.code.attribute = EAC_MEAN_LUMINOUS_INTENSITY;
					pv.set_meaning(elemType);
					av.attribute_value_type = EDCS_AVT_REAL;//Contains a number
					av.value.real_value.numeric_value_type = EDCS_NVT_SINGLE_VALUE;
					av.value.real_value.unit = EUC_UNITLESS;
					av.value.real_value.unit_scale = ESC_UNI;
					av.value.real_value.value.single_value = fill_data_dbl[3][properfeature_i];
				}
			
				//Define ColorIndex as Double (column 17)
			
				else if(propervalue_i == 9){
				
					elemType.code_type = SE_PROPCODTYP_VARIABLE; //Defined as Variable
					elemType.code.variable = SE_VARCOD_CLR_INDEX;
					pv.set_meaning(elemType);
					av.attribute_value_type = EDCS_AVT_REAL;//Contains a number
					av.value.real_value.numeric_value_type = EDCS_NVT_SINGLE_VALUE;
					av.value.real_value.unit = EUC_UNITLESS;
					av.value.real_value.unit_scale = ESC_UNI;
					av.value.real_value.value.single_value = fill_data_dbl[4][properfeature_i];
				}
			
				pv.set_value(av);
			
			
				parent.addComponent(new_pvobj);
			}
			
		}
		
		//set up Coordinate Location of each star
		
		if ( level_type[level] == SE_CLS_DRM_EI_3D_LOCATION)
		{ 
			
			seDRMEI3DLocation my_coord(new_obj);
			
			SRM_Long_Float ra = fill_data_dbl[0][properfeature_i]*15; //convert RA from hours to degree
			SRM_Long_Float dec = fill_data_dbl[1][properfeature_i];
			SRM_Long_Float radius = fill_data_dbl[2][properfeature_i];
			
			my_coord.set_coordinate_right_ascension(ra);
			my_coord.set_coordinate_declination(dec);
			my_coord.set_coordinate_radius(radius);
			
		}
		
			parent.addComponent(new_obj);
		
			add_env_root(transmittal, new_obj, level+1, loc_index, nlines_str, fill_data_str, nlines_dbl, fill_data_dbl);
	}
  
	
}



/*
 *-----------------------------------------------------------------------------
 *-----------------------------------------------------------------------------
 */
int main( int argc, char **argv )
{
    if ( argc > 1){
	
	}
    try
    {
        seWorkspace             wksp; //include .h
        seTransmittal           transmittal; //include .h
        seDRMTransmittalRoot    root_obj;
        seDRMModelLibrary       my_mlib;

		
		
		
		printf(BaseToolString, ToolName, ToolVersionString, SedrisVersionString);
        printf("Creating \"%s\" transmittal...\n\n", OUT_FNAME);
        wksp.createTransmittal(OUT_FNAME, transmittal);
		
        transmittal.createObject(root_obj);
        root_obj.set_name(OUT_FNAME);
		transmittal.setRootObject(root_obj);
		
	
	/* add a model library */
		
        transmittal.createObject(my_mlib);
        root_obj.addComponent(my_mlib);
		
		printf("Model Library added!\n");
		
		
		
    /*Load  MESA models*/
		
        char *file_name, argom[200];
        char *file_name2, argom2[200];
        char *file_name3, argom3[200];
        char write_filename[200];
		unsigned int nmesa_data[nModels];
        unsigned int nProfHead_data[nModels];
		
        unsigned int startmesa_data[nModels];
        

        
		double *mesa_data[nModels][7];
        double *ProfHead_data[nModels][10];
        int num,h;
        FILE *f_plot;
        int k=0;
        
        file_name = "history.data";
        file_name2 = "profile100_original.data";
        file_name3 = "profile100.data";
       
        
        
        for (k=0;k<=2;k++) {
            
            if (k==0){
               // printf("%d\n",k);
                
               nmesa_data[k]=loadfile(file_name,6,col,mesa_data[k],0);  //return the lines number of the raw in each column
               printf("Number of Raws = %d\n",nmesa_data[k]);
               printf("second column first row \n");
               printf("%f\n",mesa_data[0][1][0]);
            }
            if (k==1){
                //printf("%d\n",k);
                //Origianl profile file -> file_name3
                nmesa_data[k]=loadfile(file_name2,10,col2,mesa_data[k],0);  //return the lines number of the raw in each column
                nmesa_data[k]=1; // Take just the first row
                printf("Second File \n");
                printf("Number of Raws = %d\n",nmesa_data[k]);
                printf("second column first row \n");
                printf("%f\n",mesa_data[1][1][0]);
            }
            
            //Read all the file opening it and put inside the file_name3 the pointer of the file skipping the first row
            
            if (k==2){
                //printf("%d\n",k);
	      FILE *in;
	      int skip;
	      char buffer[MAXLINESIZE];
	      if (1) {
		printf("Hello\n");
		// skip the first four lines
		in=fopen(file_name2,"r");
		for (skip=0;skip<7;skip++) {
		  fgets(buffer,MAXLINESIZE-1,in);
		  printf("%s",buffer);
		}
	        nmesa_data[k]=loadfile_fileptr(in,13,col3,mesa_data[k],0);  //return the lines number of the raw in each column
		fclose(in);
	      } else {
		nmesa_data[k]=loadfile(file_name2,13,col3,mesa_data[k],1);  //return the lines number of the raw in each column

	      }
	      //                Nmesa_data[k]= nmesa_data[k] - 1;
                printf("Third File \n");
                printf("Number of Raws = %d\n",nmesa_data[k]);
                printf("second column first row \n");
                printf("%f\n",mesa_data[2][1][0]);
            }
        }
   
        
        
        //Associate the SEDRIS ATTRIBUTE TO THE MESA FILE FIRST
        
        
        
        //Teff = mesa_data[k][MESA_Teff]; // Temperature
        
        //LogL_Lsun = mesa_data[k][MESA_LogL]; // Luminosity
        
        /*Add MESA models in Model Library*/
        
        add_next_level(transmittal, my_mlib, 0, 0, nmesa_data[0],nmesa_data[1],nmesa_data[2],mesa_data[0],mesa_data[1],mesa_data[2],"Model1");
        
       
		
		
		
		/* add Environment Root */	

		//add_env_root(transmittal, root_obj, 0, 0, 0, 0, nread_dbl, data_dbl);
        //add_env_root(transmittal, root_obj, 0, 0, 0, 0, nread_dbl, data_dbl);
		
		printf("Environment Root added!\n");
		
	
	

		/* free all of the individual char arrays */
		
		/*for (k=0; k<nModels; k++) {
			
			for (j=0;j<7;j++) {
			
				free((void *) mesa_data[k][j]);
			}
			
		}*/
		
	transmittal.close();
		
        printf("Done!\n");
		
	}
	catch ( seException &e )
    {
        fprintf(stderr, "Error - %s\n", e.getWhat());
    }
    return 0;
} /* end main */
