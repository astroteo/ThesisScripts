/* Include files */

#include "test_hill_sim_sfun.h"
#include "test_hill_sim_sfun_debug_macros.h"
#include "c1_test_hill_sim.h"
#include "c2_test_hill_sim.h"
#include "c3_test_hill_sim.h"
#include "c4_test_hill_sim.h"
#include "c5_test_hill_sim.h"
#include "c6_test_hill_sim.h"
#include "c7_test_hill_sim.h"
#include "c8_test_hill_sim.h"
#include "c9_test_hill_sim.h"
#include "c10_test_hill_sim.h"
#include "c11_test_hill_sim.h"
#include "c12_test_hill_sim.h"
#include "c13_test_hill_sim.h"
#include "c14_test_hill_sim.h"
#include "c15_test_hill_sim.h"
#include "c17_test_hill_sim.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */
uint32_T _test_hill_simMachineNumber_;

/* Function Declarations */

/* Function Definitions */
void test_hill_sim_initializer(void)
{
}

void test_hill_sim_terminator(void)
{
}

/* SFunction Glue Code */
unsigned int sf_test_hill_sim_method_dispatcher(SimStruct *simstructPtr,
  unsigned int chartFileNumber, const char* specsCksum, int_T method, void *data)
{
  if (chartFileNumber==1) {
    c1_test_hill_sim_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==2) {
    c2_test_hill_sim_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==3) {
    c3_test_hill_sim_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==4) {
    c4_test_hill_sim_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==5) {
    c5_test_hill_sim_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==6) {
    c6_test_hill_sim_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==7) {
    c7_test_hill_sim_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==8) {
    c8_test_hill_sim_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==9) {
    c9_test_hill_sim_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==10) {
    c10_test_hill_sim_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==11) {
    c11_test_hill_sim_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==12) {
    c12_test_hill_sim_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==13) {
    c13_test_hill_sim_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==14) {
    c14_test_hill_sim_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==15) {
    c15_test_hill_sim_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==17) {
    c17_test_hill_sim_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  return 0;
}

extern void sf_test_hill_sim_uses_exported_functions(int nlhs, mxArray * plhs[],
  int nrhs, const mxArray * prhs[])
{
  plhs[0] = mxCreateLogicalScalar(0);
}

unsigned int sf_test_hill_sim_process_check_sum_call( int nlhs, mxArray * plhs[],
  int nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[20];
  if (nrhs<1 || !mxIsChar(prhs[0]) )
    return 0;

  /* Possible call to get the checksum */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"sf_get_check_sum"))
    return 0;
  plhs[0] = mxCreateDoubleMatrix( 1,4,mxREAL);
  if (nrhs>1 && mxIsChar(prhs[1])) {
    mxGetString(prhs[1], commandName,sizeof(commandName)/sizeof(char));
    commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
    if (!strcmp(commandName,"machine")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2208077481U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2980739902U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(4095585121U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(786556834U);
    } else if (!strcmp(commandName,"exportedFcn")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0U);
    } else if (!strcmp(commandName,"makefile")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3386569575U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(623549423U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(519468450U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(841276536U);
    } else if (nrhs==3 && !strcmp(commandName,"chart")) {
      unsigned int chartFileNumber;
      chartFileNumber = (unsigned int)mxGetScalar(prhs[2]);
      switch (chartFileNumber) {
       case 1:
        {
          extern void sf_c1_test_hill_sim_get_check_sum(mxArray *plhs[]);
          sf_c1_test_hill_sim_get_check_sum(plhs);
          break;
        }

       case 2:
        {
          extern void sf_c2_test_hill_sim_get_check_sum(mxArray *plhs[]);
          sf_c2_test_hill_sim_get_check_sum(plhs);
          break;
        }

       case 3:
        {
          extern void sf_c3_test_hill_sim_get_check_sum(mxArray *plhs[]);
          sf_c3_test_hill_sim_get_check_sum(plhs);
          break;
        }

       case 4:
        {
          extern void sf_c4_test_hill_sim_get_check_sum(mxArray *plhs[]);
          sf_c4_test_hill_sim_get_check_sum(plhs);
          break;
        }

       case 5:
        {
          extern void sf_c5_test_hill_sim_get_check_sum(mxArray *plhs[]);
          sf_c5_test_hill_sim_get_check_sum(plhs);
          break;
        }

       case 6:
        {
          extern void sf_c6_test_hill_sim_get_check_sum(mxArray *plhs[]);
          sf_c6_test_hill_sim_get_check_sum(plhs);
          break;
        }

       case 7:
        {
          extern void sf_c7_test_hill_sim_get_check_sum(mxArray *plhs[]);
          sf_c7_test_hill_sim_get_check_sum(plhs);
          break;
        }

       case 8:
        {
          extern void sf_c8_test_hill_sim_get_check_sum(mxArray *plhs[]);
          sf_c8_test_hill_sim_get_check_sum(plhs);
          break;
        }

       case 9:
        {
          extern void sf_c9_test_hill_sim_get_check_sum(mxArray *plhs[]);
          sf_c9_test_hill_sim_get_check_sum(plhs);
          break;
        }

       case 10:
        {
          extern void sf_c10_test_hill_sim_get_check_sum(mxArray *plhs[]);
          sf_c10_test_hill_sim_get_check_sum(plhs);
          break;
        }

       case 11:
        {
          extern void sf_c11_test_hill_sim_get_check_sum(mxArray *plhs[]);
          sf_c11_test_hill_sim_get_check_sum(plhs);
          break;
        }

       case 12:
        {
          extern void sf_c12_test_hill_sim_get_check_sum(mxArray *plhs[]);
          sf_c12_test_hill_sim_get_check_sum(plhs);
          break;
        }

       case 13:
        {
          extern void sf_c13_test_hill_sim_get_check_sum(mxArray *plhs[]);
          sf_c13_test_hill_sim_get_check_sum(plhs);
          break;
        }

       case 14:
        {
          extern void sf_c14_test_hill_sim_get_check_sum(mxArray *plhs[]);
          sf_c14_test_hill_sim_get_check_sum(plhs);
          break;
        }

       case 15:
        {
          extern void sf_c15_test_hill_sim_get_check_sum(mxArray *plhs[]);
          sf_c15_test_hill_sim_get_check_sum(plhs);
          break;
        }

       case 17:
        {
          extern void sf_c17_test_hill_sim_get_check_sum(mxArray *plhs[]);
          sf_c17_test_hill_sim_get_check_sum(plhs);
          break;
        }

       default:
        ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0.0);
      }
    } else if (!strcmp(commandName,"target")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2515920432U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3908508645U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2530489944U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2924353608U);
    } else {
      return 0;
    }
  } else {
    ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(4158895491U);
    ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3433748098U);
    ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1554242913U);
    ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1336271974U);
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_test_hill_sim_autoinheritance_info( int nlhs, mxArray * plhs[],
  int nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[32];
  char aiChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]) )
    return 0;

  /* Possible call to get the autoinheritance_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_autoinheritance_info"))
    return 0;
  mxGetString(prhs[2], aiChksum,sizeof(aiChksum)/sizeof(char));
  aiChksum[(sizeof(aiChksum)/sizeof(char)-1)] = '\0';

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        if (strcmp(aiChksum, "Dt9t43ZibJg9Q830ilaPcC") == 0) {
          extern mxArray *sf_c1_test_hill_sim_get_autoinheritance_info(void);
          plhs[0] = sf_c1_test_hill_sim_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 2:
      {
        if (strcmp(aiChksum, "dSxiyxhoC9B9ZPxT5r1FL") == 0) {
          extern mxArray *sf_c2_test_hill_sim_get_autoinheritance_info(void);
          plhs[0] = sf_c2_test_hill_sim_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 3:
      {
        if (strcmp(aiChksum, "qeIpyoXjADvIjVAro9fseG") == 0) {
          extern mxArray *sf_c3_test_hill_sim_get_autoinheritance_info(void);
          plhs[0] = sf_c3_test_hill_sim_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 4:
      {
        if (strcmp(aiChksum, "kBgpzCkSK6MwExxHFcg5RD") == 0) {
          extern mxArray *sf_c4_test_hill_sim_get_autoinheritance_info(void);
          plhs[0] = sf_c4_test_hill_sim_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 5:
      {
        if (strcmp(aiChksum, "e3MVNLCwSeBYQf3Q1bYj2") == 0) {
          extern mxArray *sf_c5_test_hill_sim_get_autoinheritance_info(void);
          plhs[0] = sf_c5_test_hill_sim_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 6:
      {
        if (strcmp(aiChksum, "e3MVNLCwSeBYQf3Q1bYj2") == 0) {
          extern mxArray *sf_c6_test_hill_sim_get_autoinheritance_info(void);
          plhs[0] = sf_c6_test_hill_sim_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 7:
      {
        if (strcmp(aiChksum, "huGOsAXfP87hikxBDcfshF") == 0) {
          extern mxArray *sf_c7_test_hill_sim_get_autoinheritance_info(void);
          plhs[0] = sf_c7_test_hill_sim_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 8:
      {
        if (strcmp(aiChksum, "Mj5HKuTYlGwOfi1Ziet7uB") == 0) {
          extern mxArray *sf_c8_test_hill_sim_get_autoinheritance_info(void);
          plhs[0] = sf_c8_test_hill_sim_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 9:
      {
        if (strcmp(aiChksum, "6FGv3FCeGqDujgJNkeZ0zD") == 0) {
          extern mxArray *sf_c9_test_hill_sim_get_autoinheritance_info(void);
          plhs[0] = sf_c9_test_hill_sim_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 10:
      {
        if (strcmp(aiChksum, "Mt0vRI09jaUJFBOgHk0LEB") == 0) {
          extern mxArray *sf_c10_test_hill_sim_get_autoinheritance_info(void);
          plhs[0] = sf_c10_test_hill_sim_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 11:
      {
        if (strcmp(aiChksum, "rC43iWV0DM0cQYMHgO1uhH") == 0) {
          extern mxArray *sf_c11_test_hill_sim_get_autoinheritance_info(void);
          plhs[0] = sf_c11_test_hill_sim_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 12:
      {
        if (strcmp(aiChksum, "zGRCC5HaYAKYKT69O5JXh") == 0) {
          extern mxArray *sf_c12_test_hill_sim_get_autoinheritance_info(void);
          plhs[0] = sf_c12_test_hill_sim_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 13:
      {
        if (strcmp(aiChksum, "YSdLoE9uyGdH1Q5rvyQYMB") == 0) {
          extern mxArray *sf_c13_test_hill_sim_get_autoinheritance_info(void);
          plhs[0] = sf_c13_test_hill_sim_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 14:
      {
        if (strcmp(aiChksum, "zlgIWMJPGkwMiRCJSvNnTB") == 0) {
          extern mxArray *sf_c14_test_hill_sim_get_autoinheritance_info(void);
          plhs[0] = sf_c14_test_hill_sim_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 15:
      {
        if (strcmp(aiChksum, "WG6ut1EsSACocMxs7LMPlE") == 0) {
          extern mxArray *sf_c15_test_hill_sim_get_autoinheritance_info(void);
          plhs[0] = sf_c15_test_hill_sim_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 17:
      {
        if (strcmp(aiChksum, "7JHkTr2Jfy800CjMIZ8MuG") == 0) {
          extern mxArray *sf_c17_test_hill_sim_get_autoinheritance_info(void);
          plhs[0] = sf_c17_test_hill_sim_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_test_hill_sim_get_eml_resolved_functions_info( int nlhs, mxArray
  * plhs[], int nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[64];
  if (nrhs<2 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the get_eml_resolved_functions_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_eml_resolved_functions_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        extern const mxArray
          *sf_c1_test_hill_sim_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c1_test_hill_sim_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 2:
      {
        extern const mxArray
          *sf_c2_test_hill_sim_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c2_test_hill_sim_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 3:
      {
        extern const mxArray
          *sf_c3_test_hill_sim_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c3_test_hill_sim_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 4:
      {
        extern const mxArray
          *sf_c4_test_hill_sim_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c4_test_hill_sim_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 5:
      {
        extern const mxArray
          *sf_c5_test_hill_sim_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c5_test_hill_sim_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 6:
      {
        extern const mxArray
          *sf_c6_test_hill_sim_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c6_test_hill_sim_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 7:
      {
        extern const mxArray
          *sf_c7_test_hill_sim_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c7_test_hill_sim_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 8:
      {
        extern const mxArray
          *sf_c8_test_hill_sim_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c8_test_hill_sim_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 9:
      {
        extern const mxArray
          *sf_c9_test_hill_sim_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c9_test_hill_sim_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 10:
      {
        extern const mxArray
          *sf_c10_test_hill_sim_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c10_test_hill_sim_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 11:
      {
        extern const mxArray
          *sf_c11_test_hill_sim_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c11_test_hill_sim_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 12:
      {
        extern const mxArray
          *sf_c12_test_hill_sim_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c12_test_hill_sim_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 13:
      {
        extern const mxArray
          *sf_c13_test_hill_sim_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c13_test_hill_sim_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 14:
      {
        extern const mxArray
          *sf_c14_test_hill_sim_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c14_test_hill_sim_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 15:
      {
        extern const mxArray
          *sf_c15_test_hill_sim_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c15_test_hill_sim_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 17:
      {
        extern const mxArray
          *sf_c17_test_hill_sim_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c17_test_hill_sim_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_test_hill_sim_third_party_uses_info( int nlhs, mxArray * plhs[],
  int nrhs, const mxArray * prhs[] )
{
  char commandName[64];
  char tpChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the third_party_uses_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  mxGetString(prhs[2], tpChksum,sizeof(tpChksum)/sizeof(char));
  tpChksum[(sizeof(tpChksum)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_third_party_uses_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        if (strcmp(tpChksum, "p0whV3tFDQKf8CVuy7Nh6D") == 0) {
          extern mxArray *sf_c1_test_hill_sim_third_party_uses_info(void);
          plhs[0] = sf_c1_test_hill_sim_third_party_uses_info();
          break;
        }
      }

     case 2:
      {
        if (strcmp(tpChksum, "ZUpMNVlLLORUYOllxMRgfF") == 0) {
          extern mxArray *sf_c2_test_hill_sim_third_party_uses_info(void);
          plhs[0] = sf_c2_test_hill_sim_third_party_uses_info();
          break;
        }
      }

     case 3:
      {
        if (strcmp(tpChksum, "wN4GAzBP5nj64wAY1Cd8aD") == 0) {
          extern mxArray *sf_c3_test_hill_sim_third_party_uses_info(void);
          plhs[0] = sf_c3_test_hill_sim_third_party_uses_info();
          break;
        }
      }

     case 4:
      {
        if (strcmp(tpChksum, "alWxm0YUNIMSWe7ib54K7D") == 0) {
          extern mxArray *sf_c4_test_hill_sim_third_party_uses_info(void);
          plhs[0] = sf_c4_test_hill_sim_third_party_uses_info();
          break;
        }
      }

     case 5:
      {
        if (strcmp(tpChksum, "VNIU2KUAjpdfiQCWOeG2F") == 0) {
          extern mxArray *sf_c5_test_hill_sim_third_party_uses_info(void);
          plhs[0] = sf_c5_test_hill_sim_third_party_uses_info();
          break;
        }
      }

     case 6:
      {
        if (strcmp(tpChksum, "VNIU2KUAjpdfiQCWOeG2F") == 0) {
          extern mxArray *sf_c6_test_hill_sim_third_party_uses_info(void);
          plhs[0] = sf_c6_test_hill_sim_third_party_uses_info();
          break;
        }
      }

     case 7:
      {
        if (strcmp(tpChksum, "LN2KjxY7hXuK1UIfjsGmcC") == 0) {
          extern mxArray *sf_c7_test_hill_sim_third_party_uses_info(void);
          plhs[0] = sf_c7_test_hill_sim_third_party_uses_info();
          break;
        }
      }

     case 8:
      {
        if (strcmp(tpChksum, "sVXdLiO3cmBabCb4yoSmzG") == 0) {
          extern mxArray *sf_c8_test_hill_sim_third_party_uses_info(void);
          plhs[0] = sf_c8_test_hill_sim_third_party_uses_info();
          break;
        }
      }

     case 9:
      {
        if (strcmp(tpChksum, "8X2ew74QSX1h6e4TuAPnhG") == 0) {
          extern mxArray *sf_c9_test_hill_sim_third_party_uses_info(void);
          plhs[0] = sf_c9_test_hill_sim_third_party_uses_info();
          break;
        }
      }

     case 10:
      {
        if (strcmp(tpChksum, "aAOlwhXBZFCmkmZXjgtmXE") == 0) {
          extern mxArray *sf_c10_test_hill_sim_third_party_uses_info(void);
          plhs[0] = sf_c10_test_hill_sim_third_party_uses_info();
          break;
        }
      }

     case 11:
      {
        if (strcmp(tpChksum, "eQzREDZpwWCURyH4DOUCmG") == 0) {
          extern mxArray *sf_c11_test_hill_sim_third_party_uses_info(void);
          plhs[0] = sf_c11_test_hill_sim_third_party_uses_info();
          break;
        }
      }

     case 12:
      {
        if (strcmp(tpChksum, "G4ncdIOiWyYdA5wrqTxnbD") == 0) {
          extern mxArray *sf_c12_test_hill_sim_third_party_uses_info(void);
          plhs[0] = sf_c12_test_hill_sim_third_party_uses_info();
          break;
        }
      }

     case 13:
      {
        if (strcmp(tpChksum, "W5KsT1MrGCMRR4kZKqcWCG") == 0) {
          extern mxArray *sf_c13_test_hill_sim_third_party_uses_info(void);
          plhs[0] = sf_c13_test_hill_sim_third_party_uses_info();
          break;
        }
      }

     case 14:
      {
        if (strcmp(tpChksum, "WzTRCLcmOm7qK8LsU71WYG") == 0) {
          extern mxArray *sf_c14_test_hill_sim_third_party_uses_info(void);
          plhs[0] = sf_c14_test_hill_sim_third_party_uses_info();
          break;
        }
      }

     case 15:
      {
        if (strcmp(tpChksum, "VK9KEuOAeSTh4h3tQ9nTIF") == 0) {
          extern mxArray *sf_c15_test_hill_sim_third_party_uses_info(void);
          plhs[0] = sf_c15_test_hill_sim_third_party_uses_info();
          break;
        }
      }

     case 17:
      {
        if (strcmp(tpChksum, "v5J5KGSenGCKU3rBrfLczB") == 0) {
          extern mxArray *sf_c17_test_hill_sim_third_party_uses_info(void);
          plhs[0] = sf_c17_test_hill_sim_third_party_uses_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

unsigned int sf_test_hill_sim_jit_fallback_info( int nlhs, mxArray * plhs[], int
  nrhs, const mxArray * prhs[] )
{
  char commandName[64];
  char tpChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the jit_fallback_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  mxGetString(prhs[2], tpChksum,sizeof(tpChksum)/sizeof(char));
  tpChksum[(sizeof(tpChksum)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_jit_fallback_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        if (strcmp(tpChksum, "p0whV3tFDQKf8CVuy7Nh6D") == 0) {
          extern mxArray *sf_c1_test_hill_sim_jit_fallback_info(void);
          plhs[0] = sf_c1_test_hill_sim_jit_fallback_info();
          break;
        }
      }

     case 2:
      {
        if (strcmp(tpChksum, "ZUpMNVlLLORUYOllxMRgfF") == 0) {
          extern mxArray *sf_c2_test_hill_sim_jit_fallback_info(void);
          plhs[0] = sf_c2_test_hill_sim_jit_fallback_info();
          break;
        }
      }

     case 3:
      {
        if (strcmp(tpChksum, "wN4GAzBP5nj64wAY1Cd8aD") == 0) {
          extern mxArray *sf_c3_test_hill_sim_jit_fallback_info(void);
          plhs[0] = sf_c3_test_hill_sim_jit_fallback_info();
          break;
        }
      }

     case 4:
      {
        if (strcmp(tpChksum, "alWxm0YUNIMSWe7ib54K7D") == 0) {
          extern mxArray *sf_c4_test_hill_sim_jit_fallback_info(void);
          plhs[0] = sf_c4_test_hill_sim_jit_fallback_info();
          break;
        }
      }

     case 5:
      {
        if (strcmp(tpChksum, "VNIU2KUAjpdfiQCWOeG2F") == 0) {
          extern mxArray *sf_c5_test_hill_sim_jit_fallback_info(void);
          plhs[0] = sf_c5_test_hill_sim_jit_fallback_info();
          break;
        }
      }

     case 6:
      {
        if (strcmp(tpChksum, "VNIU2KUAjpdfiQCWOeG2F") == 0) {
          extern mxArray *sf_c6_test_hill_sim_jit_fallback_info(void);
          plhs[0] = sf_c6_test_hill_sim_jit_fallback_info();
          break;
        }
      }

     case 7:
      {
        if (strcmp(tpChksum, "LN2KjxY7hXuK1UIfjsGmcC") == 0) {
          extern mxArray *sf_c7_test_hill_sim_jit_fallback_info(void);
          plhs[0] = sf_c7_test_hill_sim_jit_fallback_info();
          break;
        }
      }

     case 8:
      {
        if (strcmp(tpChksum, "sVXdLiO3cmBabCb4yoSmzG") == 0) {
          extern mxArray *sf_c8_test_hill_sim_jit_fallback_info(void);
          plhs[0] = sf_c8_test_hill_sim_jit_fallback_info();
          break;
        }
      }

     case 9:
      {
        if (strcmp(tpChksum, "8X2ew74QSX1h6e4TuAPnhG") == 0) {
          extern mxArray *sf_c9_test_hill_sim_jit_fallback_info(void);
          plhs[0] = sf_c9_test_hill_sim_jit_fallback_info();
          break;
        }
      }

     case 10:
      {
        if (strcmp(tpChksum, "aAOlwhXBZFCmkmZXjgtmXE") == 0) {
          extern mxArray *sf_c10_test_hill_sim_jit_fallback_info(void);
          plhs[0] = sf_c10_test_hill_sim_jit_fallback_info();
          break;
        }
      }

     case 11:
      {
        if (strcmp(tpChksum, "eQzREDZpwWCURyH4DOUCmG") == 0) {
          extern mxArray *sf_c11_test_hill_sim_jit_fallback_info(void);
          plhs[0] = sf_c11_test_hill_sim_jit_fallback_info();
          break;
        }
      }

     case 12:
      {
        if (strcmp(tpChksum, "G4ncdIOiWyYdA5wrqTxnbD") == 0) {
          extern mxArray *sf_c12_test_hill_sim_jit_fallback_info(void);
          plhs[0] = sf_c12_test_hill_sim_jit_fallback_info();
          break;
        }
      }

     case 13:
      {
        if (strcmp(tpChksum, "W5KsT1MrGCMRR4kZKqcWCG") == 0) {
          extern mxArray *sf_c13_test_hill_sim_jit_fallback_info(void);
          plhs[0] = sf_c13_test_hill_sim_jit_fallback_info();
          break;
        }
      }

     case 14:
      {
        if (strcmp(tpChksum, "WzTRCLcmOm7qK8LsU71WYG") == 0) {
          extern mxArray *sf_c14_test_hill_sim_jit_fallback_info(void);
          plhs[0] = sf_c14_test_hill_sim_jit_fallback_info();
          break;
        }
      }

     case 15:
      {
        if (strcmp(tpChksum, "VK9KEuOAeSTh4h3tQ9nTIF") == 0) {
          extern mxArray *sf_c15_test_hill_sim_jit_fallback_info(void);
          plhs[0] = sf_c15_test_hill_sim_jit_fallback_info();
          break;
        }
      }

     case 17:
      {
        if (strcmp(tpChksum, "v5J5KGSenGCKU3rBrfLczB") == 0) {
          extern mxArray *sf_c17_test_hill_sim_jit_fallback_info(void);
          plhs[0] = sf_c17_test_hill_sim_jit_fallback_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

unsigned int sf_test_hill_sim_updateBuildInfo_args_info( int nlhs, mxArray *
  plhs[], int nrhs, const mxArray * prhs[] )
{
  char commandName[64];
  char tpChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the updateBuildInfo_args_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  mxGetString(prhs[2], tpChksum,sizeof(tpChksum)/sizeof(char));
  tpChksum[(sizeof(tpChksum)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_updateBuildInfo_args_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        if (strcmp(tpChksum, "p0whV3tFDQKf8CVuy7Nh6D") == 0) {
          extern mxArray *sf_c1_test_hill_sim_updateBuildInfo_args_info(void);
          plhs[0] = sf_c1_test_hill_sim_updateBuildInfo_args_info();
          break;
        }
      }

     case 2:
      {
        if (strcmp(tpChksum, "ZUpMNVlLLORUYOllxMRgfF") == 0) {
          extern mxArray *sf_c2_test_hill_sim_updateBuildInfo_args_info(void);
          plhs[0] = sf_c2_test_hill_sim_updateBuildInfo_args_info();
          break;
        }
      }

     case 3:
      {
        if (strcmp(tpChksum, "wN4GAzBP5nj64wAY1Cd8aD") == 0) {
          extern mxArray *sf_c3_test_hill_sim_updateBuildInfo_args_info(void);
          plhs[0] = sf_c3_test_hill_sim_updateBuildInfo_args_info();
          break;
        }
      }

     case 4:
      {
        if (strcmp(tpChksum, "alWxm0YUNIMSWe7ib54K7D") == 0) {
          extern mxArray *sf_c4_test_hill_sim_updateBuildInfo_args_info(void);
          plhs[0] = sf_c4_test_hill_sim_updateBuildInfo_args_info();
          break;
        }
      }

     case 5:
      {
        if (strcmp(tpChksum, "VNIU2KUAjpdfiQCWOeG2F") == 0) {
          extern mxArray *sf_c5_test_hill_sim_updateBuildInfo_args_info(void);
          plhs[0] = sf_c5_test_hill_sim_updateBuildInfo_args_info();
          break;
        }
      }

     case 6:
      {
        if (strcmp(tpChksum, "VNIU2KUAjpdfiQCWOeG2F") == 0) {
          extern mxArray *sf_c6_test_hill_sim_updateBuildInfo_args_info(void);
          plhs[0] = sf_c6_test_hill_sim_updateBuildInfo_args_info();
          break;
        }
      }

     case 7:
      {
        if (strcmp(tpChksum, "LN2KjxY7hXuK1UIfjsGmcC") == 0) {
          extern mxArray *sf_c7_test_hill_sim_updateBuildInfo_args_info(void);
          plhs[0] = sf_c7_test_hill_sim_updateBuildInfo_args_info();
          break;
        }
      }

     case 8:
      {
        if (strcmp(tpChksum, "sVXdLiO3cmBabCb4yoSmzG") == 0) {
          extern mxArray *sf_c8_test_hill_sim_updateBuildInfo_args_info(void);
          plhs[0] = sf_c8_test_hill_sim_updateBuildInfo_args_info();
          break;
        }
      }

     case 9:
      {
        if (strcmp(tpChksum, "8X2ew74QSX1h6e4TuAPnhG") == 0) {
          extern mxArray *sf_c9_test_hill_sim_updateBuildInfo_args_info(void);
          plhs[0] = sf_c9_test_hill_sim_updateBuildInfo_args_info();
          break;
        }
      }

     case 10:
      {
        if (strcmp(tpChksum, "aAOlwhXBZFCmkmZXjgtmXE") == 0) {
          extern mxArray *sf_c10_test_hill_sim_updateBuildInfo_args_info(void);
          plhs[0] = sf_c10_test_hill_sim_updateBuildInfo_args_info();
          break;
        }
      }

     case 11:
      {
        if (strcmp(tpChksum, "eQzREDZpwWCURyH4DOUCmG") == 0) {
          extern mxArray *sf_c11_test_hill_sim_updateBuildInfo_args_info(void);
          plhs[0] = sf_c11_test_hill_sim_updateBuildInfo_args_info();
          break;
        }
      }

     case 12:
      {
        if (strcmp(tpChksum, "G4ncdIOiWyYdA5wrqTxnbD") == 0) {
          extern mxArray *sf_c12_test_hill_sim_updateBuildInfo_args_info(void);
          plhs[0] = sf_c12_test_hill_sim_updateBuildInfo_args_info();
          break;
        }
      }

     case 13:
      {
        if (strcmp(tpChksum, "W5KsT1MrGCMRR4kZKqcWCG") == 0) {
          extern mxArray *sf_c13_test_hill_sim_updateBuildInfo_args_info(void);
          plhs[0] = sf_c13_test_hill_sim_updateBuildInfo_args_info();
          break;
        }
      }

     case 14:
      {
        if (strcmp(tpChksum, "WzTRCLcmOm7qK8LsU71WYG") == 0) {
          extern mxArray *sf_c14_test_hill_sim_updateBuildInfo_args_info(void);
          plhs[0] = sf_c14_test_hill_sim_updateBuildInfo_args_info();
          break;
        }
      }

     case 15:
      {
        if (strcmp(tpChksum, "VK9KEuOAeSTh4h3tQ9nTIF") == 0) {
          extern mxArray *sf_c15_test_hill_sim_updateBuildInfo_args_info(void);
          plhs[0] = sf_c15_test_hill_sim_updateBuildInfo_args_info();
          break;
        }
      }

     case 17:
      {
        if (strcmp(tpChksum, "v5J5KGSenGCKU3rBrfLczB") == 0) {
          extern mxArray *sf_c17_test_hill_sim_updateBuildInfo_args_info(void);
          plhs[0] = sf_c17_test_hill_sim_updateBuildInfo_args_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

void sf_test_hill_sim_get_post_codegen_info( int nlhs, mxArray * plhs[], int
  nrhs, const mxArray * prhs[] )
{
  unsigned int chartFileNumber = (unsigned int) mxGetScalar(prhs[0]);
  char tpChksum[64];
  mxGetString(prhs[1], tpChksum,sizeof(tpChksum)/sizeof(char));
  tpChksum[(sizeof(tpChksum)/sizeof(char)-1)] = '\0';
  switch (chartFileNumber) {
   case 1:
    {
      if (strcmp(tpChksum, "p0whV3tFDQKf8CVuy7Nh6D") == 0) {
        extern mxArray *sf_c1_test_hill_sim_get_post_codegen_info(void);
        plhs[0] = sf_c1_test_hill_sim_get_post_codegen_info();
        return;
      }
    }
    break;

   case 2:
    {
      if (strcmp(tpChksum, "ZUpMNVlLLORUYOllxMRgfF") == 0) {
        extern mxArray *sf_c2_test_hill_sim_get_post_codegen_info(void);
        plhs[0] = sf_c2_test_hill_sim_get_post_codegen_info();
        return;
      }
    }
    break;

   case 3:
    {
      if (strcmp(tpChksum, "wN4GAzBP5nj64wAY1Cd8aD") == 0) {
        extern mxArray *sf_c3_test_hill_sim_get_post_codegen_info(void);
        plhs[0] = sf_c3_test_hill_sim_get_post_codegen_info();
        return;
      }
    }
    break;

   case 4:
    {
      if (strcmp(tpChksum, "alWxm0YUNIMSWe7ib54K7D") == 0) {
        extern mxArray *sf_c4_test_hill_sim_get_post_codegen_info(void);
        plhs[0] = sf_c4_test_hill_sim_get_post_codegen_info();
        return;
      }
    }
    break;

   case 5:
    {
      if (strcmp(tpChksum, "VNIU2KUAjpdfiQCWOeG2F") == 0) {
        extern mxArray *sf_c5_test_hill_sim_get_post_codegen_info(void);
        plhs[0] = sf_c5_test_hill_sim_get_post_codegen_info();
        return;
      }
    }
    break;

   case 6:
    {
      if (strcmp(tpChksum, "VNIU2KUAjpdfiQCWOeG2F") == 0) {
        extern mxArray *sf_c6_test_hill_sim_get_post_codegen_info(void);
        plhs[0] = sf_c6_test_hill_sim_get_post_codegen_info();
        return;
      }
    }
    break;

   case 7:
    {
      if (strcmp(tpChksum, "LN2KjxY7hXuK1UIfjsGmcC") == 0) {
        extern mxArray *sf_c7_test_hill_sim_get_post_codegen_info(void);
        plhs[0] = sf_c7_test_hill_sim_get_post_codegen_info();
        return;
      }
    }
    break;

   case 8:
    {
      if (strcmp(tpChksum, "sVXdLiO3cmBabCb4yoSmzG") == 0) {
        extern mxArray *sf_c8_test_hill_sim_get_post_codegen_info(void);
        plhs[0] = sf_c8_test_hill_sim_get_post_codegen_info();
        return;
      }
    }
    break;

   case 9:
    {
      if (strcmp(tpChksum, "8X2ew74QSX1h6e4TuAPnhG") == 0) {
        extern mxArray *sf_c9_test_hill_sim_get_post_codegen_info(void);
        plhs[0] = sf_c9_test_hill_sim_get_post_codegen_info();
        return;
      }
    }
    break;

   case 10:
    {
      if (strcmp(tpChksum, "aAOlwhXBZFCmkmZXjgtmXE") == 0) {
        extern mxArray *sf_c10_test_hill_sim_get_post_codegen_info(void);
        plhs[0] = sf_c10_test_hill_sim_get_post_codegen_info();
        return;
      }
    }
    break;

   case 11:
    {
      if (strcmp(tpChksum, "eQzREDZpwWCURyH4DOUCmG") == 0) {
        extern mxArray *sf_c11_test_hill_sim_get_post_codegen_info(void);
        plhs[0] = sf_c11_test_hill_sim_get_post_codegen_info();
        return;
      }
    }
    break;

   case 12:
    {
      if (strcmp(tpChksum, "G4ncdIOiWyYdA5wrqTxnbD") == 0) {
        extern mxArray *sf_c12_test_hill_sim_get_post_codegen_info(void);
        plhs[0] = sf_c12_test_hill_sim_get_post_codegen_info();
        return;
      }
    }
    break;

   case 13:
    {
      if (strcmp(tpChksum, "W5KsT1MrGCMRR4kZKqcWCG") == 0) {
        extern mxArray *sf_c13_test_hill_sim_get_post_codegen_info(void);
        plhs[0] = sf_c13_test_hill_sim_get_post_codegen_info();
        return;
      }
    }
    break;

   case 14:
    {
      if (strcmp(tpChksum, "WzTRCLcmOm7qK8LsU71WYG") == 0) {
        extern mxArray *sf_c14_test_hill_sim_get_post_codegen_info(void);
        plhs[0] = sf_c14_test_hill_sim_get_post_codegen_info();
        return;
      }
    }
    break;

   case 15:
    {
      if (strcmp(tpChksum, "VK9KEuOAeSTh4h3tQ9nTIF") == 0) {
        extern mxArray *sf_c15_test_hill_sim_get_post_codegen_info(void);
        plhs[0] = sf_c15_test_hill_sim_get_post_codegen_info();
        return;
      }
    }
    break;

   case 17:
    {
      if (strcmp(tpChksum, "v5J5KGSenGCKU3rBrfLczB") == 0) {
        extern mxArray *sf_c17_test_hill_sim_get_post_codegen_info(void);
        plhs[0] = sf_c17_test_hill_sim_get_post_codegen_info();
        return;
      }
    }
    break;

   default:
    break;
  }

  plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
}

void test_hill_sim_debug_initialize(struct SfDebugInstanceStruct* debugInstance)
{
  _test_hill_simMachineNumber_ = sf_debug_initialize_machine(debugInstance,
    "test_hill_sim","sfun",0,16,0,0,0);
  sf_debug_set_machine_event_thresholds(debugInstance,
    _test_hill_simMachineNumber_,0,0);
  sf_debug_set_machine_data_thresholds(debugInstance,
    _test_hill_simMachineNumber_,0);
}

void test_hill_sim_register_exported_symbols(SimStruct* S)
{
}

static mxArray* sRtwOptimizationInfoStruct= NULL;
mxArray* load_test_hill_sim_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct==NULL) {
    sRtwOptimizationInfoStruct = sf_load_rtw_optimization_info("test_hill_sim",
      "test_hill_sim");
    mexMakeArrayPersistent(sRtwOptimizationInfoStruct);
  }

  return(sRtwOptimizationInfoStruct);
}

void unload_test_hill_sim_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct!=NULL) {
    mxDestroyArray(sRtwOptimizationInfoStruct);
    sRtwOptimizationInfoStruct = NULL;
  }
}
