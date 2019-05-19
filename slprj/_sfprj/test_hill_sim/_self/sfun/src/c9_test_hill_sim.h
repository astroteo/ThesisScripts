#ifndef __c9_test_hill_sim_h__
#define __c9_test_hill_sim_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc9_test_hill_simInstanceStruct
#define typedef_SFc9_test_hill_simInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c9_sfEvent;
  boolean_T c9_isStable;
  boolean_T c9_doneDoubleBufferReInit;
  uint8_T c9_is_active_c9_test_hill_sim;
  real_T (*c9_R_c_Moon)[3];
  real_T (*c9_a_Moon)[3];
} SFc9_test_hill_simInstanceStruct;

#endif                                 /*typedef_SFc9_test_hill_simInstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c9_test_hill_sim_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c9_test_hill_sim_get_check_sum(mxArray *plhs[]);
extern void c9_test_hill_sim_method_dispatcher(SimStruct *S, int_T method, void *
  data);

#endif
