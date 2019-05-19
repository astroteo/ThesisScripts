#ifndef __c17_test_hill_sim_h__
#define __c17_test_hill_sim_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc17_test_hill_simInstanceStruct
#define typedef_SFc17_test_hill_simInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c17_sfEvent;
  boolean_T c17_isStable;
  boolean_T c17_doneDoubleBufferReInit;
  uint8_T c17_is_active_c17_test_hill_sim;
  real_T *c17_a_sat;
  real_T (*c17_a_c)[3];
  real_T (*c17_a_c_SAT)[3];
} SFc17_test_hill_simInstanceStruct;

#endif                                 /*typedef_SFc17_test_hill_simInstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c17_test_hill_sim_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c17_test_hill_sim_get_check_sum(mxArray *plhs[]);
extern void c17_test_hill_sim_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
