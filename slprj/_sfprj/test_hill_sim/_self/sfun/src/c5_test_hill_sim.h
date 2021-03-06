#ifndef __c5_test_hill_sim_h__
#define __c5_test_hill_sim_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc5_test_hill_simInstanceStruct
#define typedef_SFc5_test_hill_simInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c5_sfEvent;
  boolean_T c5_isStable;
  boolean_T c5_doneDoubleBufferReInit;
  uint8_T c5_is_active_c5_test_hill_sim;
  real_T (*c5_r_t)[3];
  real_T (*c5_W)[3];
  real_T (*c5_v_t)[3];
} SFc5_test_hill_simInstanceStruct;

#endif                                 /*typedef_SFc5_test_hill_simInstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c5_test_hill_sim_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c5_test_hill_sim_get_check_sum(mxArray *plhs[]);
extern void c5_test_hill_sim_method_dispatcher(SimStruct *S, int_T method, void *
  data);

#endif
