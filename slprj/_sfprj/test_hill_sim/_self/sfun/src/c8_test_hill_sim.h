#ifndef __c8_test_hill_sim_h__
#define __c8_test_hill_sim_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc8_test_hill_simInstanceStruct
#define typedef_SFc8_test_hill_simInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c8_sfEvent;
  boolean_T c8_isStable;
  boolean_T c8_doneDoubleBufferReInit;
  uint8_T c8_is_active_c8_test_hill_sim;
  real_T (*c8_R_Sun_c)[3];
  real_T (*c8_R_c)[3];
  real_T *c8_e_f;
  real_T (*c8_a_Sun)[3];
  real_T *c8_m_cubesat;
} SFc8_test_hill_simInstanceStruct;

#endif                                 /*typedef_SFc8_test_hill_simInstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c8_test_hill_sim_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c8_test_hill_sim_get_check_sum(mxArray *plhs[]);
extern void c8_test_hill_sim_method_dispatcher(SimStruct *S, int_T method, void *
  data);

#endif
