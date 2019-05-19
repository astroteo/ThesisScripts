#ifndef __c12_test_hill_sim_h__
#define __c12_test_hill_sim_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc12_test_hill_simInstanceStruct
#define typedef_SFc12_test_hill_simInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c12_sfEvent;
  boolean_T c12_isStable;
  boolean_T c12_doneDoubleBufferReInit;
  uint8_T c12_is_active_c12_test_hill_sim;
  real_T *c12_rho;
  real_T (*c12_V)[3];
  real_T (*c12_a_drag)[3];
  real_T *c12_m_cubesat;
} SFc12_test_hill_simInstanceStruct;

#endif                                 /*typedef_SFc12_test_hill_simInstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c12_test_hill_sim_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c12_test_hill_sim_get_check_sum(mxArray *plhs[]);
extern void c12_test_hill_sim_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
