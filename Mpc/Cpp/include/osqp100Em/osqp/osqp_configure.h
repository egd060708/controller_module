#ifndef OSQP_CONFIGURE_H
# define OSQP_CONFIGURE_H


/* OSQP_ENABLE_DEBUG */
/* #undef OSQP_ENABLE_DEBUG */

/* Operating system */
/* #undef IS_LINUX */
/* #undef IS_MAC */
#define IS_WINDOWS

/* Algebra backend in use - Only one of the following is defined at compilation time */
#define OSQP_ALGEBRA_BUILTIN
/* #undef OSQP_ALGEBRA_MKL */
/* #undef OSQP_ALGEBRA_CUDA */

/* Enable code generation */
/* #undef OSQP_CODEGEN */

/* Enable profiler annotations */
/* #undef OSQP_PROFILER_ANNOTATIONS */

/* Enable derivative computation in the solver */
#define OSQP_ENABLE_DERIVATIVES

/* OSQP_EMBEDDED_MODE */
#define OSQP_EMBEDDED_MODE (2)

/* Header file containing custom memory allocators */
/* #undef OSQP_CUSTOM_MEMORY */

/* OSQP_ENABLE_PRINTING */
/* #undef OSQP_ENABLE_PRINTING */

/* Header file containing custom printing functions */
/* #undef OSQP_CUSTOM_PRINTING */

/* OSQP_ENABLE_PROFILING */
/* #undef OSQP_ENABLE_PROFILING */

/* OSQP_ENABLE_INTERRUPT */
/* #undef OSQP_ENABLE_INTERRUPT */

/* OSQP_USE_FLOAT */
/* #undef OSQP_USE_FLOAT */

/* OSQP_USE_LONG */
/* #undef OSQP_USE_LONG */

/* OSQP_PACK_SETTINGS (advanced usage) */
/* #undef OSQP_PACK_SETTINGS */

#endif /* ifndef OSQP_CONFIGURE_H */
