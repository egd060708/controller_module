/*
 * This file is part of QDLDL, a library for performing the LDL^T factorization
 * of a symmetric indefinite matrix.
 *
 * QDLDL is part of the OSQP project, and is available at https://github.com/osqp/qdldl.
 *
 * Copyright 2018, Paul Goulart, Bartolomeo Stellato, Goran Banjac, Ian McInerney, The OSQP developers
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * SPDX-License-Identifier: Apache-2.0
 * SPDX-ExternalRef: PACKAGE_MANAGER purl pkg:github/osqp/qdldl
 */
#ifndef QDLDL_TYPES_H
#define QDLDL_TYPES_H

#ifdef __cplusplus
extern "C" {
#endif /* ifdef __cplusplus */

// For the QDLDL_INT_TYPE_MAX
#include <limits.h>

// QDLDL integer and float types

typedef int    QDLDL_int;   /* for indices */
typedef double  QDLDL_float; /* for numerical values  */
typedef unsigned char   QDLDL_bool;  /* for boolean values  */

//Maximum value of the signed type QDLDL_int.
#define QDLDL_INT_MAX INT_MAX

/*
 * Configuration options
 */

/* When defined, QDLDL is using floats instead of doubles */
/* #undef QDLDL_FLOAT */

/* When defined, QDLDL is using long long instead of int types */
/* #undef QDLDL_LONG */

#ifdef __cplusplus
}
#endif /* ifdef __cplusplus */

#endif /* ifndef QDLDL_TYPES_H */
