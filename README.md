# Welcome to MMath!
MMath is a single-header math library made in C for use with OpenGL and other similar graphics APIs.

###### Quick refernces:
- [Features](#features)
- [On the to-do list](#on-the-to-do-list)
- [How to use MMath](#how-to-use-mmath)

###### Additional information:

MMath has been lightly tested with both MSVC and GCC, but problems may still exist. Please report any bugs found in the [issues section](https://github.com/8bitslime/MMath/issues) so they can be dealt with quickly.

Some documentation is provided at the top of [`MMath.h`](./MMath.h), but more will be added in the near future.

---

### Features
- Vectors
- Square matrices
- Quaternions
- Transformations
- Easy appending to:
	- [vectors](https://github.com/8bitslime/MMath/blob/master/MMath.h#L29)
    - [matrices](https://github.com/8bitslime/MMath/blob/master/MMath.h#L48)

---

### On the to-do list
- Variadic functions
- Rectangular matrices
- SIMD optimizations
- Various profilings

---

### How to use MMath
To add MMath to your project, simply put the [`MMath.h`](./MMath.h) header file in your project's include directory and use `#include "MMath.h"` anywhere math is required.

If you require *double* precision, add the line `#define MMATH_DOUBLE` before including [`MMath.h`](./MMath.h).
