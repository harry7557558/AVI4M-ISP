Use octatree to triangulate implicit surfaces. Directly build hierarchy structure that can be used for accelerating raytracing and store it in a linear `int16` buffer.

### Parameters:

`vec3 P0, P1`: the corners of the grid

`ivec3 SEARCH_DIF`: outermost grid size

`int PLOT_DEPTH`: depth of the hierarchy

`int EDGE_ROUNDING`: determines the precision of coordinates

All vertices in the mesh have coordinates in `[0, (SEARCH_DIF<<PLOT_DEPTH)*EDGE_ROUNDING]`.

### Buffer structure specification:

Designed to be placed in GLSL `highp uint` array.

Started with `prod(SEARCH_DIF)` pointers, the start of the top layer (grid), in flattened `[z][y][x]` order;

A block in a middle layer contains `8` pointers, the children in the next layer;

A block in the bottom layer contains an integer `n`, the number of triangles, followed by `3*n` pointers to their vertices;

A list of vertices, each contains 3 integers, the XYZ coordinates;

Color?? The above subjects to change.

All integers are 16-bit unsigned. All pointers are 32-bit little endian. A null pointer is represented by `0x0000 0x0000`.

```glsl
const ivec3 VERTEX_LIST[8] = {
    ivec3(0,0,0), ivec3(0,1,0), ivec3(1,1,0), ivec3(1,0,0),
    ivec3(0,0,1), ivec3(0,1,1), ivec3(1,1,1), ivec3(1,0,1)
};
```
