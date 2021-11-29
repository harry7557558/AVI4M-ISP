Use octatree to triangulate implicit surfaces. Directly build hierarchy structure that can be used for accelerating raytracing and store it in a linear `uint8` buffer.

### Parameters:

`vec3 P0, P1`: the corners of the grid

`ivec3 SEARCH_DIF`: outermost grid size

`int PLOT_DEPTH`: depth of the hierarchy

`int EDGE_ROUNDING`: between 2 and 255, determines the precision of the coordinates

All vertices in the geometry have absolute coordinates in `[0, (SEARCH_DIF<<PLOT_DEPTH)*EDGE_ROUNDING]` and relative coordinates in `[0, EDGE_ROUNDING]`.

### Buffer structure specification:

Designed to be placed in GLSL `highp uint` array.

Started with `prod(SEARCH_DIF)` 32-bit pointers, the start of the top layer (grid), in flattened `[z][y][x]` order;

A block in a middle layer contains `8` 32-bit pointers, the children in the next layer;

A block in the bottom layer contains a 8-bit integer `n`, the number of triangles, followed by `n` triangles. Each triangle contains `3×3=9` 8-bit integers, the coordinates of the vertices, followed by `3` 8-bit integers, the RGB color.

All integers and pointers are little endian. A null pointer is represented by `0x00000000`.

```glsl
const ivec3 VERTEX_LIST[8] = {
    ivec3(0,0,0), ivec3(0,1,0), ivec3(1,1,0), ivec3(1,0,0),
    ivec3(0,0,1), ivec3(0,1,1), ivec3(1,1,1), ivec3(1,0,1)
};
```
