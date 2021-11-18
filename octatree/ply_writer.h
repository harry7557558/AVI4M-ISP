#pragma once

#include <cstdio>

#include "meshtype.h"


bool WritePLY(const char* filename,
	const vec3 *vertices, int vn, const ivec3 *faces, int fn,
	const vec3 *vertice_cols = nullptr) {

	FILE* fp = fopen(filename, "wb");
	if (!fp) return false;

	// ply header
	fprintf(fp, "ply\n");
	fprintf(fp, "format binary_little_endian 1.0\n");
	fprintf(fp, "element vertex %d\n", vn);
	fprintf(fp, "property float x\nproperty float y\nproperty float z\n");
	if (vertice_cols)
		fprintf(fp, "property uchar red\nproperty uchar green\nproperty uchar blue\n");
	fprintf(fp, "element face %d\n", fn);
	fprintf(fp, "property list uchar int vertex_indices\n");
	fprintf(fp, "end_header\n");

	// vertices
	for (int i = 0; i < vn; i++) {
		fwrite((float*)&vertices[i], 4, 3, fp);
		if (vertice_cols) {
			vec3 c = saturate(vertice_cols[i]);
			uint8_t b[3] = {
				(uint8_t)(255.99*c.x),
				(uint8_t)(255.99*c.y),
				(uint8_t)(255.99*c.z)
			};
			fwrite(b, 1, 3, fp);
		}
	}

	// faces
	for (int i = 0; i < fn; i++) {
		fputc(3, fp);
		fwrite((int*)&faces[i], 4, 3, fp);
	}

	return fclose(fp) == 0;
}
