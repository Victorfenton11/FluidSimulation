#version 330

smooth in vec3 fragNorm;	// Interpolated model-space normal
smooth in vec3 fragColor;  // color

out vec3 outCol;	// Final pixel color

void main() {
	// Visualize normals as colors
	outCol = fragColor;  // normalize(fragNorm) * 0.5f + vec3(0.5f);
}
