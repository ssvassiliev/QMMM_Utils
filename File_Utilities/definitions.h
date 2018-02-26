typedef struct {
        int nb_points; // nb of points
        float* x; // x-coordinates
        float* y; // y-coordinates
	float* z; // z-coordinates
        short R; // radius
        float c_x; // x-coordinate of center of the sphere
        float c_y; // y-coordinate of center of the sphere
	float c_z; // z-coordinate of center of the sphere
	int* in_histo; // array in which to hold the counts of in waters
	int* out_histo; // array in which to hold the counts of out waters
} SPHERE;

