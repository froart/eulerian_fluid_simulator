#include "fluid_dynamics.hpp"
#include "opengl_setup.hpp"
#include <math.h>
#include <iostream>
#include <omp.h>
#include <memory>
#include <algorithm>
#include <vector>

using namespace std;

#define I(x,y) ((x)+(w_)*(y))
#define SWAP_ARRAYS(a,b) {float* tmp=a;a=b;b=tmp;}

Fluid::Fluid(float* image,
						 int width, 
						 int height, 
						 float dt,
						 int iterations) 
						 : d_(image),
						   w_(width), 	
							 h_(height),
							 dt_(dt),
							 it_(iterations) {
	cell_num_ = (width + 2) * (height + 2);
	fill_n(v_.begin(), cell_num_, 0.0);	// x-velocity grid 
	fill_n(u_.begin(), cell_num_, 0.0); // y-velocity grid
	fill_n(v1_.begin(), cell_num_, 0.0); // new x-velocity grid
	fill_n(u1_.begin(), cell_num_, 0.0); // new y-velocity grid
	fill_n(p_.begin(), cell_num_, 0.0); // pressure grid
	fill_n(s_.begin(), cell_num_, 0.0); // scalers grid: 0.0 for wall, 1 for fluid 
	fill_n(m_.begin(), cell_num_, 0.0); // TODO what is it for?
};

void Fluid::addSmoke(int x, int y, float amount) {
	d_[I(x,y)] += amount;
}

void Fluid::addWind(int x, int y, float x_amount, float y_amount) {
	v_[I(x,y)] += x_amount;
	u_[I(x,y)] += y_amount;
}

void Fluid::evaluate() {
	for(int i = 0; i < cell_num_; i++)
		p_[i] = 0.0; 
	project();
	extrapolate();
	advect_velocity();
	advect_smoke();
}

void Fluid::project() {
}

void Fluid::extrapolate() {
}

void Fluid::advect_velocity() {
}

void Fluid::advect_smoke() {
}

Fluid::~Fluid() {
	delete[] d_;
}
