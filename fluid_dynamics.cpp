#include "fluid_dynamics.hpp"
#include "opengl_setup.hpp"
#include <math.h>
#include <iostream>
#include <omp.h>
#include <memory>

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
	int cell_num = (width + 2) * (height + 2);
	unique_ptr<float[]> v_(new float[cell_num]); // x-velocity grid
	unique_ptr<float[]> u_(new float[cell_num]); // y-velocity grid
	unique_ptr<float[]> u1_(new float[cell_num]); // new x-velocity grid
	unique_ptr<float[]> v1_(new float[cell_num]); // new y-velocity grid
	unique_ptr<float[]> p_(new float[cell_num]); // pressure grid
	unique_ptr<float[]> s_(new float[cell_num]); // scalers grid: 0 for wall, 1 for fluid
	unique_ptr<float[]> m_(new float[cell_num]); // TODO what is it for?
};

void Fluid::addSmoke(int x, int y, float amount) {
	d_[I(x,y)] += amount;
}

void Fluid::addWind(int x, int y, float x_amount, float y_amount) {
	v_[I(x,y)] += x_amount;
	u_[I(x,y)] += y_amount;
}

void Fluid::evaluate() {
}

Fluid::~Fluid() {
	delete[] d_;
}
