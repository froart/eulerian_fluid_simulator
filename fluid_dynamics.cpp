#include "fluid_dynamics.hpp"
#include "opengl_setup.hpp"
#include <math.h>
#include <iostream>
#include <omp.h>
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
	vector<float> v(cell_num_, 0.0);
	vector<float> u(cell_num_, 0.0);
	vector<float> v1(cell_num_, 0.0);
	vector<float> u1(cell_num_, 0.0);
	vector<float> p(cell_num_, 0.0);
	vector<float> s(cell_num_, 0.0);
	vector<float> m(cell_num_, 0.0);
	v_ = v;
	u_ = u;
	v1_ = v1;
	v1_ = v1;
	p_ = p;
	s_ = s;
	m_ = m;
};

void Fluid::addSmoke(int x, int y, float amount) {
	d_[I(x,y)] += amount;
}

void Fluid::addWind(int x, int y, float x_amount, float y_amount) {
	v_[I(x,y)] += x_amount;
	u_[I(x,y)] += y_amount;
}

void Fluid::evaluate() {
	fill_n(p_.begin(), cell_num_, 0.0);
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
