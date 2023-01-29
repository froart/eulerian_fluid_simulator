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
						 int w_, 
						 int h_, 
						 float dt,
						 float dens,
						 int iterations) 
						 : m_(image),
							 dens_(dens),
							 dt_(dt),
							 it_(iterations) {
	int num_x = w_ + 2;
	int num_y = h_ + 2;
	int cell_num = num_x * num_y;
	vector<float> u(cell_num_, 0.0);
	vector<float> v(cell_num_, 0.0);
	vector<float> u1(cell_num_, 0.0);
	vector<float> v1(cell_num_, 0.0);
	vector<float> p(cell_num_, 0.0);
	vector<float> s(cell_num_, 0.0);
	u_ = u;
	v_ = v;
	u1_ = u1;
	v1_ = v1;
	p_ = p;
	s_ = s;
	over_relaxation_ = 1.9;
};

void Fluid::addSmoke(int x, int y, float amount) {
}

void Fluid::addWind(int x, int y, float x_amount, float y_amount) {
}

void Fluid::evaluate() {
	fill_n(p_.begin(), cell_num_, 0.0);
	project();
	extrapolate();
	advect_velocity();
	advect_smoke();
}

void Fluid::project() { // force imcompressibility
	for(int k; k < it_; k++)
		for(int j = 1; j < h_-1; j++)
			for(int i = 1; i < w_-1; i++) {
				if(s_[I(i,j)] == 0.0) continue; // if evaluating at wall

				float sx0 = s_[I(i-1,j)];
				float sx1 = s_[I(i+1,j)];
				float sy0 = s_[I(i,j-1)];
				float sy1 = s_[I(i,j+1)];
				float s = sx0 + sx1 + sy0 + sy1;
				if(s == 0.0) continue; // if nothing to compute here...

				float div = u_[I(i+1,j)] - u_[I(i,j)] + v_[I(i,j+1)] - v_[I(i,j)];
				// computing the pressure (not necessary for the simulation
				// FIXME  should be minus???
				float p = -(div / s) * over_relaxation_;
				// TODO where does this equation come from???
				p_[I(i,j)] += p * (dens_ * h_) / dt_; 
				// fix the fluid to be imcompressible
				u_[I(i,j)] -= sx0 * p;
				u_[I(i+1,j)] += sx1 * p;
				v_[I(i,j)] -= sy0 * p;
				v_[I(i,j+1)] += sy1 * p;
			}
}

void Fluid::extrapolate() { // border conditions
	for(int i = 0; i < w_; ++i) {
		u_[I(i,0)] = u_[I(i,1)]; // left wall
		u_[I(i,h_-1)] = u_[I(i,h_-2)]; // right wall
	}
	for(int j = 0; j < h_; ++j) {
		v_[I(0,j)] = v_[I(1,j)]; // bottom wall
		v_[I(w_-1,j)] = v_[I(w_-2,j)]; // top wall
	}
}

void Fluid::advect_velocity() {
}

void Fluid::advect_smoke() {
}
