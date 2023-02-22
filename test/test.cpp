#include "../opengl_setup.hpp"
#include "../fluid_dynamics.hpp"
#include <iostream>
#include <omp.h>
#include <memory>
#include <unistd.h>
#include <cmath>
#include <tuple>

using namespace std;

int width = 300+2;
int height = 300+2;
int cell_num = width * height;
Fluid f(1.0, width, height, 1/60.0, 100, 5);			
float* image = f.setImage();
int brush_size = 30;
float speed = 100;

void circular_vortex(float R, float x0, float y0, float speed) {
	// initial
	float u0 = 0;
	float v0 = speed;
	float xp = R*cos(0)+x0;
	float yp = R*sin(0)+y0;
	float s = 2*numbers::pi/360;
	f.addWind(xp, yp, u0, v0);
	// circular vortex
	for(int i = 1; i < 360; ++i) {		
		float x = R*cos(i*s)+x0;
		float y = R*sin(i*s)+y0;
		float dx = x - xp;
		float dy = y - yp;
		float a = atan(dy/dx); // angle of a tangent
		a = i < 180 ? -a : a;
		float u = speed*cos(a);
		float v = speed*sin(a);
		f.addWind(x, y, u, v);
		xp = x;
		yp = y;
  }
}

void loop_code() {
	// post fluid code here
//for(int k = -10; k < 10; ++k) {
//	f.addWind(width/3+k, height/3, 0, speed); 
//	f.addWind(2*width/3, height/3+k, -speed, 0); 
//	f.addWind(width/3, 2*height/3+k, speed, 0); 
//	f.addWind(2*width/3+k, 2*height/3, 0, -speed); 
//	f.addWind(width/2+15, height/2+k, speed, 0); 
//	f.addWind(width/2-15, height/2+k, -speed, 0); 
//	f.addWind(width/2+k, height/2+15, 0, speed); 
//	f.addWind(width/2+k, height/2-15, 0, -speed); 
	circular_vortex(20, width/2, height/2, speed);
	f.evaluate();
	return;
}



int main() {
	runSimulation("Eulerian Fluid Simulator");
	return 0;
}

void mouse(int button, int state, int x, int y){
	if(button == GLUT_RIGHT_BUTTON)
		f.clearImage();
	if(button == GLUT_LEFT_BUTTON) {
		for(int i = -brush_size/2; i <= brush_size/2; ++i)
			for(int j = -brush_size/2; j <= brush_size/2; ++j)
				f.addSmoke(x+i, height-y+j, 1.0);
	}
}

void keyboard(unsigned char c, int x, int y) {
  if(c == 27) {
		exit(0);
	}
}
