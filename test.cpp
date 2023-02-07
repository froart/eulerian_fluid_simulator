#include "opengl_setup.hpp"
#include "fluid_dynamics.hpp"
#include <iostream>
#include <omp.h>
#include <memory>
#include <unistd.h>

using namespace std;

int width = 200+2;
int height = 200+2;
int cell_num = width * height;
Fluid f(1.0, width, height, 1/120.0, 100, 5);			
float* image = f.setImage();
int brush_size = 10;

void loop_code() {
	// post fluid code here
	float speed = 10;
 	for(int k = -5; k < 5; ++k) {
  	f.addWind(width/3+k, height/3, 0, speed); 
  	f.addWind(2*width/3, height/3+k, -speed, 0); 
  	f.addWind(width/3, 2*height/3+k, speed, 0); 
  	f.addWind(2*width/3+k, 2*height/3, 0, -speed); 
  }
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
