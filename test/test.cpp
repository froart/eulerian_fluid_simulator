#include <vt_framebuffer.hpp>
#include "../include/eulerian_fluid_dynamics.hpp"
#include <string>

int main() 
{
  const unsigned int width  = 300;
  const unsigned int height = 300;

  vt::FrameBuffer framebuffer( "Gas simulation", width, height );
  Gas gas( 0.1, width, height, 0.1, 10.0, 10 ); 

  framebuffer.bind( gas.getSmokeField() );

  for( int j = height / 2; j < height / 2 + 10; ++j )
  {
    for( int i = width / 2; i < width / 2 + 10; ++i )
    {
      gas.addSmoke( j, i, 1.0 );
    }
  }

  while( !framebuffer.requestedToExit() )
  {
    gas.evaluate();
    framebuffer.update();
    for( int i = width / 2; i < width / 2 + 10; ++i )
      gas.addWind( height / 2, i, -1000.0, 1100.0 );
  }

  return 0;
}
