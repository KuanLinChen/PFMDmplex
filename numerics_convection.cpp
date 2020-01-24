#include "numerics_convection.hpp"
CCentJST::CCentJST() 
{
  
}
void CCentJST::Init( CGeometry *mm, CVariable *vv )
{
  s = new CSysSolve() ;

  m = mm ;
  var = vv ;
  s->Init( &dmCell, &m->gVec, m->cEndInterior) ;

}

