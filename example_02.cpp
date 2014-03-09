#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#pragma warning(disable: 4996)

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif

#include <time.h>
#include <math.h>

#undef Success
#include <Eigen/Core>
#include <Eigen/Dense>

#include "Cimg.h"

//****************************************************
// macros and inline functions
//****************************************************

#define PI 3.14159265  // Should be used from mathlib
#define EPS 1e-10       
#define DELETE_OBJECT(x)	if ((x)) { delete (x); (x) = NULL; }   // delete a class object
#define DELETE_ARRAY(x)		if ((x)) { delete [] (x); (x) = NULL; } // delete an array
#define FOR(i,n) for( int i=0; i<n; i++ )                           // for loop 
#define FOR_u(i, n) for (size_t i = 0; i < n; i++)                  // for loop 
#define SQUARE(x) ((x)*(x))
inline float sqr(float x) { return x*x; }
typedef unsigned char uchar; 
using namespace std;
using namespace Eigen;
using namespace cimg_library;


//****************************************************
// Basic Classes
//****************************************************

class CViewport {
  public:
    int m_w, m_h; // width and height
};

//****************************************************
// color
//****************************************************

class CColor {
public: 
	CColor() { 
		m_rgb[0] = 0.0f; m_rgb[1] = 0.0f; m_rgb[2] = 0.0f; 
	}

	CColor(float _r, float _g, float _b) { 
		m_rgb[0] = _r;  m_rgb[1] = _g; m_rgb[2] = _b; 
	}
	
	CColor Add(CColor _c) {
		return CColor(m_rgb[0]+_c.m_rgb[0], m_rgb[1]+_c.m_rgb[1], m_rgb[2]+_c.m_rgb[2]);
	}

	void Print() {
		printf("(R, G, B) = (%3.3f, %3.3f, %3.3f)\n", m_rgb[0], m_rgb[1], m_rgb[2]);
	}

public: 
	float m_rgb[3];
};


//****************************************************
// 3d vector
//****************************************************

class CVec3 {
public: 
	CVec3() {
		m_x = 0.0f; m_y = 0.0f; m_z = 0.0f; 
	}

	CVec3(float _x, float _y, float _z) {
		m_x = _x; m_y = _y; m_z = _z; 
	}

	CVec3 Add(CVec3 _v) {
		return CVec3(m_x + _v.m_x, m_y + _v.m_y, m_z + _v.m_z);
	}

	CVec3 Subtract(CVec3 _v) {
		return CVec3(m_x - _v.m_x, m_y - _v.m_y, m_z - _v.m_z);
	}

	CVec3 Scale(float _s) { 
		return CVec3(m_x * _s, m_y * _s, m_z * _s);
	}

	float Dot(CVec3 _v) {
		return m_x * _v.m_x + m_y * _v.m_y + m_z * _v.m_z; 
	}

	CVec3 Cross(CVec3 _v) {
		float x = m_y * _v.m_z - m_z * _v.m_y; // u2v3-u3v2 
		float y = m_z * _v.m_x - m_x * _v.m_z; // u3v1-u1v3 
		float z = m_x * _v.m_y - m_y * _v.m_x; // v1v2-u2v1
		return CVec3(x, y, z);
	}

	float Norm() {
		return sqrt(SQUARE(m_x) + SQUARE(m_y) + SQUARE(m_z)); 
	}

	CVec3 UnitVec() {
		float n = Norm();
		return Scale(1/(n+EPS));
	}

	void Print() {
		printf("(x, y, z) = (%3.3f, %3.3f, %3.3f)\n", m_x, m_y, m_z);
	}
public: 
	float m_x, m_y, m_z; 
};

//****************************************************
// Light
//****************************************************
class CLight {
public: 
	enum ELightType {
		Point,
		Directional,
	};

	CLight() {}

	CLight(CVec3 _dir, CColor _color, ELightType _type) {
		m_dir = _dir; m_color = _color; m_type =  _type; 
	}

	void Print() {
		if (m_type == Point) 
			printf("Point Light ");
		else 
			printf("Directional Light ");
		printf("(x, y, z, R, G, B) = (%2.2f, %2.2f, %2.2f, %2.2f, %2.2f, %2.2f)\n", 
			m_dir.m_x, m_dir.m_y, m_dir.m_z, m_color.m_rgb[0], m_color.m_rgb[1], m_color.m_rgb[2]);
	}

public: 
	ELightType m_type;  // light type
	CVec3 m_dir;        //  location/position/direction
	CColor m_color;    // light color
};

//****************************************************
// shape 
//****************************************************
class C3DModel {
public: 
	//virtual CShape() = 0;
	virtual CVec3 Norm(CVec3 _pos) = 0;
	virtual bool IsVisible(float _x, float _y) = 0; 
	virtual float Z(float _x, float _y) = 0;
	virtual CVec3 Center() = 0;
	virtual float Radius() = 0;
public: 

};

class CSphere : public C3DModel {
public: 
	CSphere() {}
	CSphere(CVec3 _center, float _radius) {
		m_center = _center; 
		m_radius = _radius;
	}

	CVec3 Norm(CVec3 _pos) {
		return _pos.Subtract(m_center).UnitVec();
	}

	bool IsVisible(float _x, float _y) {
		return SQUARE(_x-m_center.m_x) + SQUARE(_y-m_center.m_y) < SQUARE(m_radius);
	}

	float Z(float _x, float _y) {
		return sqrt(SQUARE(m_radius) - SQUARE(_x-m_center.m_x) - SQUARE(_y-m_center.m_y))+m_center.m_z; 
	}

	float Radius() { return m_radius;}
	CVec3 Center() { return m_center; }
public: 
	CVec3 m_center; 
	float m_radius; 
};

// class Screen {
// public:
// 	Screen(){ screen.resize(256, 256); m_center[0]=128; m_center[1]=128;}
// 	Screen(int _w, int _h){
// 		m_w
// 	}
	
// public:
// 	int m_center[2];
// 	int m_w, m_h;

// 	pngwriter screen(1,1,0,);

// }



//****************************************************
// Global Variables
//****************************************************

CViewport g_viewport;
CColor g_Ka, g_Kd, g_Ks;  // parameters 
vector<CLight> g_lights;
vector<C3DModel*> g_models; 
float g_v; 
bool g_isMultiple = false; 
bool g_isToon = false;  
bool g_isSave = true; 

CVec3 d; //anisotropic fibrics direction vector
bool g_isAniso = false;


//USING_PART_OF_NAMESPACE_EIGEN

int main(int argc, char *argv[]){
	Matrix3f m3;
	m3 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	Matrix4f m4=Matrix4f::Identity();
	Vector4i v4(1,2,3,4);
	cout<<"m3\n" << m3 << "\nm4:\n"
		<<m4 << "\nv4:\n" <<v4<<endl;

	CImg<unsigned char> image("lena.pgm"), visu(500,400,1,3,0);
  const unsigned char red[] = { 255,0,0 }, green[] = { 0,255,0 }, blue[] = { 0,0,255 };
  image.blur(2.5);
  CImgDisplay main_disp(image,"Click a point"), draw_disp(visu,"Intensity profile");
  while (!main_disp.is_closed() && !draw_disp.is_closed()) {
    main_disp.wait();
    if (main_disp.button() && main_disp.mouse_y()>=0) {
      const int y = main_disp.mouse_y();
      visu.fill(0).draw_graph(image.get_crop(0,y,0,0,image.width()-1,y,0,0),red,1,1,0,255,0);
      visu.draw_graph(image.get_crop(0,y,0,1,image.width()-1,y,0,1),green,1,1,0,255,0);
      visu.draw_graph(image.get_crop(0,y,0,2,image.width()-1,y,0,2),blue,1,1,0,255,0).display(draw_disp);
      }
    }
  return 0;

}
