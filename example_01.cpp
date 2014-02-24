#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#pragma warning(disable: 4996)
// #include <conio.h>
#include "stb_image_write.h"
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
//****************************************************
// Simple init function
//****************************************************

void initScene(){  // set models
	C3DModel* sphere1 = new CSphere(CVec3(0.5f, 0.5f, 0.5f), 1.0/ 3.0);
	g_models.push_back(sphere1);
	if (g_isMultiple) {

	}
  // Nothing to do here for this simple example.

}


//****************************************************
// reshape viewport if the window is resized
//****************************************************
void myReshape(int w, int h) {
  g_viewport.m_w = w;
  g_viewport.m_h = h;

  glViewport (0,0,g_viewport.m_w,g_viewport.m_h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0, g_viewport.m_w, 0, g_viewport.m_h);

}


//****************************************************
// A routine to set a pixel by drawing a GL point.  This is not a
// general purpose routine as it assumes a lot of stuff specific to
// this example.
//****************************************************

void setPixel(int x, int y, GLfloat r, GLfloat g, GLfloat b) {
  glColor3f(r, g, b);
  glVertex2f(x, y);   // The 0.5 is to target pixel
  // centers 
  // Note: Need to check for gap
  // bug on inst machines.
}


// TO-DO: RenderSpecularTerm
// TO-DO: point light 
CColor RenderAmbientTerm(CLight _light, CColor _Ka) {
	CColor c;
	FOR (k, 3) {
		c.m_rgb[k] = _Ka.m_rgb[k] * _light.m_color.m_rgb[k];
	}
	return c; 
}

CColor RenderDiffuseTerm(CLight _light, CColor _Kd, CVec3 _n) {
	CColor c;
	float dot = _n.Dot(_light.m_dir.UnitVec());
	/*if (dot > 0.0) {
		_light.m_dir.UnitVec().Print();
		_n.Print();
	}*/
	dot = max(dot, 0.0f);
	//printf("dot = %2.2f\n", dot);
	FOR (k, 3)
		c.m_rgb[k] = _Kd.m_rgb[k] * dot * _light.m_color.m_rgb[k];

	return c; 
}

CColor RenderSpecularTerm(CLight _light, CColor _Ks, CVec3 _n, CVec3 _v, float _p) {
	CColor c;
	_light.m_dir = _light.m_dir.UnitVec();
    CVec3 r = _n.Scale(2*_light.m_dir.Dot(_n)).Subtract(_light.m_dir);
	//r = r.UnitVec();
    //r.Print();
    //printf("\t%f\n",r.Norm());
    float dot = pow(max(r.Dot(_v.UnitVec()),0.0f), _p);

	FOR (k, 3) 
		c.m_rgb[k] = _Ks.m_rgb[k] * dot * _light.m_color.m_rgb[k];

	return c; 
}


CColor aniso_specular(CLight _light, CVec3 _n, CVec3 _v, CVec3 _d){
	CVec3 T = _n.UnitVec().Cross(_d.UnitVec().Cross(_n.UnitVec()).UnitVec());

	CColor K;
	return K;

}

CColor Render(CLight _light, CVec3 _n, CVec3 _v, CVec3 _d, float _p, CColor _Ka, CColor _Kd, CColor _Ks, bool Aniso) {
	CColor c_all(0.0f, 0.0f, 0.0f);

	CColor c_ambient = RenderAmbientTerm(_light, _Ka); 
	c_all = c_all.Add(c_ambient);
	
	CColor c_diffuse = RenderDiffuseTerm(_light, _Kd, _n);
	c_all = c_all.Add(c_diffuse);

	if (Aniso){
		_Ks=aniso_specular(_light, _n, _v, _d);
	}

	CColor c_specular = RenderSpecularTerm(_light, _Ks, _n, _v, _p);
	c_all = c_all.Add(c_specular);
	
	return c_all; 
}



//****************************************************
// function that saves image to disk
//***************************************************
void Save() {
	printf("save image\n");
	int w = g_viewport.m_w;
	int h = g_viewport.m_h;
	uchar* buffer = new uchar[w * h * 3];
	glReadPixels (0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, buffer);
	stbi_write_png ("shading.png", w, h, 3, buffer + (w * 3 * (h - 1)), -3 * w);
	DELETE_ARRAY(buffer);
}

//****************************************************
// function that does the actual drawing of stuff
//***************************************************
void myDisplay() {
  glClear(GL_COLOR_BUFFER_BIT);				// clear the color buffer
  glMatrixMode(GL_MODELVIEW);			        // indicate we are specifying camera transformations
  glLoadIdentity();				        // make sure transformation is "zero'd"
  glBegin(GL_POINTS);
    
    
  CVec3 Origin = g_models[0]->Center();
  Origin.Print();
  CVec3 Eye(0.5f, 0.5f, 1.5f); 

  
  // Start drawing
  FOR (w, g_viewport.m_w) {
	  FOR (h, g_viewport.m_h) {
		  //w = 67; h = 191;
		  float x = w/(float)g_viewport.m_w; 
		  float y = h/(float)g_viewport.m_h; 
		
		  CColor color(0.0, 0.0, 0.0);

		  

		  FOR_u (k, g_models.size()) { // draw k-th model
			  C3DModel* model = g_models[k];

			  if (model->IsVisible(x, y)) { // visibility test
				 
				  float z = model->Z(x, y);   // compute z: bug could be 2 z 
				  CVec3 n = model->Norm(CVec3(x, y, z));  // compute norm
				  CVec3 v = Eye.Subtract(CVec3(x, y, z)); // compute eye
				  FOR_u (l, g_lights.size()) {
					  CLight light = g_lights[l];
					  light.m_dir = light.m_dir.Scale(1/model->Radius());
					  light.m_dir = light.m_dir.Add(Origin);
					  if (light.m_type == CLight::Point){
						
						  light.m_dir = light.m_dir.Subtract(CVec3(x, y, z));
                      } else 
						  light.m_dir = light.m_dir.Scale(-1);



					  color = color.Add(Render(light, n, v, d, g_v, g_Ka, g_Kd, g_Ks, g_isAniso));
				  }
			  }
		  }
		  setPixel(w, h, color.m_rgb[0], color.m_rgb[1], color.m_rgb[2]);
	  }
  }

  glEnd();
  if (g_isSave)
	Save();                   // save image to disk  
  glFlush();
  glutSwapBuffers();		// swap buffers (we earlier set double buffer)

  // if (getch() == ' ')     // press key to quit               
	 //  exit(0);
}

void TestParams() {
	//CLight light(CVec3(-2.0f, 2.0f, 2.0f), CColor(1.0f, 1.0f, 1.0f), CLight::Point);
	//g_lights.push_back(light);
	///*g_Ka.m_rgb[0] = 0.5f; 
	//g_Ka.m_rgb[1] = 0.5f; 
	//g_Ka.m_rgb[2] = 0.5f; */

	//g_Kd.m_rgb[0] = 1.0f; 
	//g_Kd.m_rgb[1] = 1.0f; 
	//g_Kd.m_rgb[2] = 1.0f; 
}


//****************************************************
// parse command lines
//****************************************************
void ParseParams(int argc, char* argv[]) {
	int count = 0; 
	/*for (int i = 1; i < argc; i++)
		printf("%s ", argv[i]);*/
	//printf("\n");
	while (count < argc) {
		//printf("%s ", argv[count]);
		if (strcmp(argv[count], "-ka") == 0) {
			printf("%s\n", argv[count+1]);
			FOR (k, 3)
				sscanf(argv[count+k+1], "%f", &g_Ka.m_rgb[k]);
			count+=3;
		}

		if (strcmp(argv[count], "-kd") == 0) {
			FOR (k, 3)
				sscanf(argv[count+k+1], "%f", &g_Kd.m_rgb[k]);
			count += 3;
		}

		if (strcmp(argv[count], "-ks") == 0) {
			FOR (k, 3)
				sscanf(argv[count+k+1], "%f", &g_Ks.m_rgb[k]);
			count += 3;
		}

		if (strcmp(argv[count], "-sp") == 0) {
			sscanf(argv[count+1], "%f", &g_v);
			count += 1;
		}

		if (strcmp(argv[count], "-pl") == 0) {
			float x, y, z; 
			sscanf(argv[count+1], "%f", &x);
			sscanf(argv[count+2], "%f", &y);
			sscanf(argv[count+3], "%f", &z);
			CColor c; 
			sscanf(argv[count+4], "%f", &c.m_rgb[0]);
			sscanf(argv[count+5], "%f", &c.m_rgb[1]);
			sscanf(argv[count+6], "%f", &c.m_rgb[2]);
			CLight light(CVec3(x, y, z), c, CLight::Point);
			g_lights.push_back(light);
			count += 6;
		}

		if (strcmp(argv[count], "-dl") == 0) {
			float x, y, z; 
			sscanf(argv[count+1], "%f", &x);
			sscanf(argv[count+2], "%f", &y);
			sscanf(argv[count+3], "%f", &z);
			CColor c; 
			sscanf(argv[count+4], "%f", &c.m_rgb[0]);
			sscanf(argv[count+5], "%f", &c.m_rgb[1]);
			sscanf(argv[count+6], "%f", &c.m_rgb[2]);
			CLight light(CVec3(x, y, z), c, CLight::Directional);
			g_lights.push_back(light);
			count += 6;
		}

		if (strcmp(argv[count], "-m") == 0) 
			g_isMultiple = true; 

		if (strcmp(argv[count], "-t") == 0)
			g_isToon = true; 


		if (strcmp(argv[count], "-s") == 0)
			g_isSave = true; 

		count++;

		if (strcmp(argv[count], "-a") == 1){
			g_isAniso = true;
			
				sscanf(argv[count+1], "%f", &d.m_x);
				sscanf(argv[count+2], "%f", &d.m_y);
				sscanf(argv[count+3], "%f", &d.m_z);
				count += 3;
		}
	}

	printf("--------------------------------------------\n");
	printf("parameter list\n");
	printf("Ambient Term Ka:\n"); // (%2.2f %2.2f %2.2f)", g_Ka.m_rgb[0], g_Ka.m_rgb[1], g_Ka.m_rgb[2]);
	g_Ka.Print();
	printf("Diffuse Term Kd:\n");// (%2.2f %2.2f %2.2f)", g_Kd.m_rgb[0], g_Kd.m_rgb[1], g_Kd.m_rgb[2]);
	g_Kd.Print();
	printf("Specular Term Ks:\n"); //, g_Ks.m_rgb[0], g_Ks.m_rgb[1], g_Ks.m_rgb[2]);
	g_Ks.Print();
	printf("Power v: %2.2f\n", g_v); //, g_Ks.m_rgb[0], g_Ks.m_rgb[1], g_Ks.m_rgb[2]);
	int nLights = (int)g_lights.size();
	printf("(%d) lights:\n", nLights);
	FOR (n, nLights)
		g_lights[n].Print();

	if (g_isToon)
		printf("enable toon shading\n");

	printf("--------------------------------------------\n");
}

//****************************************************
// the usual stuff, nothing exciting here
//****************************************************
int main(int argc, char *argv[]) {
    //This initializes glut
  //TestParams(); // for debug 
  glutInit(&argc, argv);

  //This tells glut to use a double-buffered window with red, green, and blue channels 
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

  // Initialize the viewport size
  g_viewport.m_w = 400;
  g_viewport.m_h = 400;

  //The size and position of the window
  glutInitWindowSize(g_viewport.m_w, g_viewport.m_h);
  glutInitWindowPosition(0,0);
  glutCreateWindow(argv[0]);

  initScene();							// quick function to set up scene
  
  // parse parameters 
  ParseParams(argc, argv);

  glutDisplayFunc(myDisplay);				// function to run when its time to draw something
  glutReshapeFunc(myReshape);				// function to run when the window gets resized
  glutMainLoop();							// infinite loop that will keep drawing and resizing
  // and whatever else

  // clear global variables
  g_lights.clear(); 
  FOR_u (i, g_models.size())
	  DELETE_OBJECT(g_models[i]); 
  g_models.clear(); 

  return 0;
}








