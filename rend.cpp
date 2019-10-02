/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include <algorithm>

#define PI (float) 3.14159265358979323846


void copyMatrix(GzMatrix a, GzMatrix b) {
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			a[i][j] = b[i][j];
}

void multiMatrix(GzMatrix a, GzMatrix c, GzMatrix b) {
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			a[i][j] = c[i][0] * b[0][j] + c[i][1] * b[1][j] + c[i][2] * b[2][j] + c[i][3] * b[3][j];
}

void norm(GzCoord v) {

	float n = sqrt((v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]));

	v[0] = v[0] / n;
	v[1] = v[1] / n;
	v[2] = v[2] / n;
}

int GzRender::GzRotXMat(float degree, GzMatrix mat)
{
/* HW 3.1
// Create rotate matrix : rotate along x axis
// Pass back the matrix using mat value
*/
	degree *= (PI / 180);

	mat[0][0] = 1;
	mat[1][1] = cos(degree);
	mat[1][2] = -sin(degree);
	mat[2][1] = sin(degree);
	mat[2][2] = cos(degree);
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRender::GzRotYMat(float degree, GzMatrix mat)
{
/* HW 3.2
// Create rotate matrix : rotate along y axis
// Pass back the matrix using mat value
*/
	degree *= (PI / 180);

	mat[1][1] = 1;
	mat[0][0] = cos(degree);
	mat[0][2] = sin(degree);
	mat[2][0] = -sin(degree);
	mat[2][2] = cos(degree);
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRender::GzRotZMat(float degree, GzMatrix mat)
{
/* HW 3.3
// Create rotate matrix : rotate along z axis
// Pass back the matrix using mat value
*/

	degree *= (PI / 180);

	mat[2][2] = 1;
	mat[0][0] = cos(degree);
	mat[0][1] = -sin(degree);
	mat[1][0] = sin(degree);
	mat[1][1] = cos(degree);
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRender::GzTrxMat(GzCoord translate, GzMatrix mat)
{
/* HW 3.4
// Create translation matrix
// Pass back the matrix using mat value
*/
	mat[0][0] = 1;
	mat[1][1] = 1;
	mat[2][2] = 1;
	mat[3][3] = 1;
	mat[0][3] = translate[0];
	mat[1][3] = translate[1];
	mat[2][3] = translate[2];
	return GZ_SUCCESS;
}


int GzRender::GzScaleMat(GzCoord scale, GzMatrix mat)
{
/* HW 3.5
// Create scaling matrix
// Pass back the matrix using mat value
*/
	mat[0][0] = scale[0];
	mat[1][1] = scale[1];
	mat[2][2] = scale[2];
	mat[3][3] = 1;

	return GZ_SUCCESS;
}


GzRender::GzRender(int xRes, int yRes)
{
/* HW1.1 create a framebuffer for MS Windows display:
 -- set display resolution
 -- allocate memory for framebuffer : 3 bytes(b, g, r) x width x height
 -- allocate memory for pixel buffer
 */
	xres = xRes;
	yres = yRes;
	pixelbuffer = new GzPixel[xres * yres];
	framebuffer = new char[xRes * yRes * 3];
	

/* HW 3.6
- setup Xsp and anything only done once 
- init default camera 
*/ 

	matlevel = 0;

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++)
			Xsp[i][j] = 0;
	}

	Xsp[0][0] = Xsp[0][3] = xres / 2;
	Xsp[1][1] = -yres / 2;
	Xsp[1][3] = yres / 2;
	Xsp[2][2] = INT_MAX;
	Xsp[3][3] = 1;

	m_camera.position[0] = DEFAULT_IM_X;
	m_camera.position[1] = DEFAULT_IM_Y;
	m_camera.position[2] = DEFAULT_IM_Z;

	m_camera.lookat[0] = 0.0;
	m_camera.lookat[1] = 0.0;
	m_camera.lookat[2] = 0.0;

	m_camera.worldup[0] = 0.0;
	m_camera.worldup[1] = 1.0;
	m_camera.worldup[2] = 0.0;

	m_camera.FOV = DEFAULT_FOV;

}

GzRender::~GzRender()
{
/* HW1.2 clean up, free buffer memory */
	delete[] pixelbuffer;
	delete[] framebuffer;
}

int GzRender::GzDefault()
{
/* HW1.3 set pixel buffer to some default values - start a new frame */

	for (int i = 0; i < xres * yres; i++) {
		pixelbuffer[i] = { 128 * 16, 112 * 16, 96 * 16, 0, INT_MAX };
	}
	return GZ_SUCCESS;

}

int GzRender::GzBeginRender()
{
/* HW 3.7 
- setup for start of each frame - init frame buffer color,alpha,z
- compute Xiw and projection xform Xpi from camera definition 
- init Ximage - put Xsp at base of stack, push on Xpi and Xiw 
- now stack contains Xsw and app can push model Xforms when needed 
*/ 

	matlevel = 0;

	GzCoord z = { m_camera.lookat[0] - m_camera.position[0],  m_camera.lookat[1] - m_camera.position[1] , m_camera.lookat[2] - m_camera.position[2] };
	norm(z);
	float productZU = z[0] * m_camera.worldup[0] + z[1] * m_camera.worldup[1] + z[2] * m_camera.worldup[2];
	GzCoord c = { m_camera.worldup[0] - z[0] * productZU, m_camera.worldup[1] - z[1] * productZU, m_camera.worldup[2] - z[2] * productZU };
	norm(c);
	GzCoord x = { c[1] * z[2] - c[2] * z[1], c[2] * z[0] - c[0] * z[2], c[0] * z[1] - c[1] * z[0] };
	norm(x);

	GzMatrix Xiw = {
		{ x[0], x[1], x[2], -(x[0] * m_camera.position[0] + x[1] * m_camera.position[1] + x[2] * m_camera.position[2])},
		{ c[0], c[1], c[2], -(c[0] * m_camera.position[0] + c[1] * m_camera.position[1] + c[2] * m_camera.position[2])},
		{ z[0], z[1], z[2], -(z[0] * m_camera.position[0] + z[1] * m_camera.position[1] + z[2] * m_camera.position[2])},
		{ 0, 0, 0, 1}
	};

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++)
			m_camera.Xiw[i][j] = Xiw[i][j];
	}

	GzMatrix Xpi = {
		{ 1,0,0,0 },
		{ 0,1,0,0 },
		{ 0,0,tan(m_camera.FOV * PI / 180 / 2),0 },
		{ 0,0,tan(m_camera.FOV * PI / 180 / 2),1 }
	};

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++)
			m_camera.Xpi[i][j] = Xpi[i][j];
	}


	Ximage[matlevel - 1];

	GzPushMatrix(Xsp);
	GzPushMatrix(Xpi);
	GzPushMatrix(Xiw);

	return GZ_SUCCESS;
}

int GzRender::GzPutCamera(GzCamera camera)
{
/* HW 3.8 
/*- overwrite renderer camera structure with new camera definition
*/

	m_camera.position[0] = camera.position[0];
	m_camera.position[1] = camera.position[1];
	m_camera.position[2] = camera.position[2];

	m_camera.lookat[0] = camera.lookat[0];
	m_camera.lookat[1] = camera.lookat[1];
	m_camera.lookat[2] = camera.lookat[2];

	m_camera.worldup[0] = camera.worldup[0];
	m_camera.worldup[1] = camera.worldup[1];
	m_camera.worldup[2] = camera.worldup[2];

	m_camera.FOV = camera.FOV;

	return GZ_SUCCESS;	
}

int GzRender::GzPushMatrix(GzMatrix	matrix)
{
/* HW 3.9 
- push a matrix onto the Ximage stack
- check for stack overflow
*/
	if (matlevel == 0)
		 copyMatrix(Ximage[0], matrix);
	
	else
		multiMatrix(Ximage[matlevel], Ximage[matlevel - 1], matrix);

	matlevel++;
	return GZ_SUCCESS;
}

int GzRender::GzPopMatrix()
{
/* HW 3.10
- pop a matrix off the Ximage stack
- check for stack underflow
*/
	if (matlevel > 0)
		matlevel--;
	return GZ_SUCCESS;
}

int GzRender::GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
/* HW1.4 write pixel values into the buffer */

	if (i >= 0 && i < xres && j >= 0 && j < yres) {
		if (z < pixelbuffer[j * xres + i].z)
			pixelbuffer[j * xres + i] = { min(4095, max(r, 0)), min(4095, max(g, 0)), min(4095, max(b, 0)), min(4095, max(a, 0)), z };
	}
	return GZ_SUCCESS;
}


int GzRender::GzGet(int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
/* HW1.5 retrieve a pixel information from the pixel buffer */

	if (i >= 0 && i < xres && j >= 0 && j < yres) {
		*r = pixelbuffer[i * yres + j].red;
		*g = pixelbuffer[i * yres + j].green;
		*b = pixelbuffer[i * yres + j].blue;
		*a = pixelbuffer[i * yres + j].alpha;
		*z = pixelbuffer[i * yres + j].z;
	}

	return GZ_SUCCESS;
}


int GzRender::GzFlushDisplay2File(FILE* outfile)
{
/* HW1.6 write image to ppm file -- "P6 %d %d 255\r" */

	fprintf(outfile, "P6 %d %d 255\r", xres, yres);
	for (int i = 0; i < xres * yres; i++) {
		char r = pixelbuffer[i].red >> 4;
		char g = pixelbuffer[i].green >> 4;
		char b = pixelbuffer[i].blue >> 4;
		fwrite(&r, 1, 1, outfile);
		fwrite(&g, 1, 1, outfile);
		fwrite(&b, 1, 1, outfile);
	}
	
	return GZ_SUCCESS;
}

int GzRender::GzFlushDisplay2FrameBuffer()
{
/* HW1.7 write pixels to framebuffer: 
	- put the pixels into the frame buffer
	- CAUTION: when storing the pixels into the frame buffer, the order is blue, green, and red 
	- NOT red, green, and blue !!!
*/
	for (int i = 0; i < xres * yres; i++) {
		framebuffer[i * 3 + 2] = (char)(pixelbuffer[i].red >> 4);
		framebuffer[i * 3 + 1] = (char)(pixelbuffer[i].green >> 4);
		framebuffer[i * 3] = (char)(pixelbuffer[i].blue >> 4);
	}

	return GZ_SUCCESS;
}


/***********************************************/
/* HW2 methods: implement from here */

int GzRender::GzPutAttribute(int numAttributes, GzToken	*nameList, GzPointer *valueList) 
{
/* HW 2.1
-- Set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
-- In later homeworks set shaders, interpolaters, texture maps, and lights
*/
	for (int i = 0; i < numAttributes; i++) {
		if (nameList[i] == GZ_RGB_COLOR) {
			float* color = (float*)(valueList[i]);
			flatcolor[0] = color[0];
			flatcolor[1] = color[1];
			flatcolor[2] = color[2];
		}

	}
	return GZ_SUCCESS;
}

int GzRender::GzPutTriangle(int numParts, GzToken *nameList, GzPointer *valueList)
/* numParts - how many names and values */
{
/* HW 2.2
-- Pass in a ((GzCoord*)valueList[i])angle description with tokens and values corresponding to
      GZ_NULL_TOKEN:		do nothing - no values
      GZ_POSITION:		3 vert positions in model space
-- Invoke the rastrizer/scanline framework
-- Return error code
*/

	for (int i = 0; i < numParts; i++) {
		if (nameList[i] == GZ_POSITION) {
			point value[3];
			float W1;
			for (int j = 0; j < 3; ++j) {
				// xform verticies
				value[j].x = (double)Ximage[matlevel - 1][0][0] * ((GzCoord*)valueList[i])[j][0] + Ximage[matlevel - 1][0][1] * ((GzCoord*)valueList[i])[j][1] + Ximage[matlevel - 1][0][2] * ((GzCoord*)valueList[i])[j][2]
					+ Ximage[matlevel - 1][0][3];
				value[j].y = (double)Ximage[matlevel - 1][1][0] * ((GzCoord*)valueList[i])[j][0] + Ximage[matlevel - 1][1][1] * ((GzCoord*)valueList[i])[j][1] + Ximage[matlevel - 1][1][2] * ((GzCoord*)valueList[i])[j][2]
					+ Ximage[matlevel - 1][1][3] * 1.0;
				value[j].z = (double)Ximage[matlevel - 1][2][0] * ((GzCoord*)valueList[i])[j][0] + Ximage[matlevel - 1][2][1] * ((GzCoord*)valueList[i])[j][1] + Ximage[matlevel - 1][2][2] * ((GzCoord*)valueList[i])[j][2]
					+ Ximage[matlevel - 1][2][3] * 1.0;
				W1 = (double)Ximage[matlevel - 1][3][0] * ((GzCoord*)valueList[i])[j][0] + Ximage[matlevel - 1][3][1] * ((GzCoord*)valueList[i])[j][1] + Ximage[matlevel - 1][3][2] * ((GzCoord*)valueList[i])[j][2]
					+ Ximage[matlevel - 1][3][3] * 1.0;
				value[j].x /= W1;
				value[j].y /= W1;
				value[j].z /= W1;
			}

			for (int j = 0; j < 3; j++) {
				if (value[j].z < 0)
					break;
			}

			sort_y(value);

			float A, B, C, D;
			float x10, x20, y10, y20, z10, z20;
			x10 = (value[1].x - value[0].x);
			y10 = (value[1].y - value[0].y);
			z10 = (value[1].z - value[0].z);
			x20 = (value[2].x - value[0].x);
			y20 = (value[2].y - value[0].y);
			z20 = (value[2].z - value[0].z);
			A = ((y10 * z20) - (y20 * z10));
			B = ((x20 * z10) - (x10 * z20));
			C = ((x10 * y20) - (x20 * y10));
			D = ((value[0].y * (B)) + (value[0].x * (A)) + (value[0].z * (C))) * (-1);

			float xt;
			if (value[0].x == value[2].x)
				xt = value[0].x;
			else
				xt = (value[0].x - value[2].x) * (value[1].y - value[2].y) / (value[0].y - value[2].y) + value[2].x;
			float x1 = min(xt, value[1].x);
			float x2 = max(xt, value[1].x);
			for (int i = ((int)(value[0].y)) + 1; i < ((int)(value[1].y)) + 1; i++) {
				int xl = (int)(((float)i - value[0].y) / (value[1].y - value[0].y) * (x1 - value[0].x) + value[0].x) + 1;
				int xr = (int)(((float)i - value[0].y) / (value[1].y - value[0].y) * (x2 - value[0].x) + value[0].x) + 1;
				for (int j = xl; j < xr; j++) {
					float z;
					z = -(A * j + B * i + D) / C;
					GzPut(j, i, ctoi(flatcolor[RED]), ctoi(flatcolor[GREEN]), ctoi(flatcolor[BLUE]), 1, (int)z);
				}
			}

			for (int i = ((int)(value[1].y)) + 1; i < ((int)(value[2].y)) + 1; i++) {
				int xl = (int)(((float)i - value[1].y) / (value[2].y - value[1].y) * (value[2].x - x1) + x1) + 1;
				int xr = (int)(((float)i - value[1].y) / (value[2].y - value[1].y) * (value[2].x - x2) + x2) + 1;
				for (int j = xl; j < xr; j++) {
					float z;
					z = -(A * j + B * i + D) / C;
					GzPut(j, i, ctoi(flatcolor[RED]), ctoi(flatcolor[GREEN]), ctoi(flatcolor[BLUE]), 1, (int)z);
				}
			}
			


			

		}
	}
	return GZ_SUCCESS;
}


int cmp(point x, point y) {
	if (x.y != y.y)
		return x.y < y.y;
	return x.x < y.x;
}


void GzRender::sort_y(point* list) {
	std::sort(list, list + 3, cmp);
}