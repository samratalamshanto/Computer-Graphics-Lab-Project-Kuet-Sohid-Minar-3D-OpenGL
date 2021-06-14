#include<GL/gl.h>
#include<GL/glu.h>
#include<GL/glut.h>
#include<stdlib.h>
#include<stdio.h>
#include<windows.h>
#include<math.h>
#include<iostream>
#include<bits/stdc++.h>
#include "BmpLoader.h"
#include <vector>




using namespace std;

const int wWidth = 1920/2;
const int wHeight = 1080/2;
const float ar = float(wWidth)/float(wHeight);
bool look_point= true;


//main
GLfloat eyeX = 50;
GLfloat eyeY = 120;
GLfloat eyeZ = -10;

GLfloat lookX = 0;
GLfloat lookY = 120;
GLfloat lookZ = 0;

//upper
GLfloat eyeX1 = 100;
GLfloat eyeY1 = 500;
GLfloat eyeZ1 = -100;

GLfloat lookX1 = 0;
GLfloat lookY1 = 120;
GLfloat lookZ1 = 0;




/*
GLfloat eyeX;
GLfloat eyeY;
GLfloat eyeZ;

GLfloat lookX ;
GLfloat lookY ;
GLfloat lookZ ;
*/

GLfloat  axis_x=0.0, axis_y=0.0;
float rotX = 0, rotY = 0, rotZ = 0, theta = 0,theta1 = 0;

bool bRotate= false;
//bool light0 = false, light1 = false, light2 = false, light3= false;
bool light0 = true, light1 = true, light2 = true, light3= true;
bool amb = true, dif = true, spec = true, em = true;
bool Bool_day = true;
int rotate_sky =1;


unsigned int ID[25];
int grass = 0, floor1 = 1, floor2= 2, floor3 = 3, floor4=4, flowers = 5, flowers1 = 6, flowers2 = 7, red = 8,white=9, water=10, alpona=11, black=12,night=13;
int wood=14,leaf=15,building=16,field=17, day=18,pool=19,poster1=20,badminton=21;
float col = 0;

const double PI = 3.14159265389;


/* GLUT callback Handlers */


int anglex= 0, angley = 0, anglez = 0;          //rotation angles
int window;
int wired=0;
int shcpt=1;
int animat = 0;
const int L=5;
const int dgre=3;
int ncpt=L+1;
int clikd=0;
const int nt = 40;				//number of slices along x-direction
const int ntheta = 100;
const int L1=9;
const int dgre1=5;
int ncpt1=L1+1;
const int ntheta1 = 20;


GLfloat ctrlpoints[L+1][3] =
{

    { 0.0, 0.0, 0.0},
    { 0.0, 0.54, 0.0},
    { 0.0, 1.01, 0.0},
    { 0.54, 1.01, 0.0},
    { 0.55, 0.55, 0.0},
    { 0.55, 0.0, 0.0}

};


double ex=0, ey=0, ez=15, lx=0,ly=0,lz=0, hx=0,hy=1,hz=0;

float wcsClkDn[3],wcsClkUp[3];
///////////////////////////////
class point1
{
public:
    point1()
    {
        x=0;
        y=0;
    }
    int x;
    int y;
} clkpt[2];
int flag=0;
GLint viewport[4]; //var to hold the viewport info
GLdouble modelview[16]; //var to hold the modelview info
GLdouble projection[16]; //var to hold the projection matrix info

//////////////////////////
void scsToWcs(float sx,float sy, float wcsv[3] );
void processMouse(int button, int state, int x, int y);
void matColor(float kdr, float kdg, float kdb,  float shiny, int frnt_Back=0, float ambFactor=1.0, float specFactor=1.0);
///////////////////////////

void scsToWcs(float sx,float sy, float wcsv[3] )
{

    GLfloat winX, winY, winZ; //variables to hold screen x,y,z coordinates
    GLdouble worldX, worldY, worldZ; //variables to hold world x,y,z coordinates

    //glGetDoublev( GL_MODELVIEW_MATRIX, modelview ); //get the modelview info
    glGetDoublev( GL_PROJECTION_MATRIX, projection ); //get the projection matrix info
    glGetIntegerv( GL_VIEWPORT, viewport ); //get the viewport info

    winX = sx;
    winY = (float)viewport[3] - (float)sy;
    winZ = 0;

    //get the world coordinates from the screen coordinates
    gluUnProject( winX, winY, winZ, modelview, projection, viewport, &worldX, &worldY, &worldZ);
    wcsv[0]=worldX;
    wcsv[1]=worldY;
    wcsv[2]=worldZ;


}
void processMouse(int button, int state, int x, int y)
{
    if(button==GLUT_LEFT_BUTTON && state==GLUT_DOWN)
    {
        if(flag!=1)
        {
            flag=1;
            clkpt[0].x=x;
            clkpt[0].y=y;
        }


        scsToWcs(clkpt[0].x,clkpt[0].y,wcsClkDn);
        cout<<"\nD: "<<x<<" "<<y<<" wcs: "<<wcsClkDn[0]<<" "<<wcsClkDn[1];
    }
    else if(button==GLUT_LEFT_BUTTON && state==GLUT_UP)
    {
        if (flag==1)
        {
            clkpt[1].x=x;
            clkpt[1].y=y;
            flag=0;
        }
        float wcs[3];
        scsToWcs(clkpt[1].x,clkpt[1].y,wcsClkUp);
        cout<<"\nU: "<<x<<" "<<y<<" wcs: "<<wcsClkUp[0]<<" "<<wcsClkUp[1];

        clikd=!clikd;
    }
}

//control points
long long nCr(int n, int r)
{
    if(r > n / 2) r = n - r; // because C(n, r) == C(n, n - r)
    long long ans = 1;
    int i;

    for(i = 1; i <= r; i++)
    {
        ans *= n - r + i;
        ans /= i;
    }

    return ans;
}


//polynomial interpretation for N points
void BezierCurve ( double t,  float xy[2])
{
    double y=0;
    double x=0;
    t=t>1.0?1.0:t;
    for(int i=0; i<=L; i++)
    {
        int ncr=nCr(L,i);
        double oneMinusTpow=pow(1-t,double(L-i));
        double tPow=pow(t,double(i));
        double coef=oneMinusTpow*tPow*ncr;
        x+=coef*ctrlpoints[i][0];
        y+=coef*ctrlpoints[i][1];

    }
    xy[0] = float(x);
    xy[1] = float(y);

    //return y;
}

///////////////////////
void setNormal(GLfloat x1, GLfloat y1,GLfloat z1, GLfloat x2, GLfloat y2,GLfloat z2, GLfloat x3, GLfloat y3,GLfloat z3)
{
    GLfloat Ux, Uy, Uz, Vx, Vy, Vz, Nx, Ny, Nz;

    Ux = x2-x1;
    Uy = y2-y1;
    Uz = z2-z1;

    Vx = x3-x1;
    Vy = y3-y1;
    Vz = z3-z1;

    Nx = Uy*Vz - Uz*Vy;
    Ny = Uz*Vx - Ux*Vz;
    Nz = Ux*Vy - Uy*Vx;

    glNormal3f(-Nx,-Ny,-Nz);
}

void minar_circle_Bezier()
{
    int i, j;
    float x, y, z, r;				//current coordinates
    float x1, y1, z1, r1;			//next coordinates
    float theta;

    const float startx = 0, endx = ctrlpoints[L][0];
    //number of angular slices
    const float dx = (endx - startx) / nt;	//x step size
    const float dtheta = 2*PI / ntheta;		//angular step size

    float t=0;
    float dt=1.0/nt;
    float xy[2];
    BezierCurve( t,  xy);
    x = xy[0];
    r = xy[1];
    //rotate about z-axis
    float p1x,p1y,p1z,p2x,p2y,p2z;
    for ( i = 0; i < nt; ++i )  			//step through x
    {
        theta = 0;
        t+=dt;
        BezierCurve( t,  xy);
        x1 = xy[0];
        r1 = xy[1];

        //draw the surface composed of quadrilaterals by sweeping theta
        glBegin( GL_QUAD_STRIP );
        for ( j = 0; j <= ntheta; ++j )
        {
            theta += dtheta;
            double cosa = cos( theta );
            double sina = sin ( theta );
            y = r * cosa;
            y1 = r1 * cosa;	//current and next y
            z = r * sina;
            z1 = r1 * sina;	//current and next z

            //edge from point at x to point at next x
            glVertex3f (x, y, z);

            if(j>0)
            {
                setNormal(p1x,p1y,p1z,p2x,p2y,p2z,x, y, z);
            }
            else
            {
                p1x=x;
                p1y=y;
                p1z=z;
                p2x=x1;
                p2y=y1;
                p2z=z1;

            }
            glVertex3f (x1, y1, z1);

            //forms quad with next pair of points with incremented theta value
        }
        glEnd();
        x = x1;
        r = r1;
    } //for i

}
void showControlPoints()
{
    glPointSize(5.0);
    glColor3f(1.0, 0.0, 1.0);
    glBegin(GL_POINTS);
    for (int i = 0; i <=L; i++)
        glVertex3fv(&ctrlpoints[i][0]);
    glEnd();
}


void curve()
{
    const double t = glutGet(GLUT_ELAPSED_TIME) / 5000.0;
    const double a = t*90.0;

    // glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if(wired)
    {
        glPolygonMode( GL_FRONT, GL_LINE ) ;
        glPolygonMode( GL_BACK, GL_LINE ) ;

    }
    else
    {
        glPolygonMode( GL_FRONT,GL_FILL ) ;
        glPolygonMode( GL_BACK, GL_FILL ) ;
    }

    glPushMatrix();

    if(animat)
        glRotated(a,0,0,1);

    glRotatef( anglex, 1.0, 0.0, 0.0);
    glRotatef( angley, 0.0, 1.0, 0.0);         	//rotate about y-axis
    glRotatef( anglez, 0.0, 0.0, 1.0);

    glRotatef( 90, 0.0, 0.0, 1.0);
    glTranslated(-3.5,0,0);
    glGetDoublev( GL_MODELVIEW_MATRIX, modelview ); //get the modelview info

    // void matColor(float kdr, float kdg, float kdb,  float shiny, int frnt_Back=0, float ambFactor=1.0, float specFactor=1.0);

    //  matColor(0.9,0.5,0.1,20);   // front face color
    // matColor(0.0,0.5,0.8,20,1);  // back face color

    glPushMatrix();



    // glRotatef(90,0,1,0);
    // glScalef(0.5,1,0.5);
    minar_circle_Bezier();
    glPopMatrix();






    if(shcpt)
    {
        matColor(0.0,0.0,0.9,20);
        showControlPoints();
    }

    glPopMatrix();

}

static GLfloat v_cube[8][3] =
{
    {0, 0, 0}, {0, 0, 2}, {0, 2, 0}, {0, 2, 2}, {2, 0, 0}, {2, 0, 2}, {2, 2, 0}, {2, 2, 2}
};

static GLbyte c_ind [6][4] =
{
    {0, 2, 6, 4}, {1, 3, 7, 5}, {0, 4, 5, 1}, {2, 6, 7,3}, {0, 1, 3, 2}, {4, 5, 7, 6}
};

static void getNormal3p
(GLfloat x1, GLfloat y1,GLfloat z1, GLfloat x2, GLfloat y2,GLfloat z2, GLfloat x3, GLfloat y3,GLfloat z3)
{
    GLfloat Ux, Uy, Uz, Vx, Vy, Vz, Nx, Ny, Nz;

    Ux = x2-x1;
    Uy = y2-y1;
    Uz = z2-z1;

    Vx = x3-x1;
    Vy = y3-y1;
    Vz = z3-z1;

    Nx = Uy*Vz - Uz*Vy;
    Ny = Uz*Vx - Ux*Vz;
    Nz = Ux*Vy - Uy*Vx;

    glNormal3f(Nx,Ny,Nz);
}


void cube(float R = 255, float G = 255, float B = 255, float alpha = 1)
{
    float colR = R/255, colG = G/255, colB = B/255;
    //glColor3f(1,1,1);
    GLfloat noMat[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat matAmb [] = {colR*.2, colG*.2, colB*.2, 1};
    GLfloat matDif [] = {colR, colG, colB, 1};
    GLfloat matSpec [] = {1, 1, 1, 1};
    GLfloat matEm [] = {0.1, 0.1, 0.1, 1};
    GLfloat matShin[] = {30};

    glMaterialfv(GL_FRONT, GL_AMBIENT, matAmb);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, matDif);
    glMaterialfv(GL_FRONT, GL_SPECULAR, matSpec);
    glMaterialfv(GL_FRONT, GL_SHININESS, matShin);

    if(em)
    {
        glMaterialfv(GL_FRONT, GL_EMISSION, matEm);
    }
    else
    {
        glMaterialfv(GL_FRONT, GL_EMISSION, noMat);
    }



    glBegin(GL_QUADS);
    for(GLint i = 0; i<6; i++)
    {
        getNormal3p(v_cube[c_ind[i][0]][0], v_cube[c_ind[i][0]][1], v_cube[c_ind[i][0]][2],
                    v_cube[c_ind[i][1]][0], v_cube[c_ind[i][1]][1], v_cube[c_ind[i][1]][2],
                    v_cube[c_ind[i][2]][0], v_cube[c_ind[i][2]][1], v_cube[c_ind[i][2]][2]);

        glVertex3fv(&v_cube[c_ind[i][0]][0]);
        glTexCoord2f(1,1);
        glVertex3fv(&v_cube[c_ind[i][1]][0]);
        glTexCoord2f(1,0);
        glVertex3fv(&v_cube[c_ind[i][2]][0]);
        glTexCoord2f(0,0);
        glVertex3fv(&v_cube[c_ind[i][3]][0]);
        glTexCoord2f(0,1);

        for(GLint j =0; j<4; j++)
        {
            glVertex3fv(&v_cube[c_ind[i][j]][0]);
        }
    }
    glEnd();
}


void light()
{
    //Light
    glEnable(GL_LIGHTING);
    GLfloat noLight[] = {0, 0, 0, 1};
    GLfloat lightAmb[] = {0.5, 0.5, 0.5, 1};
    GLfloat lightDif[] = {1, 1, 1, 1};
    GLfloat lightSpec[] = {1, 1, 1, 1};
    GLfloat light1Pos[] = {0, -500, 0, 1};
    GLfloat light4Pos[] = {0,500,0, 1};

    // GLfloat light1Pos[] = {90, 90, 90, 1};
    //  GLfloat light4Pos[] = {90, 90, -90, 1};
    GLfloat light2Pos[] = {683, 300, -350, 1}; //spot light
    GLfloat light3Pos[] = {-380, 400, -50, 1}; //spot light  GLfloat light2Pos[] = {15, 40, -45, 1}; //spot light
    // GLfloat light3Pos[] = {15, 40, 45, 1

    glEnable(GL_LIGHT0);  //1
    glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmb);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDif);
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpec);
    glLightfv(GL_LIGHT0, GL_POSITION, light1Pos);


    glEnable(GL_LIGHT1); //2 spot light
    glLightfv(GL_LIGHT1, GL_AMBIENT, lightAmb);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, lightDif);
    glLightfv(GL_LIGHT1, GL_SPECULAR, lightSpec);
    glLightfv(GL_LIGHT1, GL_POSITION, light2Pos);

    glEnable(GL_LIGHT2); //3 spot light
    glLightfv(GL_LIGHT2, GL_AMBIENT, lightAmb);
    glLightfv(GL_LIGHT2, GL_DIFFUSE, lightDif);
    glLightfv(GL_LIGHT2, GL_SPECULAR, lightSpec);
    glLightfv(GL_LIGHT2, GL_POSITION, light3Pos);

    glEnable(GL_LIGHT3); //4
    glLightfv(GL_LIGHT3, GL_AMBIENT, lightAmb);
    glLightfv(GL_LIGHT3, GL_DIFFUSE, lightDif);
    glLightfv(GL_LIGHT3, GL_SPECULAR, lightSpec);
    glLightfv(GL_LIGHT3, GL_POSITION, light4Pos);

    //GLfloat spotDirection[] = {0, -18, 20, 1};

    //2    {90, 36, -120, 1};
    //3    {125, 36, 180, 1};

    GLfloat spotDirection[] = {0, -1, 0, 1};   //2
    GLfloat spotCutOff[] = {60};

    glLightfv(GL_LIGHT2, GL_SPOT_DIRECTION, spotDirection);
    glLightfv(GL_LIGHT2, GL_SPOT_CUTOFF, spotCutOff);

    GLfloat spotDirection2[] = {0, -1, 0, 1}; //3
    GLfloat spotCutOff2[] = {60};

    glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, spotDirection2);
    glLightfv(GL_LIGHT1, GL_SPOT_CUTOFF, spotCutOff2);


    if(light0)
    {
        glEnable(GL_LIGHT0);
    }
    else
    {
        glDisable(GL_LIGHT0);
    }

    if(light1)
    {
        glEnable(GL_LIGHT1);
    }
    else
    {
        glDisable(GL_LIGHT1);
    }

    if(light2)
    {
        glEnable(GL_LIGHT2);
    }
    else
    {
        glDisable(GL_LIGHT2);
    }
    if(light3)
    {
        glEnable(GL_LIGHT3);
    }
    else
    {
        glDisable(GL_LIGHT3);
    }


    if(amb)
    {
        glLightfv(GL_LIGHT0,GL_AMBIENT,lightAmb);
        glLightfv(GL_LIGHT1,GL_AMBIENT,lightAmb);
        glLightfv(GL_LIGHT2,GL_AMBIENT,lightAmb);
        glLightfv(GL_LIGHT3,GL_AMBIENT,lightAmb);

    }
    else
    {
        glLightfv(GL_LIGHT0,GL_AMBIENT,noLight);
        glLightfv(GL_LIGHT1,GL_AMBIENT,noLight);
        glLightfv(GL_LIGHT2,GL_AMBIENT,noLight);
        glLightfv(GL_LIGHT3,GL_AMBIENT,noLight);
    }

    if(dif)
    {
        glLightfv(GL_LIGHT0,GL_DIFFUSE,lightDif);
        glLightfv(GL_LIGHT1,GL_DIFFUSE,lightDif);
        glLightfv(GL_LIGHT2,GL_DIFFUSE,lightDif);
        glLightfv(GL_LIGHT3,GL_DIFFUSE,lightDif);
    }
    else
    {
        glLightfv(GL_LIGHT0,GL_DIFFUSE,noLight);
        glLightfv(GL_LIGHT1,GL_DIFFUSE,noLight);
        glLightfv(GL_LIGHT2,GL_DIFFUSE,noLight);
        glLightfv(GL_LIGHT3,GL_DIFFUSE,noLight);
    }

    if(spec)
    {
        glLightfv(GL_LIGHT0,GL_SPECULAR,lightSpec);
        glLightfv(GL_LIGHT1,GL_SPECULAR,lightSpec);
        glLightfv(GL_LIGHT2,GL_SPECULAR,lightSpec);
        glLightfv(GL_LIGHT3,GL_SPECULAR,lightSpec);
    }
    else
    {
        glLightfv(GL_LIGHT0,GL_SPECULAR,noLight);
        glLightfv(GL_LIGHT1,GL_SPECULAR,noLight);
        glLightfv(GL_LIGHT2,GL_SPECULAR,noLight);
        glLightfv(GL_LIGHT3,GL_SPECULAR,noLight);
    }

}

void axes()
{
    float length = 30;
    float width = 0.3;

    //X-axis
    glPushMatrix();
    glTranslatef(length/2, 0, 0);
    glScalef(length, width, width);
    glTranslatef(-0.5, -0.5, -0.5);
    cube(255, 0, 0);
    glPopMatrix();

    //Y-axis
    glPushMatrix();
    glTranslatef(0, length/2, 0);
    glScalef(width, length, width);
    glTranslatef(-0.5, -0.5, -0.5);
    cube(0, 255, 0);
    glPopMatrix();

    //Z-axis
    glPushMatrix();
    glTranslatef(0, 0, length/2);
    glScalef(width, width, length);
    glTranslatef(-0.5, -0.5, -0.5);
    cube(0, 0, 255);
    glPopMatrix();
}



static void resize(int width, int height)
{
    //const float ar = (float) width / (float) height;
    glViewport(0, 0, width, width/ar);

}
void cylinder(float w, float h, float l,float r,float g,float b)
{
    int numberOfCube = 360;
    float width = w;
    float height= h;
    float length = l;
    for(double i=0; i<=360; i=i+1.1)
    {
        glPushMatrix();
        glRotated(i, 0, 1, 0);
        glScalef(width,height,length);
        glTranslatef(-0.5, -0.5, -0.5);
        cube(r,g,b);
        glPopMatrix();
    }
}

void main_floor()
{
    glBindTexture(GL_TEXTURE_2D, ID[white]);//main floor
    glPushMatrix();
    glScalef(900,1,1040);
    glTranslatef(-1,-1,-1);
    cube(	139, 69, 19);  //165, 42, 42
    glPopMatrix();

}


void circle_3D(GLdouble radius)
{
    GLUquadric *qobj = gluNewQuadric();
    gluQuadricTexture(qobj, GL_TRUE);

    glRotatef(270, 1, 0, 0);
    gluSphere(qobj, radius, 20, 20);
    gluDeleteQuadric(qobj);

}
//r = r/255.0,g = g/255.0,b = b/255.0;
// glColor3f(r,g,b);
void cylinder_3D(GLdouble height,GLdouble rad,GLdouble rad2)
{
    GLUquadric *qobj = gluNewQuadric();
    gluQuadricTexture(qobj, GL_TRUE);
    glRotatef(90, 1, 0, 0);

    gluCylinder(qobj,  rad, rad2, height, 20, 20);
    gluDeleteQuadric(qobj);

}
void sub_tree()
{
    glBindTexture(GL_TEXTURE_2D, ID[wood]);
    glPushMatrix();
    glTranslatef(0,40,0);
    cylinder_3D(25,1,1);//base
    glPopMatrix();

    glBindTexture(GL_TEXTURE_2D, ID[leaf]);
    glPushMatrix();
    glTranslatef(0,40,0);
    circle_3D(10); //leaf
    glPopMatrix();

}

void sub_tree_upper()
{
    glPushMatrix();
    glTranslatef(0,90,0);
    glRotatef(90,0,1,0);
    glScalef(2,1,2);
    sub_tree();  //1
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,88,0);
    glRotatef(10,1,0,0);
    glScalef(2,1,2);
    sub_tree();  //1
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,88,0);
    glRotatef(-10,1,0,0);
    glScalef(2,1,2);
    sub_tree();  //1
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,90,0);
    glRotatef(10,0,0,1);
    glScalef(2,1,2);
    sub_tree();  //1
    glPopMatrix();


    glPushMatrix();
    glTranslatef(0,90,0);
    glRotatef(-10,0,0,1);
    glScalef(2,1,2);
    sub_tree();  //1
    glPopMatrix();
}
void tree()
{

//11111111111111111111111111111111111111


    glPushMatrix();
    glTranslatef(0,95,-8);
    glRotatef(55,1,0,0);
    glScalef(0.7,0.7,0.7);
    sub_tree();  //1
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,95,8);
    glRotatef(-55,1,0,0);
    glScalef(0.7,0.7,0.7);
    sub_tree();  //1
    glPopMatrix();


//3333333333333333333333333

    glPushMatrix();
    sub_tree_upper();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,10,0);
    glScalef(0.8,1,0.8);
    sub_tree_upper();
    glPopMatrix();
    glBindTexture(GL_TEXTURE_2D, ID[wood]);
    glPushMatrix();
    glTranslatef(0,120,0);
    cylinder_3D(80,4,8);  //tree base
    glPopMatrix();

}

floor_grass()

{

    /* int i,j=0;
     int w=5;
     int h=5;
     for(i=0; i<200; i+w)
     {
         for(j=0; j<200; j+h)
         {
             glBindTexture(GL_TEXTURE_2D, ID[grass]);
             glPushMatrix();
             glTranslatef(i,0,j);
             glScalef(w+i,5,h+j);
             glTranslatef(-1,-1,-1);
             cube(0, 255, 0);
             glPopMatrix();
         }
     }
    */

    glPushMatrix();
    glScalef(250,5,250);
    glTranslatef(-1,-1,-1);
    cube(0, 255, 0);
    glPopMatrix();


}
floor_num1()
{
    glPushMatrix();
    glTranslatef(0,12,0);
    //glTranslatef(-25,1,-20);
    glScalef(100,5,150);
    glTranslatef(-1,-1,-1);
    cube(128, 0, 0);
    glPopMatrix();





}
floor_num2()
{
    glPushMatrix();
    glTranslatef(-3.9,21,0);
    glScalef(96,5,142);
    glTranslatef(-1,-1,-1);
    cube(40,40,40);
    glPopMatrix();
}
floor_num3()
{

    glPushMatrix();
    glTranslatef(-8.8,30,0);
    glScalef(90,5,132);
    glTranslatef(-1,-1,-1);
    cube(128, 10, 10);
    glPopMatrix();


}
pond()
{
    glPushMatrix();
    glTranslatef(-355,1,-30);
    glScalef(200,5,252);
    glTranslatef(-1,-1,-1);
    cube(255, 255, 255);
    glPopMatrix();
}

alponaa()
{
    glPushMatrix();
    glTranslatef(40.8,35,0);

    glScalef(40,0.5,75);
    glTranslatef(-1,-1,-1);
    cube(255,255,255);
    glPopMatrix();

}
floor_set()
{
    glBindTexture(GL_TEXTURE_2D, ID[grass]);
    floor_grass();
    glBindTexture(GL_TEXTURE_2D, ID[water]);
    pond();


    glPushMatrix();   //main_stage
    glScalef(1,1,1.25);
    glBindTexture(GL_TEXTURE_2D, ID[floor4]);
    floor_num1();




    glBindTexture(GL_TEXTURE_2D, ID[floor2]);
    floor_num2();
    glBindTexture(GL_TEXTURE_2D, ID[floor3]);
    floor_num3();
    glPopMatrix();
    glBindTexture(GL_TEXTURE_2D, ID[alpona]);  //alpona
    alponaa();

    /*
    */
}

minar_base()
{
    glBindTexture(GL_TEXTURE_2D, ID[white]);

    glPushMatrix();  //minar base1
    glTranslatef(-70,35,-15);
    glScalef(3,60,3);
    cube(255,255,255);
    glPopMatrix();

    glPushMatrix();//minar base2
    glTranslatef(-70,35,15);
    glScalef(3,60,3);
    cube(255,255,255);
    glPopMatrix();

    glPushMatrix();//minar upper base
    glTranslatef(-70,152,-14.99);
    glScalef(3,2,18);
    cube(255,255,255);
    glPopMatrix();

    //glBindTexture(GL_TEXTURE_2D, ID[black]);
    glPushMatrix();//minar thin_base1
    glTranslatef(-68,35,-2);
    glScalef(0.8,59,0.8);
    cube(255,255,255);
    glPopMatrix();

    //glBindTexture(GL_TEXTURE_2D, ID[black]);
    glPushMatrix();//minar thin_base1
    glTranslatef(-68,35,6);
    glScalef(0.8,59,0.8);
    cube(255,255,255);
    glPopMatrix();
}

mid_minar_base()
{
    glPushMatrix();  //side_mid_base1
    glRotatef(25,0,1,0);
    glTranslatef(0,0,12);
    glScalef(1,.9,.9);
    minar_base();
    glPopMatrix();

    glPushMatrix();  //side_mid_base2
    glRotatef(-25,0,1,0);
    glTranslatef(0,0,-12);
    glScalef(1,.9,.9);
    minar_base();
    glPopMatrix();
}

small_minar_base()
{
    glPushMatrix(); //small_minar_base1
    glRotatef(-25,0,1,0);
    glTranslatef(0,0,-47);
    glScalef(1,.8,.8);
    minar_base();
    glPopMatrix();

    glPushMatrix(); //small_minar_base2
    glRotatef(25,0,1,0);
    glTranslatef(0,0,48);
    glScalef(1,.8,.8);
    minar_base();
    glPopMatrix();



}

upper_minar_base()
{
    glBindTexture(GL_TEXTURE_2D, ID[white]);

    glPushMatrix();  //minar base1
    glTranslatef(-70,35,-15);
    glScalef(3,30,3);
    cube(255,255,255);
    glPopMatrix();

    glPushMatrix();//minar base2
    glTranslatef(-70,35,15);
    glScalef(3,30,3);
    cube(255,255,255);
    glPopMatrix();

    glPushMatrix();//minar upper base
    glTranslatef(-70,92,-14.99);
    glScalef(3,2,18);
    cube(255,255,255);
    glPopMatrix();

    //  glBindTexture(GL_TEXTURE_2D, ID[black]);
    glPushMatrix();//minar thin_base1
    glTranslatef(-68,35,-2);
    glScalef(0.8,29,0.8);
    cube(255,255,255);
    glPopMatrix();

    //glBindTexture(GL_TEXTURE_2D, ID[black]);
    glPushMatrix();//minar thin_base1
    glTranslatef(-68,35,6);
    glScalef(0.8,29,0.8);
    cube(255,255,255);
    glPopMatrix();

}

upper_minar_base_main()
{
    glPushMatrix();
    //glTranslatef(-50,0,0);
    glTranslatef(-44,82,0);
    glRotatef(-45,0,0,1);
    upper_minar_base();
    glPopMatrix();

}

minar_set()
{
    glPushMatrix();
    glScalef(1,1,1.5);
    minar_base();
    glPopMatrix();

    glPushMatrix();
    glScalef(1,1,1.5);
    upper_minar_base_main();
    glPopMatrix();


    glPushMatrix();
    glScalef(1,1,1.5);
    mid_minar_base();
    glPopMatrix();

    glPushMatrix();
    glScalef(1,1,1.5);
    small_minar_base();
    glPopMatrix();


    glEnable(GL_COLOR_MATERIAL);
    glBindTexture(GL_TEXTURE_2D, ID[red]);
    glPushMatrix();
    glColor3f(1,0,0);
    glTranslatef(-130,90,0);
    glRotatef(85,0,0,1);
    glScalef(95,5,95);               //curve-------------------------
    curve();
    glPopMatrix();
    glColor3f(1,1,1);
    glDisable(GL_COLOR_MATERIAL);



}
void road()
{
    glBindTexture(GL_TEXTURE_2D, ID[floor4]);
    glPushMatrix();
    glTranslatef(100,5,-50);
    glScalef(75,1,35);
    cube(165, 42, 42);
    glPopMatrix();

    glBindTexture(GL_TEXTURE_2D, ID[white]);

    glPushMatrix();
    glTranslatef(100,5.5,-55);
    glScalef(75,1,3);
    cube(255,255,255);
    glPopMatrix();


    glPushMatrix();
    glTranslatef(100,5.5,15);
    glScalef(75,1,3);
    cube(255,255,255);
    glPopMatrix();


    glPushMatrix();
    glTranslatef(180,-60,-150);
    glScalef(1.5,1.5,1.5);
    tree();   //sohid minar tree
    glPopMatrix();

    glPushMatrix();
    glTranslatef(180,-60,130);
    glScalef(1.5,1.5,1.5);
    tree();   //sohid minar tree2
    glPopMatrix();


}



void circle(GLfloat x, GLfloat y, GLfloat radiusss)
{
    glBegin(GL_POLYGON);                        // Middle circle
    double radius = radiusss;
    double ori_x = x;                         // the origin or center of circle
    double ori_y = y;
    for (int i = 0; i <= 360; i++)
    {
        double angle = 2 * 3.1416 * i / 360;
        double x = cos(angle) * radius;
        double y = sin(angle) * radius;
        glColor3f(1.0f, 0.0f, 0.0f);
        glVertex2d(ori_x + x, ori_y + y);
    }
    glEnd();
}

void flower_stage()
{
    glPushMatrix();

    glEnable(GL_COLOR_MATERIAL);

    glBindTexture(GL_TEXTURE_2D, ID[flowers2]);  //flowers1
    glPushMatrix();
    glColor3f(1,1,1);
    glTranslatef(-10,50,30);
    glRotatef(-45,0,0,1);
    glScalef(15,1,15);
    curve();
    glPopMatrix();


    glBindTexture(GL_TEXTURE_2D, ID[flowers2]); //flowers2
    glPushMatrix();
    glColor3f(1,1,1);
    glTranslatef(-10,50,-20);
    glRotatef(-45,0,0,1);
    glScalef(15,1,15);
    curve();
    glPopMatrix();
    glColor3f(1,1,1);
    glDisable(GL_COLOR_MATERIAL);

    glEnable(GL_TEXTURE_2D);
    glPushMatrix();    //flower_texture1
    glTranslatef(-8,1,-8);
    glTranslatef(-10,50,30);
    glRotatef(-45,0,0,1);
    glScalef(7.5,1,7.5);
    cube(255,255,255);
    glPopMatrix();

    glPushMatrix();    //flower_texture2

    glTranslatef(-8,1,-8);
    glTranslatef(-10,50,-20);
    glRotatef(-45,0,0,1);
    glScalef(7.5,1,7.5);
    cube(255,255,255);

    glPopMatrix();

    glBindTexture(GL_TEXTURE_2D, ID[floor3]);
    glPushMatrix();//flower stager
    glTranslatef(-40,40,-35);
    glScalef(10,7,40);
    cube(15,0,0);
    glPopMatrix();
    glPopMatrix();
}

void spot_light_position()
{
    glBindTexture(GL_TEXTURE_2D, ID[black]);   //spot_light_tube1
    glPushMatrix();
    glTranslatef(90,0,0);
    glRotatef(-35,0,1,0);
    glTranslatef(50,35,135);
    glScalef(2,5,5);
    glRotatef(-25,0,0,1);
    cube(0,0,0);
    glPopMatrix();

    glBindTexture(GL_TEXTURE_2D, ID[black]);   //spot_light_tube2
    glPushMatrix();
    glTranslatef(120,0,0);
    glRotatef(35,0,1,0);
    glTranslatef(0,30,-165);
    glScalef(2,5,5);
    glRotatef(25,0,0,1);
    cube(0,0,0);
    glPopMatrix();
}

void sohid_minar()
{
    floor_set();
    road();
    glPushMatrix();
    glScalef(1,1.5,1.15);
    glTranslatef(-10,-10,0);
    minar_set();
    glPopMatrix();
    glPushMatrix();
    glTranslatef(-10,0,0);
    flower_stage();
    glPopMatrix();
    spot_light_position();

}



void goal_bar()
{
    glPushMatrix();  //goal bar1
    glTranslatef(0,0,-25);
    glScalef(2,25,2);
    cube(255,255,255);
    glPopMatrix();

    glPushMatrix();//goal bar2
    glTranslatef(0,0,25);
    glScalef(2,25,2);
    cube(255,255,255);
    glPopMatrix();

    glPushMatrix();//goal mid bar
    glTranslatef(0,50,-24.25);
    glScalef(2,2,26.5);
    cube(255,255,255);
    glPopMatrix();

}

void debox_line()
{
    glPushMatrix(); //bar line2
    glTranslatef(-180,0,0);
    glTranslatef(0,6,0);
    glScalef(4,5,55);
    glTranslatef(-1,-1,-1);
    cube(255, 255, 255);
    glPopMatrix();


    glPushMatrix(); //side line
    glTranslatef(-215,0,-50);
    glTranslatef(0,6,0);
    glScalef(35,5,4);
    glTranslatef(-1,-1,-1);
    cube(255, 255, 255);
    glPopMatrix();

    glPushMatrix(); //side line
    glTranslatef(-215,0,50);
    glTranslatef(0,6,0);
    glScalef(35,5,4);
    glTranslatef(-1,-1,-1);
    cube(255, 255, 255);
    glPopMatrix();
}

void gallery()
{
    glBindTexture(GL_TEXTURE_2D, ID[white]);//football field
    glPushMatrix(); //gallery line
    glTranslatef(0,0,-20);
    glTranslatef(-10,0,-210);
    glTranslatef(0,6,0);
    glScalef(220,25,30);
    glTranslatef(-1,-1,-1);
    cube(255, 255, 255);
    glPopMatrix();

    glPushMatrix(); // gallery line
    glTranslatef(0,0,-20);
    glTranslatef(0,0,-20);
    glTranslatef(-10,25,-210);
    glTranslatef(0,6,0);
    glScalef(220,25,30);
    glTranslatef(-1,-1,-1);
    cube(255, 255, 255);
    glPopMatrix();


    glPushMatrix(); //gallery line2
    glTranslatef(0,0,20);
    glTranslatef(-10,0,210);
    glTranslatef(0,6,0);
    glScalef(220,25,30);
    glTranslatef(-1,-1,-1);
    cube(255, 255, 255);
    glPopMatrix();

    glPushMatrix(); // gallery line2
    glTranslatef(0,0,20);
    glTranslatef(0,0,20);
    glTranslatef(-10,25,210);
    glTranslatef(0,6,0);
    glScalef(220,25,30);
    glTranslatef(-1,-1,-1);
    cube(255, 255, 255);
    glPopMatrix();


//treeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
    glPushMatrix();
    glTranslatef(270,-60,190);
    glScalef(1.5,1.5,1.5);
    tree();   //field tree
    glPopMatrix();

    glPushMatrix();
    glTranslatef(291,-60,-150);
    glScalef(1.5,1.5,1.5);
    tree();   //field tree
    glPopMatrix();



    glPushMatrix();
    glTranslatef(291,-60,30);
    glScalef(1.5,1.5,1.5);
    tree();   //field tree
    glPopMatrix();
//downs side tree
    glPushMatrix();
    glTranslatef(-270,-60,210);
    glScalef(1.5,1.5,1.5);
    tree();   //field tree
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-290,-60,100);
    glScalef(1.5,1.5,1.5);
    tree();   //field tree
    glPopMatrix();


    glPushMatrix();
    glTranslatef(-290,-60,-70);
    glScalef(1.5,1.5,1.5);
    tree();   //field tree
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-291,-60,-150);
    glScalef(1.5,1.5,1.5);
    tree();   //field tree
    glPopMatrix();

}
void football_field()
{

    glBindTexture(GL_TEXTURE_2D, ID[grass]);//football field

    glPushMatrix();
    glTranslatef(0,5,0);
    glScalef(400,5,280);
    glTranslatef(-1,-1,-1);
    cube(0, 255, 0);
    glPopMatrix();




    glBindTexture(GL_TEXTURE_2D, ID[field]);//football field
    glPushMatrix(); //field
    glRotatef(90,0,1,0);
    glTranslatef(0,6,0);
    glScalef(200,5,250);
    glTranslatef(-1,-1,-1);
    cube(255, 255, 255);
    glPopMatrix();

    glBindTexture(GL_TEXTURE_2D, ID[white]);//football field
    glPushMatrix();
    glTranslatef(-235,0,0);  //bar1
    goal_bar();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(235,0,0); //bar2
    goal_bar();
    glPopMatrix();



    glPushMatrix();

    gallery();
    glPopMatrix();


}

void hall_gate_stand()
{
    glPushMatrix();//middle stand of the hall
    glTranslatef(250,5,110);
    glScalef(4,25,4);
    cube(255,255,255);
    glPopMatrix();

    glPushMatrix();//middle stand of the hall
    glTranslatef(322,5,110);
    glScalef(4,25,4);
    cube(255,255,255);
    glPopMatrix();

}

void hall_building()
{

    glPushMatrix();
    glBindTexture(GL_TEXTURE_2D, ID[white]);
    glTranslatef(0,239,0);
    glScalef(90,1,90); //building1_top
    cube(105, 105, 105);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,239,0);
    glTranslatef(400,0,00); //building2_top
    glScalef(90,1,90);
    cube(105, 105, 105);
    glPopMatrix();

    glBindTexture(GL_TEXTURE_2D, ID[building]);
    glPushMatrix();
    glScalef(90,120,90); //building1
    cube(255,255,255);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(400,0,00); //building2
    glScalef(90,120,90);
    cube(255,255,255);
    glPopMatrix();

    /*
    */
    glBindTexture(GL_TEXTURE_2D, ID[white]);
    glPushMatrix();   //middle floor 1 of the hall
    glTranslatef(150,5,30);
    glScalef(140,1,40);
    cube(	255,25,02);
    glPopMatrix();


    glPushMatrix();//middle upper floor 2 of the hall
    glTranslatef(150,55,30);
    glScalef(140,1,40);
    cube(	70, 130, 180);
    glPopMatrix();


    glPushMatrix();//middle road of the hall
    glTranslatef(250,5,110);
    glScalef(40,1,140);
    cube(255,25,0);
    glPopMatrix();


    glPushMatrix();//middle upper road of the hall
    glTranslatef(250,55,110);
    glScalef(40,1,140);
    cube(	70, 130, 180);
    glPopMatrix();

    glBindTexture(GL_TEXTURE_2D, ID[white]);
    hall_gate_stand(); //hall road stand 1


    glPushMatrix(); //hall road stand 2
    glTranslatef(0,0,90);
    hall_gate_stand();
    glPopMatrix();

    glPushMatrix(); //hall road stand 3
    glTranslatef(0,0,180);
    hall_gate_stand();
    glPopMatrix();


    glPushMatrix(); //hall road stand 4
    glTranslatef(0,0,270);
    hall_gate_stand();
    glPopMatrix();


    glPushMatrix();   //middle floor stand 1
    glTranslatef(80,0,360);
    glRotatef(90,0,1,0);
    hall_gate_stand();
    glPopMatrix();


    glPushMatrix();   //middle floor stand 1
    glTranslatef(250,0,360);
    glRotatef(90,0,1,0);
    hall_gate_stand();
    glPopMatrix();


}

void hall()
{


    glPushMatrix(); //hall tree
    glTranslatef(40,-80,300);
    glScalef(1.5,1.5,1.5);
    tree();
    glPopMatrix();
    glPushMatrix(); //hall tree
    glTranslatef(40,-80,360);
    glScalef(1.5,1.5,1.5);
    tree();
    glPopMatrix();



    glPushMatrix(); //hall tree
    glTranslatef(90,-80,300);
    glScalef(1.5,1.5,1.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //hall tree2
    glTranslatef(420,-80,350);
    glScalef(1.5,1.5,1.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //hall tree2
    glTranslatef(460,-80,300);
    glScalef(1.5,1.5,1.5);
    tree();
    glPopMatrix();

    hall_building();

    glPushMatrix();
    glBindTexture(GL_TEXTURE_2D, ID[grass]);
    glTranslatef(280,0,180);
    glScalef(350,6,470); //floor of the hall
    glTranslatef(-1,-1,-1);
    cube(0,255,0);
    glPopMatrix();




}
void kuet()
{

    //KUET

    //K

    glBindTexture(GL_TEXTURE_2D, ID[white]);

    glPushMatrix();  //k1
    glTranslatef(0,30,0);
    glScalef(4,45,4);
    cube(255,255,255);
    glPopMatrix();

    glPushMatrix(); //k2
    glTranslatef(0,75,0);
    glRotatef(-45,1,0,0);
    glScalef(4,25,4);
    cube(255,255,255);
    glPopMatrix();

    glPushMatrix(); //k3
    glTranslatef(0,75,5);
    glRotatef(-145,1,0,0);
    glScalef(4,25,4);
    cube(255,255,255);
    glPopMatrix();




    //U

    glPushMatrix(); //u1
    glTranslatef(0,30,-50);
    glScalef(4,42,4);
    cube(255,255,255);
    glPopMatrix();

    glPushMatrix(); //u2
    glTranslatef(0,30,-80);
    glScalef(4,42,4);
    cube(255,255,255);
    glPopMatrix();

    glPushMatrix(); //u3
    glTranslatef(0,30,-80.5);
    glScalef(4,4,18);
    cube(255,255,255);
    glPopMatrix();



    //E

    glPushMatrix(); //E1
    glTranslatef(0,30,-130);
    glScalef(4,4,20);
    cube(255,255,255);
    glPopMatrix();


    glPushMatrix();//E2
    glTranslatef(0,30,-98);
    glScalef(4,42,4);
    cube(255,255,255);
    glPopMatrix();

    glPushMatrix();//E3
    glTranslatef(0,69,-130);
    glScalef(4,4,20);
    cube(255,255,255);
    glPopMatrix();


    glPushMatrix();//E4
    glTranslatef(0,108,-130);
    glScalef(4,4,20);
    cube(255,255,255);
    glPopMatrix();



    //T

    glPushMatrix();//T1
    glTranslatef(0,108,-190);
    glScalef(4,4,25);
    cube(255,255,255);
    glPopMatrix();


    glPushMatrix();//T2
    glTranslatef(0,30,-170);
    glScalef(4,42,4);
    cube(255,255,255);
    glPopMatrix();


}
void minar_field()
{
    sohid_minar();
    glPushMatrix();
    glTranslatef(-150,-2,500);
    football_field();
    glPopMatrix();
}

void main_roads1()
{
    glBindTexture(GL_TEXTURE_2D, ID[black]);
    glPushMatrix();
    glTranslatef(100,5,-50);
    glScalef(585,1,90);
    cube(165, 42, 42);
    glPopMatrix();

    glBindTexture(GL_TEXTURE_2D, ID[white]);

    glPushMatrix();
    glTranslatef(100,5.5,-55);
    glScalef(585,1,3);
    cube(255,255,255);
    glPopMatrix();


    glPushMatrix();
    glTranslatef(100,5.5,125);
    glScalef(630,1,3);
    cube(255,255,255);
    glPopMatrix();

//treeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
    glPushMatrix(); //roads1 tree
    glTranslatef(150,-55.5,-100);
    glScalef(1.5,1.5,1.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads1 tree
    glTranslatef(290,-55.5,-100);
    glScalef(1.5,1.5,1.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads1 tree
    glTranslatef(330,-55.5,-100);
    glScalef(1.5,1.5,1.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads1 tree
    glTranslatef(460,-55.5,-100);
    glScalef(1.5,1.5,1.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads1 tree
    glTranslatef(520,-55.5,-100);
    glScalef(1.5,1.5,1.5);
    tree();
    glPopMatrix();



    glPushMatrix(); //roads1 tree
    glTranslatef(690,-55.5,-100);
    glScalef(1.5,1.5,1.5);
    tree();
    glPopMatrix();


    glPushMatrix(); //roads1 tree
    glTranslatef(860,-55.5,-100);
    glScalef(1.5,1.5,1.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads1 tree
    glTranslatef(920,-55.5,-100);
    glScalef(1.5,1.5,1.5);
    tree();
    glPopMatrix();
//-------------------------------------------------------
    glPushMatrix(); //roads1 tree
    glTranslatef(150,-55.5,170);
    glScalef(1.5,1.5,1.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads1 tree
    glTranslatef(290,-55.5,170);
    glScalef(1.5,1.5,1.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads1 tree
    glTranslatef(330,-55.5,170);
    glScalef(1.5,1.5,1.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads1 tree
    glTranslatef(460,-55.5,170);
    glScalef(1.5,1.5,1.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads1 tree
    glTranslatef(520,-55.5,170);
    glScalef(1.5,1.5,1.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads1 tree
    glTranslatef(550,-55.5,170);
    glScalef(1.5,1.5,1.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads1 tree
    glTranslatef(690,-55.5,170);
    glScalef(1.5,1.5,1.5);
    tree();
    glPopMatrix();


    glPushMatrix(); //roads1 tree
    glTranslatef(1160,-55.5,170);
    glScalef(2.5,1.5,2.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads1 tree
    glTranslatef(1120,-55.5,170);
    glScalef(2.5,1.5,2.5);
    tree();
    glPopMatrix();


}

void main_roads2()
{
    glBindTexture(GL_TEXTURE_2D, ID[black]);
    glPushMatrix();
    glTranslatef(100,5,-50);
    glScalef(505,1,90);
    cube(165, 42, 42);
    glPopMatrix();

    glBindTexture(GL_TEXTURE_2D, ID[white]);

    glPushMatrix();
    glTranslatef(100,5.5,-55);
    glScalef(418,1,3);
    cube(255,255,255);
    glPopMatrix();


    glPushMatrix();
    glTranslatef(100,5.5,125);
    glScalef(415,1,3);
    cube(255,255,255);
    glPopMatrix();



    glPushMatrix(); //roads2 tree
    glTranslatef(400,-55.5,-85);
    glScalef(1.5,1.5,1.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads2 tree
    glTranslatef(720,-55.5,-95);
    glScalef(2.5,1.5,2.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads2 tree
    glTranslatef(580,-55.5,-85);
    glScalef(2.5,1.5,2.5);
    tree();
    glPopMatrix();



}
void main_roads3()
{
    glBindTexture(GL_TEXTURE_2D, ID[black]);
    glPushMatrix();
    glTranslatef(100,5,-50);
    glScalef(505,1,90);
    cube(165, 42, 42);
    glPopMatrix();

    glBindTexture(GL_TEXTURE_2D, ID[white]);

    glPushMatrix();
    glTranslatef(100,5.5,-55);
    glScalef(505,1,3);
    cube(255,255,255);
    glPopMatrix();


    glPushMatrix();
    glTranslatef(56,0,0);
    glTranslatef(100,5.5,125);
    glScalef(480,1,3);
    cube(255,255,255);
    glPopMatrix();


    glPushMatrix(); //roads3 tree
    glTranslatef(180,-55.5,-95);
    glScalef(2.5,1.5,2.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads3 tree
    glTranslatef(280,-55.5,-95);
    glScalef(2.5,1.5,2.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads3 tree
    glTranslatef(360,-55.5,-95);
    glScalef(2.5,1.5,2.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads3 tree
    glTranslatef(780,-55.5,-95);
    glScalef(2.5,1.5,2.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads3 tree
    glTranslatef(980,-55.5,-95);
    glScalef(2.5,1.5,2.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads3 tree
    glTranslatef(680,-55.5,-95);
    glScalef(2.5,1.5,2.5);
    tree();
    glPopMatrix();
//---------------------------


    glPushMatrix(); //roads3 tree
    glTranslatef(190,-55.5,185);
    glScalef(2.5,1.5,2.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads3 tree
    glTranslatef(295,-55.5,185);
    glScalef(2.5,1.5,2.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads3 tree
    glTranslatef(390,-55.5,185);
    glScalef(2.5,1.5,2.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads3 tree
    glTranslatef(880,-55.5,185);
    glScalef(2.5,1.5,2.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads3 tree
    glTranslatef(980,-55.5,185);
    glScalef(2.5,1.5,2.5);
    tree();
    glPopMatrix();

    glPushMatrix(); //roads3 tree
    glTranslatef(780,-55.5,185);
    glScalef(2.5,1.5,2.5);
    tree();
    glPopMatrix();

}

void roads()
{
    glPushMatrix(); //main_roads1
    glTranslatef(-22,0,1120);
    glRotatef(90,0,1,0);
    main_roads1();
    glPopMatrix();

    glPushMatrix(); //main_roads2
    glTranslatef(-1000,0,-280);
    //glRotatef(90,0,1,0);
    main_roads2();
    glPopMatrix();

    glPushMatrix(); //main_roads3
    glTranslatef(-100,0,-220);
    glRotatef(45,0,1,0);
    main_roads3();
    glPopMatrix();


}

void pool_bed()
{
    glBindTexture(GL_TEXTURE_2D, ID[white]);
    glPushMatrix();
    glTranslatef(0,30,0);
    glPushMatrix();
    glScalef(35,2,20);
    cube(	30, 144, 255);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(5,2,0);
    glRotatef(135,0,0,1);
    glScalef(15,2,20);
    cube(	30, 144, 255);
    glPopMatrix();


    glPushMatrix(); //leg1
    glTranslatef(0,-19,0);
    glScalef(2,10,20);
    cube(0,0,0);
    glPopMatrix();

    glPushMatrix(); //leg2
    glTranslatef(65,-19,0);
    glScalef(2,10,20);
    cube(0,0,0);
    glPopMatrix();


    glPopMatrix();
}
void ladder()
{
    glBindTexture(GL_TEXTURE_2D, ID[white]);
    glPushMatrix();

    glTranslatef(0,20,0);
    glPushMatrix();
    glScalef(2,60,2);
    cube(0,0,0);
    glPopMatrix();


    glPushMatrix();
    glTranslatef(0,0,30);
    glScalef(2,60,2);
    cube(0,0,0);
    glPopMatrix();


    for(int i =5; i<120; i=i+15)
    {
        glPushMatrix();
        glTranslatef(0,i,0);
        glScalef(4,2,15);
        cube(0,0,0);
        glPopMatrix();
    }

    glPopMatrix();


}
void pool_jumper()
{
    glBindTexture(GL_TEXTURE_2D, ID[white]);

    glPushMatrix();

    glPushMatrix(); //swiming pool jumper
    glColor3f(105/255.0, 105/255.0, 105/255.0);
    glScalef(40,50,40);
    cube(105, 105, 105);
    glPopMatrix();






    glPushMatrix(); //swiming pool jumper stand_base
    glTranslatef(0,120,0);
    glScalef(40,2,2);
    cube(255, 0, 0);
    glPopMatrix();



    glPushMatrix(); //swiming pool jumper height2
    glTranslatef(0,100,76);
    glScalef(2,10,2);
    cube(240, 248, 255);
    glPopMatrix();


    glPushMatrix(); //swiming pool jumper stand3
    glTranslatef(76,100,76);
    glScalef(2,10,2);
    cube(240, 248, 255);
    glPopMatrix();




    glPushMatrix(); //swiming pool jumper board
    glTranslatef(-120,100,10);
    glScalef(80,2,30);
    cube(255,255,255);
    glPopMatrix();



    glPushMatrix(); //swiming pool jumper stand3
    glTranslatef(76,100,0);
    glScalef(2,10,2);
    cube(240, 248, 255);
    glPopMatrix();

    glPushMatrix(); //swiming pool jumper stand2
    glTranslatef(0,100,0);
    glScalef(2,10,2);
    cube(240, 248, 255);
    glPopMatrix();

    glPushMatrix(); //swiming pool jumper stand1
    glTranslatef(0,120,76);
    glScalef(40,2,2);
    cube(255,0,0);
    glPopMatrix();


    glPushMatrix();
    glTranslatef(175,0,20);
    glRotatef(45,0,0,1);
    ladder();                       //ladder
    glPopMatrix();


    glPopMatrix();




}
void umbrella()
{

    glPushMatrix();  //umbrela
    glTranslatef(10,0,0);
    glColor3f(1,1,1);
    glColor3f(0.541, 0.169, 0.886);
    glEnable(GL_COLOR_MATERIAL);
    glPushMatrix();
    glTranslatef(0,80,0);
    glScalef(40,1,40);
    curve();
    glPopMatrix();
    glColor3f(1,1,1);
    glDisable(GL_COLOR_MATERIAL);
    glBindTexture(GL_TEXTURE_2D, ID[white]);


    glPushMatrix();
    glScalef(2,40,2);  //leg
    cube(0, 0,0);
    glPopMatrix();


    glPopMatrix();

}
void swimming_pool()
{
    glBindTexture(GL_TEXTURE_2D, ID[pool]);
    glPushMatrix(); //swiming pool
    glTranslatef(50,0,-70);
    glTranslatef(-910,0,-925);
    glScalef(410,1,295);
    cube(255,255,255);
    glPopMatrix();

    glBindTexture(GL_TEXTURE_2D, ID[white]);
    glPushMatrix();
    glTranslatef(80,0,-700);  //bed1
    glRotatef(180,0,1,0);
    pool_bed();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(70,0,-750);
    glScalef(1.5,1.5,1.5);
    umbrella();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(80,0,-790);  //bed2
    glRotatef(180,0,1,0);
    pool_bed();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(80,0,-840);  //bed3
    glRotatef(180,0,1,0);
    pool_bed();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(70,0,-690);
    glScalef(1.5,1.5,1.5);
    umbrella();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-30,0,-980);  //jumper stand
    pool_jumper();
    glPopMatrix();

}

void kuet_logo()
{
    glPushMatrix();
    glTranslated(0,120,-85);
    glRotatef(180,0,1,0);
    kuet();
    glPopMatrix();

    glBindTexture(GL_TEXTURE_2D, ID[grass]);
    glPushMatrix();
    glScalef(0.8,1,1.5);
    circle_3D(150);
    glPopMatrix();
}

void sky()
{
    glRotatef(theta1, 1,0,0);
    glPushMatrix();
    circle_3D(1300);
    glPopMatrix();

}

void poster()
{
    glBindTexture(GL_TEXTURE_2D, ID[poster1]);////////
    glPushMatrix();
    glTranslatef(0,60,0);
    glScalef(300,5,300);
    cube(255,255,255);
    glPopMatrix();
}
void poster_full()
{
    glPushMatrix();
    glTranslatef(0,5,18);
    glTranslatef(0,350,0);
    glRotatef(90,1,0,0);
    glScalef(0.3,0.3,0.3);
    poster();
    glPopMatrix();

    glBindTexture(GL_TEXTURE_2D, ID[white]);

    glPushMatrix();
    glTranslatef(5,30,30);
    glScalef(5,80,5);
    cube(0,0,255);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(165,30,30);
    glScalef(5,80,5);
    cube(0,0,255);
    glPopMatrix();
}
float rot;
static void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-3, 3, -3, 3, 2.0, 2500.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    //glViewport(0, 0, width, height);

        if(look_point==false)
        {
//            GLfloat eyeX1 = 50;
//            GLfloat eyeY = 620;
//            GLfloat eyeZ = -50;
//            GLfloat lookX = 0;
//            GLfloat lookY = 120;
//            GLfloat lookZ = 0;
            gluLookAt(eyeX1, eyeY1, eyeZ1, lookX1, lookY1, lookZ1, 0, 1, 0);
        }
        else
        {
//            GLfloat eyeX = 50;
//            GLfloat eyeY = 120;
//            GLfloat eyeZ = -50;
//            GLfloat lookX = 0;
//            GLfloat lookY = 120;
//            GLfloat lookZ = 0;
            gluLookAt( eyeX, eyeY, eyeZ,lookX, lookY, lookZ, 0, 1, 0);  //eyeX, eyeY, eyeZ,
        }



   // gluLookAt(eyeX, eyeY, eyeZ, lookX, lookY, lookZ, 0, 1, 0);  //
    glRotatef(theta, axis_x,axis_y,0);

    // glClearColor(65, 105, 225,1); // clear

    rot=2.3;
    glEnable(GL_TEXTURE_2D);
    /*
    glBindTexture(GL_TEXTURE_2D, ID[alpona]);

    */





    glClearColor(1,1,1,1); // clear
    // glBindTexture(GL_TEXTURE_2D, ID[]);





    /*

    glColor3f(1,1,1);
    //glEnable(GL_COLOR_MATERIAL);
    glBindTexture(GL_TEXTURE_2D, ID[alpona]);
    glPushMatrix();
    glTranslatef(0,50,0);

    glScalef(40,10,40);
    curve();
    glPopMatrix();
    glColor3f(1,1,1);
    //glDisable(GL_COLOR_MATERIAL);

    /*




        /*

            glPushMatrix();
            glTranslatef(600,0,-100);
            glRotatef(-90,0,1,0);
            hall();
            glPopMatrix();
*/












    //main elements
    main_floor();




    glPushMatrix();
    glTranslatef(-330,0,100);
    minar_field();
    glPopMatrix();

    roads();

    glPushMatrix();
    glTranslatef(600,0,-100);
    glRotatef(-90,0,1,0);
    hall();
    glPopMatrix();

    glPushMatrix();
    swimming_pool();
    glPopMatrix();

    glPushMatrix();
    glTranslated(670,-50,-420);
    glRotatef(-25,0,1,0);
    kuet_logo();
    glPopMatrix();


    if(Bool_day==true)
    {
        glBindTexture(GL_TEXTURE_2D, ID[day]);
        glPushMatrix();
        sky();
        glPopMatrix();
    }
    else
    {
        glBindTexture(GL_TEXTURE_2D, ID[night]);
        glPushMatrix();
        sky();
        glPopMatrix();

    }

    glPushMatrix();
    glTranslated(-105,-60,405);
    glScalef(0.6,0.5,0.6);
    glRotatef(90,0,1,0);
    poster_full();           //poster1
    glPopMatrix();

    glPushMatrix();
    glTranslated(-185,-60,-395);
    glScalef(0.6,0.5,0.6);
    poster_full();    //poster2
    glPopMatrix();

    glBindTexture(GL_TEXTURE_2D, ID[badminton]);
    glPushMatrix();
    glTranslated(165,5,530);
    glScalef(350,1,230);
    cube(	50, 205, 50); // hall field
    glPopMatrix();




    /*
    glBindTexture(GL_TEXTURE_2D, ID[alpona]);
    glPushMatrix();

    glPopMatrix();
    */








    /*--------------------------------------------------



            /*


                */


    glDisable(GL_TEXTURE_2D);
    // axes();
    light();

    glFlush();
    glutSwapBuffers();
}
void LoadTexture(const char*filename, int index)
{
    glGenTextures(1, &ID[index]);
    glBindTexture(GL_TEXTURE_2D, ID[index]);
    glPixelStorei(GL_UNPACK_ALIGNMENT, ID[index]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    BmpLoader bl(filename);
    gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGB, bl.iWidth, bl.iHeight, GL_RGB, GL_UNSIGNED_BYTE, bl.textureData );
}

void texture()
{



    LoadTexture("C:\\Users\\samra\\Downloads\\bmp\\pitch2.bmp", badminton);

    LoadTexture("C:\\Users\\samra\\Downloads\\bmp\\kuet_board.bmp", poster1);


    LoadTexture("C:\\Users\\samra\\Downloads\\bmp\\sky2.bmp", day);
    LoadTexture("C:\\Users\\samra\\Downloads\\bmp\\night2.bmp", night);

    LoadTexture("C:\\Users\\samra\\Downloads\\bmp\\grass.bmp", grass);



    LoadTexture("C:\\Users\\samra\\Downloads\\bmp\\white2.bmp", white);
    LoadTexture("C:\\Users\\samra\\Downloads\\bmp\\water.bmp", water);

    LoadTexture("C:\\Users\\samra\\Downloads\\bmp\\brick.bmp", floor1);

    LoadTexture("C:\\Users\\samra\\Downloads\\bmp\\brick4.bmp", floor2);
    // LoadTexture("C:\\Users\\samra\\Downloads\\bmp\\brick2.bmp", floor3);
    // LoadTexture("C:\\Users\\samra\\Downloads\\bmp\\brick6.bmp", floor4);

    LoadTexture("C:\\Users\\samra\\Downloads\\bmp\\water.bmp", water);
    LoadTexture("C:\\Users\\samra\\Downloads\\bmp\\alpona4.bmp", alpona);
    LoadTexture("C:\\Users\\samra\\Downloads\\bmp\\red.bmp", red); //
    //  LoadTexture("C:\\Users\\samra\\Downloads\\bmp\\flowers3.bmp", flowers1);
    LoadTexture("C:\\Users\\samra\\Downloads\\bmp\\flowers4.bmp", flowers2);
    LoadTexture("C:\\Users\\samra\\Downloads\\bmp\\black.bmp", black);
    LoadTexture("C:\\Users\\samra\\Downloads\\bmp\\building.bmp", building);
    LoadTexture("C:\\Users\\samra\\Downloads\\bmp\\tree_base2.bmp", wood);
    LoadTexture("C:\\Users\\samra\\Downloads\\bmp\\leaf2.bmp", leaf);

    LoadTexture("C:\\Users\\samra\\Downloads\\bmp\\pool2.bmp",pool);
    LoadTexture("C:\\Users\\samra\\Downloads\\bmp\\field2.bmp", field);



    //
    /*
     /*     /*
     /*
                      */


}

static void keyBoard(unsigned char key, int x, int y)
{
    float x1, z1, r, theta, dx, dz, dx_norm, dz_norm, r1=1, turn_angle_step=10, height_diff_one_less, height_diff_thresh_dist;

    x1=lookX-eyeX;
    z1=lookZ-eyeZ;
    r=sqrt(x1*x1+z1*z1);

    if(x1==0)
    {
        if(z1>0)
        {
            theta = 90;
        }
        else if(z1<0)
        {
            theta = -90;
        }
    }
    else
        theta=atan(z1/x1) * 180 / 3.1416;

    if((z1>0 && theta<0) || (z1<0 && theta>0))
        theta += 180;
    else if(z1<0 && theta<0)
    {
        theta += 360;
    }

    dx = r1*cos(theta * 3.1416 / 180);
    dz = r1*sin(theta * 3.1416 / 180);

    dx_norm = r1*cos((theta-90) * 3.1416 / 180);
    dz_norm = r1*sin((theta-90) * 3.1416 / 180);


    switch (key)
    {
    case 'R':
    case 'r':
        bRotate = !bRotate;
        axis_x=0.0;
        axis_y=1.0;
        break;


    // move look point
    case 'a':
        theta-=turn_angle_step;
        theta = theta * 3.1416 / 180;
        lookX=r*cos(theta)+eyeX;
        lookZ=r*sin(theta)+eyeZ;
        break;
    case 'w':
        lookY=lookY+3;
        break;
    case 's':
        lookY=lookY-3;
        break;
    case 'd':
        theta+=turn_angle_step;
        theta = theta * 3.1416 / 180;
        lookX=r*cos(theta)+eyeX;
        lookZ=r*sin(theta)+eyeZ;
        break;

    // Moving the camera front-back-left-right
    case 'j':
        eyeX += dx_norm*4;
        eyeZ += dz_norm*4;


        lookX += dx_norm*4;
        lookZ += dz_norm*4;
        break;
    case 'i':
        eyeX += dx*4;
        eyeZ += dz*4;
        lookX += dx*4;
        lookZ += dz*4;
        break;
    case 'k':
        eyeX -= dx*4;
        eyeZ -= dz*4;

        lookX -= dx*4;
        lookZ -= dz*4;
        break;
    case 'l':
        eyeX -= dx_norm*4;
        eyeZ -= dz_norm*4;

        lookX -= dx_norm*4;
        lookZ -= dz_norm*4;
        break;

    case '5':
        eyeY += 0.5;
        lookY +=0.5;
        break;
    case '6':
        eyeY -= 0.5;
        lookY -=0.5;
        break;
    case '\'':
        rotY++;
        if(rotY>360)
        {
            rotY = rotY - 360;
        }
        break;
    case '7':
        rotY--;
        if(rotY < 0)
        {
            rotY = rotY + 360;
        }
        break;



    case '1':
        light0 = !light0;
        break;
    case '2':
        light1 = !light1;
        break;
    case '3':
        light2 = !light2;
        break;
    case '4':
        light3 = !light3;
        break;
    case 'z':
        amb = !amb;
        break;
    case 'c':
        dif = !dif;
        break;
    case 'x':
        spec = !spec;
        break;
    case 'v':
        em = !em;
        break;
    case 'q':
        look_point=!look_point;
        /*
         GLfloat eyeX = 50;
         GLfloat eyeY = 620;
         GLfloat eyeZ = 50;
         */
        break;
    }
    glutPostRedisplay();
}

/*
static void key(unsigned char key, int x, int y)
{
    switch (key)
    {
    case 'A':
        animat=!animat;
        break;

    case 's':
    case 'S':
        shcpt=!shcpt;
        break;

    case 'w':
    case 'W':
        wired=!wired;
        break;

    case 'x':
        anglex = ( anglex + 3 ) % 360;
        break;
    case 'X':
        anglex = ( anglex - 3 ) % 360;
        break;

    case 'y':
        angley = ( angley + 3 ) % 360;
        break;
    case 'Y':
        angley = ( angley - 3 ) % 360;
        break;

    case 'z':
        anglez = ( anglez + 3 ) % 360;
        break;
    case 'Z':
        anglez = ( anglez - 3 ) % 360;
        break;


    case 'q':
    case 27 :
        glutDestroyWindow(window);
        exit(0);
        break;
    }

    glutPostRedisplay();
}
*/

bool flag_dir= true;

static void idle(void)
{
    if (bRotate == true)
    {
        theta += rot;
        if(theta > 360.0)
            theta -= 360.0*floor(theta/360.0);
    }

    if (rotate_sky)
    {
        theta1 += 2;
        if( theta1>360  )
        {
            theta1=0;
            Bool_day=!Bool_day;

        }

    }

    glutPostRedisplay();
}






const GLfloat light_ambient[]  = { 0.5f, 0.5f, 0.5f, 1.0f };
const GLfloat light_diffuse[]  = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_position[] = { 2.0f, 5.0f, 5.0f, 1.0f };

void matColor(float kdr, float kdg, float kdb,  float shiny, int frnt_Back, float ambFactor, float specFactor)
{

    const GLfloat mat_ambient[]    = { kdr*ambFactor, kdg*ambFactor, kdb*ambFactor, 1.0f };
    const GLfloat mat_diffuse[]    = { kdr, kdg, kdb, 1.0f };
    const GLfloat mat_specular[]   = { 1.0f*specFactor, 1.0f*specFactor, 1.0f*specFactor, 1.0f };
    const GLfloat high_shininess[] = { shiny };
    if(frnt_Back==0)
    {
        glMaterialfv(GL_FRONT, GL_AMBIENT,   mat_ambient);
        glMaterialfv(GL_FRONT, GL_DIFFUSE,   mat_diffuse);
        glMaterialfv(GL_FRONT, GL_SPECULAR,  mat_specular);
        glMaterialfv(GL_FRONT, GL_SHININESS, high_shininess);
    }
    else if(frnt_Back==1)
    {
        glMaterialfv(GL_BACK, GL_AMBIENT,   mat_ambient);
        glMaterialfv(GL_BACK, GL_DIFFUSE,   mat_diffuse);
        glMaterialfv(GL_BACK, GL_SPECULAR,  mat_specular);
        glMaterialfv(GL_BACK, GL_SHININESS, high_shininess);
    }
    else if(frnt_Back==2)
    {
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,   mat_ambient);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mat_diffuse);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  mat_specular);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, high_shininess);
    }

}

/* Program entry point */
void myInit()
{
    // glClearColor(.1,.1,.1,1);

    glEnable(GL_LIGHT0);
    glEnable(GL_NORMALIZE);
    //glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_LIGHTING);

    glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
}

int main(int argc, char *argv[])
{
    glutInit(&argc, argv);
    glutInitWindowSize(wWidth, wHeight);
    glutInitWindowPosition(0,0);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

    glutCreateWindow("Practice");

    texture();

    glutDisplayFunc(display);
    glutKeyboardFunc(keyBoard);
    glutIdleFunc(idle);
    glutReshapeFunc(resize);
    //glClearColor(0, 0, 0, 1);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_NORMALIZE);
    glEnable(GL_BLEND);

    myInit();


    // glutKeyboardFunc(key);
    glutMouseFunc(processMouse);
    glutIdleFunc(idle);


    printf("Use 'w' to look up, 's' to look down, 'd' to look right, and 'a' to look left.\n");
    printf("Use 'i' to move camera up, 'k' to move camera down, 'l' to move camera right, and 'j' to move camera left with the look at point fixed.\n");
    printf("Use '+' to zoom in and '-' to zoom out.\n\n\n");

    glutMainLoop();

    return 0;

}
//G:\4-2\dld lab\demo\project\New folder (2)


