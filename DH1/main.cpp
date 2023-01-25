#include <windows.h>
#include "gl.h"
#include "glu.h"

///Copyright (C) <2017>  <Eliseo Rivera> curso.tareas@gmail.com
// 


//#include "ClaseOpenGL.h"

 // start claseopengl.h
//#include "Robot.h"

 // begin Robot.h
 //#include "modelo3D.h"

 // begin modelo3D.h

//#include "Triangle3D.h"

// begin Triangle3D.h

//#include"vector3d.h"


// begin vector3d.h

#include <math.h>

class vector3d
{
    public:
    double x;
    double y;
    double z;
    vector3d(double _x, double _y, double _z);
    vector3d();
    ~vector3d();
    vector3d (const vector3d &p);
    vector3d operator+(const vector3d &p)const;
    vector3d operator*(const vector3d &p)const;
    vector3d operator-(const vector3d &p)const;
    friend vector3d operator*(const double &s, const vector3d &p);
    int operator==(vector3d &p);
    int operator!=(vector3d &p);
    double& operator[](int i);
    double dotProduct(const vector3d &p);
    double magnitude() const;
    vector3d&   normalize();
    vector3d&  operator =( const vector3d &B);

    void mostrar();
};

 // end vector3d.h
 // begin vector3d.cpp
 
 
#include <iostream>
using namespace std;
vector3d::vector3d(double _x, double _y, double _z): x(_x), y(_y), z(_z) {};
vector3d::vector3d(){ x=0;  y=0; z=0;};
vector3d::vector3d (const vector3d &p){x=p.x;y=p.y;z=p.z;};
vector3d vector3d::operator+(const vector3d &p)const
{ return vector3d(x + p.x, y + p.y, z + p.z); }
vector3d vector3d::operator*(const vector3d &p)const
{ return vector3d(y*p.z-z*p.y, z*p.x-x*p.z ,x*p.y-y*p.x); }
vector3d vector3d::operator-(const vector3d &p) const
{ return vector3d(x - p.x, y - p.y, z - p.z); }
vector3d operator*(const double &s,const vector3d &p){
    return vector3d(s * p.x, s * p.y, s * p.z);
}
int vector3d::operator==(vector3d &p)
{ return ((x == p.x) && (y == p.y) && (z == p.z));}
int vector3d::operator!=(vector3d &p)
{ return !(*this == p); }
double& vector3d::operator[](int i)
{ return ((i == 0) ? x : ((i == 1) ? y : z)); }
double vector3d:: dotProduct(const vector3d &p)
{ return (x*p.x + y*p.y + z*p.z); }

double vector3d::magnitude() const
{ return sqrt(x*x + y*y + z*z);
};
vector3d&   vector3d::normalize(){
double a=magnitude();
if (a!=0){
     x=x/a;  y=y/a; z=z/a;
};
return *this;
 }
  vector3d&  vector3d::operator =( const vector3d &B){
    x=B.x;  y=B.y; z=B.z;
 return *this;
 }
vector3d::~vector3d(){}
  void vector3d::mostrar(){
          cout<<" x = "<<x<<endl;
          cout<<" y = "<<y<<endl;
          cout<<" z = "<<z<<endl;


      }

 // end vector3d.cpp



class Triangle3D
{
public:
Triangle3D();
~Triangle3D();
vector3d vertices[3];
vector3d N;
vector3d normal();
void definirRz(float dtheta);
void dibujar();
void rotar(float dtheta);
};

 // end Triangle3D.h
 // begin Triangle3D.cpp
 
 
Triangle3D::Triangle3D()
{
    //ctor
}

Triangle3D::~Triangle3D()
{
    //dtor
}
vector3d Triangle3D::normal(){
vector3d d1,d2,n;
d1=vertices[1]-vertices[0];
d2=vertices[2]-vertices[0];
n=d1*d2;  ///devuelve el producto vectorial
n.normalize();
N=n;  ///se guarda el valor de la normal en el campo N
return n;
}

 // end Triangle3D.cpp

 // begin Matrix.h
 

#include <iostream>
class Matrix
{
public:
   double **aij;
   int n;
   int m;
   Matrix(vector3d v);
   Matrix (vector3d v, float w);
   double cofactor (int i, int j);
   Matrix inversa();
   int size () const;
   Matrix Mij(int a,int b);
   vector3d operator *(const vector3d &P);
   Matrix(int nn,int mm);
   Matrix();
   void zero(int nn, int mm);
   void identity(int nn);
   void resetIdentity();
   void definir();
   void mostrar();
   double& ij(int i,int j);
   Matrix(const Matrix &B);
   Matrix operator +( const Matrix &A);
   Matrix& operator =( const Matrix &A);
   Matrix operator -( Matrix A);
   Matrix operator *( Matrix B);
   double& entry(int i, int j);

   friend  Matrix operator*(double tt, Matrix &A);
  ~Matrix();
   double determinante();
};
 // end Matrix.h
 
 // begin Matrix.cpp
 
 
using namespace std;
Matrix::Matrix(vector3d v){
int nn=3; int mm=1;

if (nn>0&&mm>0){
  aij=new double*[nn];
  for(int i = 0; i < nn; i++)  aij[i] = new double[mm];
  n=nn;
  m=mm;

  for( int i = 0; i < n; i++)
    for(int j = 0; j < m; j++){
    //    coutincorrect <<"introduce ij ("<<i+1 <<","<<j+1<<")"<<endl;
           // cin>>ij(i,j);
aij[i][j]=0;
    };
}
ij(0,0)=v.x;
ij(1,0)=v.y;
ij(2,0)=v.z;
}
Matrix::Matrix(vector3d v,float w){
int nn=4; int mm=1;

if (nn>0&&mm>0){
  aij=new double*[nn];
  for(int i = 0; i < nn; i++)  aij[i] = new double[mm];
  n=nn;
  m=mm;

  for( int i = 0; i < n; i++)
    for(int j = 0; j < m; j++){
    //    coutincorrect <<"introduce ij ("<<i+1 <<","<<j+1<<")"<<endl;
           // cin>>ij(i,j);
aij[i][j]=0;
    };
}

ij(0,0)=v.x;
ij(1,0)=v.y;
ij(2,0)=v.z;
ij(3,0)=1;
}
double Matrix::cofactor (int i, int j){
double c;
c=pow (-1,i+j)*Mij(i,j).determinante();  //cofactor
return c;
};
Matrix Matrix::inversa(){
double d=(*this).determinante();

if (d==0) { Matrix B(1,1); cout<<" no existe la inversa"<<endl; return B;}
 if(n>2) {
    Matrix B(n,n);
     for( int i = 0; i < n; i++)
    for(int j = 0; j < n; j++){
         B.aij[j][i]=(1.0/d)*cofactor(i,j);


    }
    return B;

}
else if(n==2) {
    Matrix B(n,n);
     B.aij[0][0]=(1.0/d)*aij[1][1];

      B.aij[0][1]=-(1.0/d)*aij[0][1];
       B.aij[1][0]=-(1.0/d)*aij[1][0];
      B.aij[1][1]=(1.0/d)*aij[0][0];
    return B;

}
else if(n==1) {
    Matrix B(1,1);
     B.aij[0][0]=1/aij[0][0];


    return B;

}

};

    int Matrix::size () const{return n*m;};
    Matrix Matrix::Mij(int a,int b){ //menor de una matriz
        if (n==m&&size()>0){ // eliminar fila a y columna b

    int s=n-1;
    Matrix B(s,s); //matriz de tamaño reducido
    for (int i=0; i<a;i++)
    for (int j=0; j<b;j++)

    {
      B.aij[i][j]=aij[i][j]; //primera parte
    }

      for (int i=0; i<a;i++)
    for (int j=b+1; j<m;j++)

    {
      B.aij[i][j-1]=aij[i][j];  //se recorren las columnas
    }

        for (int i=a+1; i<n;i++)
    for (int j=0; j<b;j++)

    {
      B.aij[i-1][j]=aij[i][j];  //se recorren las columnas
    }


        for (int i=a+1; i<n;i++)
    for (int j=b+1; j<m;j++)

    {
      B.aij[i-1][j-1]=aij[i][j];   //se recorren las columnas
    }
    return B;
}
};
vector3d Matrix::operator *(const vector3d &P){//operacion multiplicacion por punto
Matrix A(3,1);
A.aij[0][0]=P.x;
A.aij[1][0]=P.y;
A.aij[2][0]=P.z;

Matrix T((*this)*A);

vector3d Pr;
Pr.x=T.aij[0][0];
Pr.y=T.aij[1][0];
Pr.z=T.aij[2][0];
return Pr;
}
Matrix::Matrix(int nn,int mm){ //constructor a partir de n y m
if (nn>0&&mm>0){
  aij=new double*[nn];
  for(int i = 0; i < nn; i++)  aij[i] = new double[mm];
  n=nn;
  m=mm;

  for( int i = 0; i < n; i++)
    for(int j = 0; j < m; j++){
    //    coutincorrect <<"introduce ij ("<<i+1 <<","<<j+1<<")"<<endl;
           // cin>>ij(i,j);
aij[i][j]=0;
    };
} ;   //numero de filas
//coutincorrect<<"\n"<<endl;
};
Matrix::Matrix(){



n=m=0;
  //numero de filas

};
void Matrix::zero(int nn, int mm){
     n=nn;
        m=mm;
     //constructor de mat


if (n>0&&m>0){



  aij=new double*[n];
  for(int i = 0; i < n; i++)  aij[i] = new double[m];


  for( int i = 0; i < n; i++)
    for(int j = 0; j < m; j++){

        ij(i,j)=0.0;;

    };
} ;   //numero de filas
//coutincorrect<<"\n"<<endl;
};
 void Matrix::resetIdentity(){

  for( int i = 0; i < n; i++)
    for(int j = 0; j < m; j++){
if(  i==j)ij(i,j)=1.0;
else ij(i,j)=0.0;

    };
 }

void Matrix::identity(int nn){
     n=nn;
        m=nn;
     //constructor de mat


if (n>0&&m>0){



  aij=new double*[n];
  for(int i = 0; i < n; i++)  aij[i] = new double[m];


  for( int i = 0; i < n; i++)
    for(int j = 0; j < m; j++){
if(  i==j)ij(i,j)=1.0;
else ij(i,j)=0.0;

    };
} ;   //numero de filas
//coutincorrect<<"\n"<<endl;
};

void Matrix::definir(){
    cout <<"Creando Matriz "<<endl;
     cout <<"introduce el numero de filas "<<endl;
     cin>> n;
         cout <<"introduce  el numero de columnas "<<endl;
              cin>> m;
     //constructor de mat


if (n>0&&m>0){



  aij=new double*[n];
  for(int i = 0; i < n; i++)  aij[i] = new double[m];


  for( int i = 0; i < n; i++)
    for(int j = 0; j < m; j++){
        cout<<"introduce ij ("<<i+1 <<","<<j+1<<")"<<endl;
        cin>>ij(i,j);

    };
} ;   //numero de filas
//coutincorrect<<"\n"<<endl;
};
void Matrix::mostrar(){
  if (this->n>0&&this->m>0){


  cout<<"\n"<<endl;
for( int i = 0; i < n; i++){
    for( int j = 0; j < m; j++){

      cout<<ij(i,j)<<"  "<<flush;
    }
    cout<<"\n"<<endl;
    }
}
}
double& Matrix::ij(int i,int j){ //devuelve la direccion
if (i>=0&&j>=0&&i<n&&j<m) return aij[i][j]; //c++ numera las ijs de una matrix a partir de cero, nosotros a partir de 1

};
Matrix::Matrix(const Matrix &B){
if(B.size()!=0){
        n=B.n;
        m=B.m;


        aij=new double*[n];
  for(int i = 0; i < n; i++)  aij[i] = new double[m];

    for( int i = 0; i < n; i++)
    for( int j = 0; j < m; j++){

       aij[i][j]=B.aij[i][j];

    }

}
};
Matrix Matrix::operator +( const Matrix &A){
if (n>0&&m>0){
 if (n==A.n&&m==A.m){

Matrix B(n,m);

for( int i = 0; i < n; i++)
    for( int j = 0; j < m; j++){

        B.aij[i][j]= aij[i][j]+A.aij[i][j];

    }

      return B;
    }
    else{Matrix B(1,1); cout<<" Tamaños diferentes"<<endl;  return B;};

    }


};
Matrix& Matrix::operator =( const Matrix &A){
if (A.n>0&&A.m>0&&n==A.n&&m==A.m){

Matrix B(n,m);
for( int i = 0; i < n; i++)
    for(int j = 0; j < m; j++){

       aij[i][j]=A.aij[i][j];
    }
      return *this;

     }   else{Matrix B(1,1); cout<<" Tamaños diferentes"<<endl;  return *this;}

};
Matrix Matrix::operator -( Matrix A){
if (this->n>0&&this->m>0){
 if (this->n==A.n&&this->m==A.m){
Matrix B(n,m);
for( int i = 0; i < this->n; i++)
    for(int j = 0; j < this->m; j++){

        B.ij(i,j)= ij(i,j)-A.ij(i,j);
    }
      return B;
    }    else{Matrix B(1,1); cout<<" Tamaños diferentes"<<endl;  return B;}
}

};
Matrix Matrix::operator *( Matrix B){ //multiplicaciion de matrices
if (size()>0&&B.size()>0){
 if (m==B.n){
Matrix A(n,B.m);

for( int i = 0; i < n; i++)
    for(int j = 0; j < B.m; j++){
            float s=0;
    for(int k=0;k<m;k++){s=s+ij(i,k)*B.ij(k,j); }
        A.ij(i,j)= s;
    }
      return A;
    }    else{Matrix B(1,1); cout<<" Tamaños diferentes"<<endl;  return B;} }


};
Matrix operator*(double tt, Matrix &A){
if (A.n>0&&A.m>0){

Matrix B(A.n,A.m);
for(int i = 0; i < A.n; i++)
    for( int j = 0; j < A.m; j++){

        B.ij(i,j)= tt*A.ij(i,j);
    }


    return B;
     }


};
Matrix::~Matrix(){


			if (n>0){
				for (int i = 0; i<n; i++){
					delete[] aij[i];  //se elimnan las columnas
				}
				delete[] aij; //eliminar las filas
				aij = NULL;
			}
	}

double Matrix::determinante(){  //se define en forma recursiva
double d=0;
if (n!=m) {cout<<"la matriz no es cuadrada"<<endl; return 0; }
 if (n>2){
for (int j=0; j<m;j++){

    d=d+aij[0][j]*pow (-1,j)*Mij(0,j).determinante();

}
return d;
} else if (n==2) {
    d=aij[0][0]*aij[1][1]-aij[0][1]*aij[1][0];
    return d;

} else if(n==1){

return aij[0][0];
};

}

double& Matrix::entry(int i, int j){
return aij[i][j];
}

 
 // end Matrix.cpp


#include <fstream>
#include <iostream>

//#include "gl.h"
//#include "glu.h"
#define PI 3.14159265
#include <string>
//Copyright (C) <2017>  <Eliseo Rivera Silva> curso.tareas@gmail.com

using namespace std;
class modelo3D
{
    public:
int ntriangles;
Triangle3D* triangulos;
Matrix Rx,Ry,Rz;  //local
modelo3D();

void leer(string nombre);
virtual ~modelo3D();
modelo3D(int ntriangulos);
void dibujar();
void rotarX();  //local x axis rotation
void rotarY();  //local y axis rotation
void rotarZ();  //local z axis rotation
void definirRz(float theta);
void definirRy(float theta);
void definirRx(float theta);

void definir_x_LocalAxisRotation(float theta);
void x_LocalAxisRotation();
void definir_y_LocalAxisRotation(float theta);
void y_LocalAxisRotation();
void definir_z_LocalAxisRotation(float theta);
void z_LocalAxisRotation();
Matrix LARx, LARy, LARz;


void trasladar(vector3d A);
ifstream archivo;

vector3d ux,uy,uz,O; //local axis nad Origin
void BodyFramedefinirRx(float theta); //rotar marco de referencia alrededor del eje X, rota todo
void BodyFramedefinirRy(float theta);
void BodyFramedefinirRz(float theta);
void BodyFrametrasladar(vector3d A);
void BodyFramerotarX();
void BodyFramerotarY();
void BodyFramerotarZ();
Matrix BFRx,BFRy,BFRz,R;
vector3d LocalMassCenter()const;
vector3d GlobalCenterMass()const;


};

 // end modelo3D.h
 
 // begin modelo3D.cpp
#include <stdint.h>


modelo3D::modelo3D()
{
    ntriangles=0;
    ux={1,0,0};uy={0,1,0};uz={0,0,1};  //local and global axis are coincident

    BFRx.identity(3);
    BFRy.identity(3);
    BFRz.identity(3);
    LARx.identity(3);
    LARy.identity(3);
    LARz.identity(3);
    R.identity(3);  //accumulated rotations
    Rz.identity(3);
    Rx.identity(3);
    Ry.identity(3);
    ///***************
  //IRx.identity(3);
    //IRy.identity(3);
    //IRz.identity(3);
    ///***************
    O={0,0,0};
}

modelo3D::~modelo3D()
{
 delete triangulos;

}
modelo3D::modelo3D(int ntriangulos){
triangulos=new Triangle3D[ntriangulos];
ntriangles=ntriangulos;
}
void modelo3D::leer(string nombre){

char head[80] = "";
///************
	archivo.open(nombre.c_str(), std::ios::in | std::ios::binary);
	if(archivo){
	archivo.read(head, 80);

	int32_t size ;
	archivo.read(reinterpret_cast<char *> (&size), sizeof(int32_t));
	ntriangles=size;
	triangulos=new Triangle3D[size];
	Triangle3D triangle;
	vector3d P0, P1, P2;
	vector3d normal;
	float p0[3], p1[3], p2[3], n[3];
	char attribute[2] = "0";
for (int i=0;i<size;i++){

		archivo.read(reinterpret_cast<char *> (&n[0]), 4);
		archivo.read(reinterpret_cast<char *> (&n[1]), 4);
		archivo.read(reinterpret_cast<char *> (&n[2]), 4);
    //  cout<<n[0]<<" , "<<n[1]<<" , "<<n[2]<<endl;
       triangulos[i].N={n[0],n[1],n[2]};
		archivo.read(reinterpret_cast<char *> (&p0[0]), 4);
		archivo.read(reinterpret_cast<char *> (&p0[1]), 4);
		archivo.read(reinterpret_cast<char *> (&p0[2]), 4);
	//	   cout<<p0[0]<<" , "<<p0[1]<<" , "<<p0[2]<<endl;
          triangulos[i].vertices[0]={p0[0],p0[1],p0[2]};


		archivo.read(reinterpret_cast<char *> (&p1[0]), 4);
		archivo.read(reinterpret_cast<char *> (&p1[1]), 4);
		archivo.read(reinterpret_cast<char *> (&p1[2]), 4);
	//	   cout<<p1[0]<<" , "<<p1[1]<<" , "<<p1[2]<<endl;
          triangulos[i].vertices[1]={p1[0],p1[1],p1[2]};

		archivo.read(reinterpret_cast<char *> (&p2[0]), 4);
		archivo.read(reinterpret_cast<char *> (&p2[1]), 4);
		archivo.read(reinterpret_cast<char *> (&p2[2]), 4);
	//	   cout<<p2[0]<<" , "<<p2[1]<<" , "<<p2[2]<<endl;

        triangulos[i].vertices[2]={p2[0],p2[1],p2[2]};
	    archivo.read(attribute, 2);
	}


	archivo.close();
}
else{
        ntriangles=0;
        cout<<"el archivo no se encuentra"<<endl;
}
}
void modelo3D::definirRx(float dtheta){


Rx.aij[0][0]=1;
Rx.aij[0][1]=0;
Rx.aij[0][2]=0;

Rx.aij[1][0]=0;
Rx.aij[1][1]=cos(dtheta);
Rx.aij[1][2]=-sin(dtheta);

Rx.aij[2][0]=0;
Rx.aij[2][1]=sin(dtheta);
Rx.aij[2][2]=cos(dtheta);

}
void modelo3D::definirRy(float dtheta){


Ry.aij[0][0]=cos(dtheta);
Ry.aij[0][1]=0;
Ry.aij[0][2]=sin(dtheta);

Ry.aij[1][0]=0;
Ry.aij[1][1]=1;
Ry.aij[1][2]=0;

Ry.aij[2][0]=-sin(dtheta);
Ry.aij[2][1]=0;
Ry.aij[2][2]=cos(dtheta);



}
void modelo3D::definirRz(float dtheta){

Rz.aij[0][0]=cos(dtheta);
Rz.aij[0][1]=-sin(dtheta);
Rz.aij[0][2]=0;

Rz.aij[1][0]=sin(dtheta);
Rz.aij[1][1]=cos(dtheta);
Rz.aij[1][2]=0;

Rz.aij[2][0]=0;
Rz.aij[2][1]=0;
Rz.aij[2][2]=1;

}

void modelo3D::dibujar(){
    glBegin(GL_TRIANGLES);
    glFrontFace(GL_FRONT_AND_BACK);

    for (int i=0;i<ntriangles;i++){

vector3d v1=triangulos[i].vertices[0];   //locales son locales
vector3d v2=triangulos[i].vertices[1];
vector3d v3=triangulos[i].vertices[2];
vector3d V1,V2, V3;
V1=v1.x*ux+v1.y*uy+v1.z*uz;  // se dibuja en base local en origen de modelo local
V2=v2.x*ux+v2.y*uy+v2.z*uz;
V3=v3.x*ux+v3.y*uy+v3.z*uz;

V1=O+V1; V2=O+V2; V3=O+V3;

vector3d d1,d2,n;
d1=V2-V1;
d2=V3-V1;
n=d1*d2;  ///devuelve el producto vectorial
n.normalize();

glNormal3f(n.x,n.y,n.z);
        glVertex3f(V1.x,V1.y,V1.z);
        glVertex3f(V2.x,V2.y,V2.z);
        glVertex3f(V3.x,V3.y,V3.z);
    }

glEnd();
};

void modelo3D::rotarZ(){  //global rotation of only model
    for(int i=0;i<ntriangles;i++){
        vector3d v1,v2,v3;
        v1=triangulos[i].vertices[0];
        v2=triangulos[i].vertices[1];
        v3=triangulos[i].vertices[2];
v1=Rz*v1;
v2=Rz*v2;
v3=Rz*v3;

        triangulos[i].vertices[0]=v1;
        triangulos[i].vertices[1]=v2;
        triangulos[i].vertices[2]=v3;
    }
    R=Rz*R;
}
void modelo3D::rotarY(){ //global rotation  of only model
    for(int i=0;i<ntriangles;i++){
        vector3d v1,v2,v3;
        v1=triangulos[i].vertices[0];
        v2=triangulos[i].vertices[1];
        v3=triangulos[i].vertices[2];
v1=Ry*v1;
v2=Ry*v2;
v3=Ry*v3;

        triangulos[i].vertices[0]=v1;
        triangulos[i].vertices[1]=v2;
        triangulos[i].vertices[2]=v3;
    }
    R=Ry*R;
}
void modelo3D::rotarX(){//global rotation  of only model
    for(int i=0;i<ntriangles;i++){
        vector3d v1,v2,v3;
        v1=triangulos[i].vertices[0];
        v2=triangulos[i].vertices[1];
        v3=triangulos[i].vertices[2];
v1=Rx*v1;
v2=Rx*v2;
v3=Rx*v3;

        triangulos[i].vertices[0]=v1;
        triangulos[i].vertices[1]=v2;
        triangulos[i].vertices[2]=v3;
    }
    R=Rx*R;
}

void modelo3D::trasladar(vector3d A){  //local traslation
 for(int i=0;i<ntriangles;i++){
        vector3d v1,v2,v3;
        v1=triangulos[i].vertices[0];
        v2=triangulos[i].vertices[1];
        v3=triangulos[i].vertices[2];


        triangulos[i].vertices[0]=v1+A;
        triangulos[i].vertices[1]=v2+A;
        triangulos[i].vertices[2]=v3+A;
    }


}

void modelo3D::BodyFramedefinirRx(float dtheta){


BFRx.aij[0][0]=1;
BFRx.aij[0][1]=0;
BFRx.aij[0][2]=0;

BFRx.aij[1][0]=0;
BFRx.aij[1][1]=cos(dtheta);
BFRx.aij[1][2]=-sin(dtheta);

BFRx.aij[2][0]=0;
BFRx.aij[2][1]=sin(dtheta);
BFRx.aij[2][2]=cos(dtheta);


}
void modelo3D::BodyFramedefinirRy(float dtheta){

BFRy.aij[0][0]=cos(dtheta);
BFRy.aij[0][1]=0;
BFRy.aij[0][2]=sin(dtheta);

BFRy.aij[1][0]=0;
BFRy.aij[1][1]=1;
BFRy.aij[1][2]=0;

BFRy.aij[2][0]=-sin(dtheta);
BFRy.aij[2][1]=0;
BFRy.aij[2][2]=cos(dtheta);


}
void modelo3D::BodyFramedefinirRz(float dtheta){

BFRz.aij[0][0]=cos(dtheta);
BFRz.aij[0][1]=-sin(dtheta);
BFRz.aij[0][2]=0;

BFRz.aij[1][0]=sin(dtheta);
BFRz.aij[1][1]=cos(dtheta);
BFRz.aij[1][2]=0;

BFRz.aij[2][0]=0;
BFRz.aij[2][1]=0;
BFRz.aij[2][2]=1;

}

void modelo3D::definir_x_LocalAxisRotation(float dtheta){
LARx.aij[0][0]=1;
LARx.aij[0][1]=0;
LARx.aij[0][2]=0;

LARx.aij[1][0]=0;
LARx.aij[1][1]=cos(dtheta);
LARx.aij[1][2]=sin(dtheta);

LARx.aij[2][0]=0;
LARx.aij[2][1]=-sin(dtheta);
LARx.aij[2][2]=cos(dtheta);
}
void modelo3D::definir_y_LocalAxisRotation(float dtheta){


LARy.aij[0][0]=cos(dtheta);
LARy.aij[0][1]=0;
LARy.aij[0][2]=sin(dtheta);

LARy.aij[1][0]=0;
LARy.aij[1][1]=1;
LARy.aij[1][2]=0;

LARy.aij[2][0]=-sin(dtheta);
LARy.aij[2][1]=0;
LARy.aij[2][2]=cos(dtheta);



}
void modelo3D::definir_z_LocalAxisRotation(float dtheta){

LARz.aij[0][0]=cos(dtheta);
LARz.aij[0][1]=-sin(dtheta);
LARz.aij[0][2]=0;

LARz.aij[1][0]=sin(dtheta);
LARz.aij[1][1]=cos(dtheta);
LARz.aij[1][2]=0;

LARz.aij[2][0]=0;
LARz.aij[2][1]=0;
LARz.aij[2][2]=1;

}


void  modelo3D::x_LocalAxisRotation(){
R=LARx*R;
vector3d Lx(1,0,0),Ly(0,1,0),Lz(0,0,1);///definir vectores
Lx=LARx*Lx;
Ly=LARx*Ly;
Lz=LARx*Lz;
ux=Lx.x*ux+Lx.y*uy+Lx.z*uz;
uy=Ly.x*ux+Ly.y*uy+Ly.z*uz;
uz=Lz.x*ux+Lz.y*uy+Lz.z*uz;
}
void  modelo3D::y_LocalAxisRotation(){
R=LARy*R;
vector3d Lx(1,0,0),Ly(0,1,0),Lz(0,0,1);///definir vectores
Lx=LARy*Lx;
Ly=LARy*Ly;
Lz=LARy*Lz;
ux=Lx.x*ux+Lx.y*uy+Lx.z*uz;
uy=Ly.x*ux+Ly.y*uy+Ly.z*uz;
uz=Lz.x*ux+Lz.y*uy+Lz.z*uz;
}
void  modelo3D::z_LocalAxisRotation(){
R=LARz*R;
vector3d Lx(1,0,0),Ly(0,1,0),Lz(0,0,1);///definir vectores
Lx=LARz*Lx;
Ly=LARz*Ly;
Lz=LARz*Lz;
ux=Lx.x*ux+Lx.y*uy+Lx.z*uz;
uy=Ly.x*ux+Ly.y*uy+Ly.z*uz;
uz=Lz.x*ux+Lz.y*uy+Lz.z*uz;
}

void modelo3D::BodyFramerotarX(){
R=BFRx*R;
ux=BFRx*ux;
uy=BFRx*uy;
uz=BFRx*uz;

};
void modelo3D::BodyFramerotarY(){
R=BFRy*R;
ux=BFRy*ux;
uy=BFRy*uy;
uz=BFRy*uz;
};
void modelo3D::BodyFramerotarZ(){
R=BFRz*R;
ux=BFRz*ux;
uy=BFRz*uy;
uz=BFRz*uz;
};
/*
void modelo3D::IRx(){
    Matrix::inversa(IRx());

}*/

void modelo3D::BodyFrametrasladar(vector3d A){
O=A;

};
vector3d modelo3D::LocalMassCenter()const
{
 vector3d center;
 for(int i=0;i<ntriangles;i++){
        vector3d v1,v2,v3;
        v1=triangulos[i].vertices[0];
        v2=triangulos[i].vertices[1];
        v3=triangulos[i].vertices[2];
center=center+(1.0/3.0)*(v1+v2+v3);

    }
    center=(1.0/ntriangles)*center;
    cout<<" Local center mass is :"<<endl;
    center.mostrar();
    return center;

}
vector3d modelo3D::GlobalCenterMass()const{
vector3d L=LocalMassCenter();



        vector3d C=O+L.x*ux+L.y*uy+L.z*uz;
    cout<<" Global center mass is: "<<endl;
        C.mostrar();
        return C;

}

 // end modelo3D.cpp

#include <vector>
#include <cstdlib>
///Copyright (C) <2017>  <Eliseo Rivera> curso.tareas@gmail.com
class Robot
{
    public:
        Robot();
        ~Robot();

        modelo3D *base;
           modelo3D *b1;
              modelo3D *b2;
                 modelo3D *b3;
                    modelo3D *b4;
                       modelo3D *b5;
                          modelo3D *b6;
                               modelo3D *gripe;

void inicializar();
void renderizar();
void configurarTH();

void AplicarTHx(float theta, vector3d d);
void AplicarTHy(float theta, vector3d d);
void AplicarTHz(float theta, vector3d d);
Matrix THx,THy,THz,TH;

std::vector<Matrix> THList;
std::vector<vector3d> Origenes;
std::vector<modelo3D*> modelos;


float theta1, theta2, theta3,theta4, theta5, theta6,theta7;
vector3d d1,d2,d3,d4,d5,d6,d7;
private :
void DefinirTHx(float theta, vector3d d);
void DefinirTHy(float theta, vector3d d);
void DefinirTHz(float theta, vector3d d);
void  Drawarrow3D( vector3d A,  vector3d B, vector3d color, double cota1,double R) ;


};

 // end robot.h
 // begin robot.cpp
 // I guess the Denavit–Hartenberg params are somewhere here
#define PI 3.14159265
///Copyright (C) <2017>  <Eliseo Rivera> curso.tareas@gmail.com
Robot::Robot()
{
   theta1=0;theta2=0;theta3=0; theta4=0; theta5=0; theta6=0,theta7=0; // I guess the default theta parameters. Where's alpha then?!
   THx.identity(4);
   THy.identity(4);
   THz.identity(4);
int d=4;
  TH.identity(4);


}

Robot::~Robot()
{

delete base;
delete b1;
delete b2;
delete b3;
delete b4;
delete b5;
delete b6;
delete gripe;
    //dtor
}

void Robot::inicializar(){
base=new modelo3D();
b1=new modelo3D();
b2=new modelo3D();
b3=new modelo3D();
b4=new modelo3D();
b5=new modelo3D();
b6=new modelo3D();
gripe=new modelo3D();

 // Read STL files
base->leer("base.STL");
b1->leer("b1.STL");
b2->leer("b2.STL");
b3->leer("b3.STL");
b4->leer("b4.STL");
b5->leer("b5.STL");
b6->leer("b6.STL");
gripe->leer("gripe.STL");


modelos.push_back(base);
modelos.push_back(b1);
modelos.push_back(b2);
modelos.push_back(b3);
modelos.push_back(b4);
modelos.push_back(b5);
modelos.push_back(b6);
modelos.push_back(gripe);

}

// Robot dimensions defined somewhere here!!!
void Robot::configurarTH(){
AplicarTHz(0,{0,0,0}); //base
THList.push_back(THz);
AplicarTHx(0,{0,0,0}); //base
THList.push_back(THx);

AplicarTHz(0,{0,0,2.5}); //b1 // fourth value here seems to be a or d value - something related to Z axis
//AplicarTHz(30,{0,0,2.5}); //b1 // changing first value simply rotates joint
THList.push_back(THz);
AplicarTHx(90,{0,0,0}); //b1 // First value here is apparently alpha! Woohooo!

//AplicarTHx(-90,{0,0,0}); //b1
THList.push_back(THx);

AplicarTHz(0,{0,0,6}); //b2 // changing second value just causes an initial mismatch of position but corrects itself on 1st movement
THList.push_back(THz);
AplicarTHx(-90,{0,0,0}); //b2

//AplicarTHx(-90,{0,0,0}); //b2
THList.push_back(THx);

AplicarTHz(0,{0,0,6}); //b3
THList.push_back(THz);
AplicarTHx(0,{12,0,0}); //b3 // Changed!!

//AplicarTHx(0,{12,0,0}); //b3
THList.push_back(THx);

AplicarTHz(0,{0,0,1}); //b4
THList.push_back(THz);
AplicarTHx(0,{12,0,0}); //b4
THList.push_back(THx);

AplicarTHz(0,{0,0,6}); //b5
THList.push_back(THz);
AplicarTHx(-90,{0,0,0}); //b5
THList.push_back(THx);

AplicarTHz(0,{0,0,6}); //b6
THList.push_back(THz);
AplicarTHx(-90,{0,0,0}); //b6
THList.push_back(THx);

AplicarTHz(0,{0,0,0.0}); //gripe
THList.push_back(THz);
AplicarTHx(0,{0,0,0}); //gripe
THList.push_back(THx);

}
void Robot::renderizar(){


TH.resetIdentity();

modelo3D *model;

for (int m=0;m<modelos.size();m++){

    model=modelos[m];
    TH=TH* THList[2*m+0]*THList[2*m+1];


vector3d ux,uy,uz,O;
ux={1,0,0};
uy={0,1,0};
uz={0,0,1};

Matrix ux4(ux,1),uy4(uy,1),uz4(uz,1),O4(O,1);


ux4=TH*ux4-TH*O4;
uy4=TH*uy4-TH*O4;
uz4=TH*uz4-TH*O4;
O4=TH*O4;


ux={ux4.aij[0][0],ux4.aij[1][0],ux4.aij[2][0]};
uy={uy4.aij[0][0],uy4.aij[1][0],uy4.aij[2][0]};
uz={uz4.aij[0][0],uz4.aij[1][0],uz4.aij[2][0]};
O={O4.aij[0][0],O4.aij[1][0],O4.aij[2][0]};

//if (m<2){
         Drawarrow3D(O,O+4*ux,{1,0.1,0.2},0.03,0.1);
         Drawarrow3D(O,O+4*uy,{.1,1,0.2},0.03,0.1);
         Drawarrow3D(O,O+4*uz,{0.1,0.2,1},0.03,0.1);
       //  }
         glColor4f(fabs(cos(m*PI/modelos.size())),fabs(sin(20*(m-5)*PI/modelos.size())),0.2,0.5);

glEnable(GL_BLEND);
 glBegin(GL_TRIANGLES);

  glFrontFace(GL_FRONT_AND_BACK);
    for (int i=0;i<model->ntriangles;i++){

vector3d v1=model->triangulos[i].vertices[0];   //posiciones locales
vector3d v2=model->triangulos[i].vertices[1];
vector3d v3=model->triangulos[i].vertices[2];
Matrix v14(v1,1),v24(v2,1),v34(v3,1);

v14=TH*v14;
v24=TH*v24;
v34=TH*v34;
v1={v14.entry(0,0),v14.entry(1,0),v14.entry(2,0)};
v2={v24.entry(0,0),v24.entry(1,0),v24.entry(2,0)};
v3={v34.entry(0,0),v34.entry(1,0),v34.entry(2,0)};



Matrix N(4,1),d14(4,1),d24(4,1);
d14=v24-v14;
d24=v34-v14;
vector3d d1,d2,n;
d1={d14.entry(0,0),d14.entry(1,0),d14.entry(2,0)};
d2={d24.entry(0,0),d24.entry(1,0),d24.entry(2,0)};
n=d1*d2;  ///devuelve el producto vectorial
n.normalize();



        glNormal3f(n.x,n.y,n.z);
        glVertex3f(v1.x,v1.y,v1.z);
        glVertex3f(v2.x,v2.y,v2.z);
        glVertex3f(v3.x,v3.y,v3.z);
    }
glEnd();
// }
 glDisable(GL_BLEND);


///DIBUJAR EJES


//}
}

}




void Robot::DefinirTHx(float dtheta, vector3d d){

THx.aij[0][0]=1;
THx.aij[0][1]=0;
THx.aij[0][2]=0;
THx.aij[0][3]=d.x;

THx.aij[1][0]=0;
THx.aij[1][1]=cos(dtheta);
THx.aij[1][2]=-sin(dtheta);
THx.aij[1][3]=d.y;

THx.aij[2][0]=0;
THx.aij[2][1]=sin(dtheta);
THx.aij[2][2]=cos(dtheta);
THx.aij[2][3]=d.z;

THx.aij[3][0]=0;
THx.aij[3][1]=0;
THx.aij[3][2]=0;
THx.aij[3][3]=1;

}
void Robot::DefinirTHy(float dtheta, vector3d d){


THy.aij[0][0]=cos(dtheta);
THy.aij[0][1]=0;
THy.aij[0][2]=sin(dtheta);

THy.aij[1][0]=0;
THy.aij[1][1]=1;
THy.aij[1][2]=0;

THy.aij[2][0]=-sin(dtheta);
THy.aij[2][1]=0;
THy.aij[2][2]=cos(dtheta);

THy.aij[3][0]=0;
THy.aij[3][1]=0;
THy.aij[3][2]=0;
THy.aij[3][3]=1;

THy.aij[0][3]=d.x;
THy.aij[1][3]=d.y;
THy.aij[2][3]=d.z;
}
void Robot::DefinirTHz(float dtheta, vector3d d){

THz.aij[0][0]=cos(dtheta);
THz.aij[0][1]=-sin(dtheta);
THz.aij[0][2]=0;
THz.aij[0][3]=d.x;

THz.aij[1][0]=sin(dtheta);
THz.aij[1][1]=cos(dtheta);
THz.aij[1][2]=0;
THz.aij[1][3]=d.y;

THz.aij[2][0]=0;
THz.aij[2][1]=0;
THz.aij[2][2]=1;
THz.aij[2][3]=d.z;

THz.aij[3][0]=0;
THz.aij[3][1]=0;
THz.aij[3][2]=0;
THz.aij[3][3]=1;

}

void Robot::AplicarTHx(float theta, vector3d d){
theta=theta*PI/180.0;

DefinirTHx(theta,d);

}
void Robot::AplicarTHy(float theta, vector3d d){
theta=theta*PI/180.0;
DefinirTHy(theta,d);

}
void Robot::AplicarTHz(float theta, vector3d d){
theta=theta*PI/180.0;
DefinirTHz(theta,d);

}

void Robot::Drawarrow3D( vector3d A,  vector3d B, vector3d color, double cota1,double R)
{

double color1,color2,color3,a,b,c,d,e;



color1=color.x;//abs(color1/255);
color2=color.y;//abs(color2/255);
color3=color.z;//abs(color3/255);

glColor3f( color1,color2, color3);

vector3d n=B-A,np,vertex[10],normallight;
n.normalize();
if(n.z!=0)np={1,1,(-1/n.z)*(n.x+n.y)};
else if(n.y!=0)np={1,(-1/n.y)*(n.x+n.z),1};
else np={(-1/n.x)*(n.y+n.z),1,1};

np.normalize();
vertex[0]=R*np;
vertex[2]=R*(n*np).normalize();
vertex[1]=R*((0.5)*(vertex[2]-vertex[0])+vertex[0]).normalize();
vertex[4]=R*(n*vertex[2]).normalize();
vertex[3]=R*((0.5)*(vertex[4]-vertex[2])+vertex[2]).normalize();
vertex[6]=R*(n*vertex[4]).normalize();
vertex[5]=R*((0.5)*(vertex[6]-vertex[4])+vertex[4]).normalize();
vertex[7]=R*((0.5)*(vertex[0]-vertex[6])+vertex[6]).normalize();
vertex[8]=vertex[0];
vertex[9]=vertex[1];
int nx=8;
double d_thetha,fraccion=0.1,radioflecha=R+.7*R;
d_thetha=2.0f*PI/nx;


  ///tubos
 glBegin( GL_TRIANGLE_STRIP );

         for(int i=0;i<9;i++)
               {

normallight=n*(vertex[i-1]-vertex[i+1]);
normallight.normalize();
glNormal3f(normallight.x, normallight.y, normallight.z);
                 glVertex3f(vertex[i].x+A.x,vertex[i].y+A.y,vertex[i].z+A.z);

glVertex3f(vertex[i].x+B.x-fraccion*(B.x-A.x),vertex[i].y+B.y-fraccion*(B.y-A.y),vertex[i].z+B.z-fraccion*(B.z-A.z));

    // top face

                }

glEnd();



//flecha
vertex[0]=radioflecha*np;
vertex[2]=radioflecha*(n*np).normalize();
vertex[1]=radioflecha*((0.5)*(vertex[2]-vertex[0])+vertex[0]).normalize();
vertex[4]=radioflecha*(n*vertex[2]).normalize();
vertex[3]=radioflecha*((0.5)*(vertex[4]-vertex[2])+vertex[2]).normalize();
vertex[6]=radioflecha*(n*vertex[4]).normalize();
vertex[5]=radioflecha*((0.5)*(vertex[6]-vertex[4])+vertex[4]).normalize();
vertex[7]=radioflecha*((0.5)*(vertex[0]-vertex[6])+vertex[6]).normalize();
vertex[8]=vertex[0];
vertex[9]=vertex[1];
vector3d Ap(B-fraccion*(B-A));



 glBegin( GL_TRIANGLE_STRIP );  //flecha

         for(int i=0;i<9;i++)
               {

normallight=n*(vertex[i-1]-vertex[i+1]);
normallight.normalize();
glNormal3f(normallight.x, normallight.y, normallight.z);
                 glVertex3f(vertex[i].x+Ap.x,vertex[i].y+Ap.y,vertex[i].z+Ap.z);


glNormal3f(n.x, n.y, n.z);
glVertex3f(Ap.x+fraccion*(B-A).x,Ap.y+fraccion*(B-A).y,Ap.z+fraccion*(B-A).z);

    // top face

                }

glEnd();


}

 // end robot.cpp


 // begin claseopengl.h

//#include "modelo3D.h" // already done above



#define PI 3.14159265

class ClaseOpenGL
{
    public:
        ClaseOpenGL();
        virtual ~ClaseOpenGL();
        void inicializar();
        void preparar(float dt);
        void renderizar();
        void Resize(int width, int height);
        float theta;
        void  Drawarrow3D( vector3d A,  vector3d B, vector3d color, double cota1,double R=0.01);

        float  cameraX,  cameraZ, cameraY ,   mouseX,mouseY, camerafactor,angle;
        float Rcamera, phiCamera, thetaCamera;
        float radians;

        vector3d rojo={1,0.0,0.0};
        vector3d verde={0.0,1,0.0};
        vector3d azul={0.0,0.0,1};
        void dibujarBodyFrame(const modelo3D &modelo);

        Robot SSRMS;
};
 // end claseopengl.h


 // start claseopengl.cpp

ClaseOpenGL::ClaseOpenGL(){
theta=0.0;

};

void ClaseOpenGL::inicializar(){
SSRMS.inicializar();///cargar modelos
SSRMS.configurarTH();///Posiciones de las piezas // Positions of pieces
    /////Camera default position
Rcamera=50;
phiCamera=PI/3.0;
thetaCamera=PI/4.0;
camerafactor=0.005;
////




    glClearColor(1,1,1,1);
    glEnable( GL_LINE_SMOOTH );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_DEPTH_TEST);					// hidden surface removal
    glShadeModel(GL_SMOOTH);					// use smooth shading
    glEnable(GL_LIGHTING);

    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHT2);
    glEnable(GL_COLOR_MATERIAL );

    }


ClaseOpenGL::~ClaseOpenGL()
{
}
void ClaseOpenGL::renderizar(){


        radians =  double(PI*(angle-90.0f)/180.0f);

        // calculate the camera's position

         cameraX=Rcamera*cos(phiCamera)*cos(thetaCamera);
         cameraY=Rcamera*cos(phiCamera)*sin(thetaCamera);
         cameraZ=Rcamera*sin(phiCamera);
     glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT );

    glLoadIdentity();
    gluLookAt(cameraX ,cameraY  ,cameraZ , 0, 0,0 , 0.0, 0,fabs(cos(thetaCamera))*sin(phiCamera+0.1));

// Draw coordinate system vectors
    Drawarrow3D({0,0,0},{100,0,0},rojo,0.03,0.03);
    Drawarrow3D({0,0,0},{0,100,0},verde,0.03,0.03);
    Drawarrow3D({0,0,0},{0,0,100},azul,0.03,0.03);
/// //////////////////////////// PARTES A AGREGAR
SSRMS.renderizar();

}
void ClaseOpenGL::preparar(float dt){
  //objeto->rotarX();
  //objeto->rotarZ();

}
void ClaseOpenGL::Resize(int width, int height){

  glViewport(0, 0, width, height);		// reset the viewport to new dimensions
        glMatrixMode(GL_PROJECTION);			// set projection matrix current matrix
        glLoadIdentity();						// reset projection matrix
        // calculate aspect ratio of window
        gluPerspective(52.0f,(GLdouble)width/(GLdouble)height,1.0f,1000.0f);
        glMatrixMode(GL_MODELVIEW);				// set modelview matrix
        glLoadIdentity();						// reset modelview matrix
}
void ClaseOpenGL::Drawarrow3D( vector3d A,  vector3d B, vector3d color, double cota1,double R)
{

double color1,color2,color3,a,b,c,d,e;



color1=color.x;//abs(color1/255);
color2=color.y;//abs(color2/255);
color3=color.z;//abs(color3/255);

glColor4f( color1,color2, color3,0.8);

vector3d n=B-A,np,vertex[10],normallight;
n.normalize();
if(n.z!=0)np={1,1,(-1/n.z)*(n.x+n.y)};
else if(n.y!=0)np={1,(-1/n.y)*(n.x+n.z),1};
else np={(-1/n.x)*(n.y+n.z),1,1};

np.normalize();
vertex[0]=R*np;
vertex[2]=R*(n*np).normalize();
vertex[1]=R*((0.5)*(vertex[2]-vertex[0])+vertex[0]).normalize();
vertex[4]=R*(n*vertex[2]).normalize();
vertex[3]=R*((0.5)*(vertex[4]-vertex[2])+vertex[2]).normalize();
vertex[6]=R*(n*vertex[4]).normalize();
vertex[5]=R*((0.5)*(vertex[6]-vertex[4])+vertex[4]).normalize();
vertex[7]=R*((0.5)*(vertex[0]-vertex[6])+vertex[6]).normalize();
vertex[8]=vertex[0];
vertex[9]=vertex[1];
int nx=8;
double d_thetha,fraccion=0.1,radioflecha=R+.7*R;
d_thetha=2.0f*PI/nx;


  ///tubos
 glBegin( GL_TRIANGLE_STRIP );

         for(int i=0;i<9;i++)
               {

normallight=n*(vertex[i-1]-vertex[i+1]);
normallight.normalize();
glNormal3f(normallight.x, normallight.y, normallight.z);
                 glVertex3f(vertex[i].x+A.x,vertex[i].y+A.y,vertex[i].z+A.z);

glVertex3f(vertex[i].x+B.x-fraccion*(B.x-A.x),vertex[i].y+B.y-fraccion*(B.y-A.y),vertex[i].z+B.z-fraccion*(B.z-A.z));

    // top face

                }

glEnd();



//flecha
vertex[0]=radioflecha*np;
vertex[2]=radioflecha*(n*np).normalize();
vertex[1]=radioflecha*((0.5)*(vertex[2]-vertex[0])+vertex[0]).normalize();
vertex[4]=radioflecha*(n*vertex[2]).normalize();
vertex[3]=radioflecha*((0.5)*(vertex[4]-vertex[2])+vertex[2]).normalize();
vertex[6]=radioflecha*(n*vertex[4]).normalize();
vertex[5]=radioflecha*((0.5)*(vertex[6]-vertex[4])+vertex[4]).normalize();
vertex[7]=radioflecha*((0.5)*(vertex[0]-vertex[6])+vertex[6]).normalize();
vertex[8]=vertex[0];
vertex[9]=vertex[1];
vector3d Ap(B-fraccion*(B-A));



 glBegin( GL_TRIANGLE_STRIP );  //flecha

         for(int i=0;i<9;i++)
               {

normallight=n*(vertex[i-1]-vertex[i+1]);
normallight.normalize();
glNormal3f(normallight.x, normallight.y, normallight.z);
                 glVertex3f(vertex[i].x+Ap.x,vertex[i].y+Ap.y,vertex[i].z+Ap.z);


glNormal3f(n.x, n.y, n.z);
glVertex3f(Ap.x+fraccion*(B-A).x,Ap.y+fraccion*(B-A).y,Ap.z+fraccion*(B-A).z);

    // top face

                }

glEnd();


}

void ClaseOpenGL::dibujarBodyFrame(const modelo3D &model){

        Drawarrow3D(model.O,model.O+80*model.ux,rojo,0.03,1);
        Drawarrow3D(model.O,model.O+80*model.uy,verde,0.03,1);
        Drawarrow3D(model.O,model.O+80*model.uz,azul,0.03,1);

}
 // end claseopengl.cpp

// begin main.cpp

LRESULT CALLBACK WindowProc(HWND, UINT, WPARAM, LPARAM);
void EnableOpenGL(HWND hwnd, HDC*, HGLRC*);
void DisableOpenGL(HWND, HDC, HGLRC);
vector3d d1,d2,d3,d4,d5,d6,d7;

float dtheta(3);



 ClaseOpenGL* Miclase=NULL;
int WINAPI WinMain(HINSTANCE hInstance,
                   HINSTANCE hPrevInstance,
                   LPSTR lpCmdLine,
                   int nCmdShow)
{

    WNDCLASSEX wcex;
    HWND hwnd;
    HDC hDC;
    HGLRC hRC;
    MSG msg;
    BOOL bQuit = FALSE;


    /* register window class */
    wcex.cbSize = sizeof(WNDCLASSEX);
    wcex.style = CS_OWNDC;
    wcex.lpfnWndProc = WindowProc;
    wcex.cbClsExtra = 0;
    wcex.cbWndExtra = 0;
    wcex.hInstance = hInstance;
    wcex.hIcon = LoadIcon(NULL, IDI_APPLICATION);
    wcex.hCursor = LoadCursor(NULL, IDC_ARROW);
    wcex.hbrBackground = (HBRUSH)GetStockObject(BLACK_BRUSH);
    wcex.lpszMenuName = NULL;
    wcex.lpszClassName = "GLSample";
    wcex.hIconSm = LoadIcon(NULL, IDI_APPLICATION);;

Miclase = new ClaseOpenGL;  //agregar código opengl

    if (!RegisterClassEx(&wcex))
        return 0;

    /* create main window */
    hwnd = CreateWindowEx(0,
                          "GLSample",
                          "Robótica. Configuración  Denavit–Hartenberg ",
                          WS_OVERLAPPEDWINDOW,
                          CW_USEDEFAULT,
                          CW_USEDEFAULT,
                          800,
                          600,
                          NULL,
                          NULL,
                          hInstance,
                          NULL);
EnableOpenGL(hwnd, &hDC, &hRC); //va primero esto
ShowWindow(hwnd, nCmdShow); // y despues mostrar ventana
UpdateWindow(hwnd);					// update the window

    /* enable OpenGL for the window */

Miclase->inicializar();
SwapBuffers(hDC);
    /* program main loop */
    while (!bQuit)
    {
        /* check for messages */
        if (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
        {
            /* handle or dispatch messages */
            if (msg.message == WM_QUIT)
            {
                bQuit = TRUE;
            }
            else
            {
                TranslateMessage(&msg);
                DispatchMessage(&msg);
            }
        }
        else
        {
            /* OpenGL animation code goes here */

            Miclase->renderizar();
            SwapBuffers(hDC);
            Miclase->preparar(0.001);

            Sleep (1);
        }
    }
delete Miclase;
Miclase=NULL;
    /* shutdown OpenGL */
    DisableOpenGL(hwnd, hDC, hRC);

    /* destroy the window explicitly */
    DestroyWindow(hwnd);

    return msg.wParam;
}

LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{int height, width;
    switch (uMsg)
    {
    case WM_CREATE:


    break;
        case WM_CLOSE:
            PostQuitMessage(0);
        break;
    case WM_SIZE:
		height = HIWORD(lParam);		// retrieve width and height
		width = LOWORD(lParam);

	Miclase->Resize(width, height);

    break;

        case WM_DESTROY:
            return 0;

        case WM_KEYDOWN:
        {
                Matrix P({0,0,0},1);
    ( Miclase->SSRMS.TH*P).mostrar();
            switch (wParam)
            {



          case 'R':

                Miclase->SSRMS.theta1=Miclase->SSRMS.theta1+dtheta;
                Miclase->SSRMS.AplicarTHz(Miclase->SSRMS.theta1, {0,0,2.5});     //b1
                Miclase->SSRMS.THList[2]=Miclase->SSRMS.THz;

                break;
            case 'F':

                Miclase->SSRMS.theta1=Miclase->SSRMS.theta1-dtheta;
                Miclase->SSRMS.AplicarTHz( Miclase->SSRMS.theta1, {0,0,2.5});     //b1
              Miclase->SSRMS.THList[2]=Miclase->SSRMS.THz;


                break;

             case 'T':

                Miclase->SSRMS.theta2=Miclase->SSRMS.theta2+dtheta;
                Miclase->SSRMS.AplicarTHz(   Miclase->SSRMS.theta2, {0,0,6});  //b2
                Miclase->SSRMS.THList[4]=Miclase->SSRMS.THz;

                break;
            case 'G':

                Miclase->SSRMS.theta2=Miclase->SSRMS.theta2-dtheta;
                Miclase->SSRMS.AplicarTHz(   Miclase->SSRMS.theta2, {0,0,6});     //b2
                Miclase->SSRMS.THList[4]=Miclase->SSRMS.THz;


                break;

         case 'Y':

                Miclase->SSRMS.theta3=Miclase->SSRMS.theta3+dtheta;
                Miclase->SSRMS.AplicarTHz(   Miclase->SSRMS.theta3, {0,0,6});     //b3
                Miclase->SSRMS.THList[6]=Miclase->SSRMS.THz;

                break;
            case 'H':

                Miclase->SSRMS.theta3=Miclase->SSRMS.theta3-dtheta;
                Miclase->SSRMS.AplicarTHz(   Miclase->SSRMS.theta3,{0,0,6});       //b3
                Miclase->SSRMS.THList[6]=Miclase->SSRMS.THz;

                break;



          case 'U':

                Miclase->SSRMS.theta4=Miclase->SSRMS.theta4+dtheta;
                Miclase->SSRMS.AplicarTHz(   Miclase->SSRMS.theta4, {0,0,1});   //b4
                Miclase->SSRMS.THList[8]=Miclase->SSRMS.THz;


                break;

              case 'J':

                Miclase->SSRMS.theta4=Miclase->SSRMS.theta4-dtheta;
                Miclase->SSRMS.AplicarTHz(   Miclase->SSRMS.theta4, {0,0,1});   //b4
                Miclase->SSRMS.THList[8]=Miclase->SSRMS.THz;


                break;


                case 'I':
                Miclase->SSRMS.theta5=Miclase->SSRMS.theta5+dtheta;
                Miclase->SSRMS.AplicarTHz(   Miclase->SSRMS.theta5, {0,0,6});    //b5
                Miclase->SSRMS.THList[10]=Miclase->SSRMS.THz;
                break;

              case 'K':

                Miclase->SSRMS.theta5=Miclase->SSRMS.theta5-dtheta;
                Miclase->SSRMS.AplicarTHz(   Miclase->SSRMS.theta5, {0,0,6});     //b5
                Miclase->SSRMS.THList[10]=Miclase->SSRMS.THz;


                break;

                case 'O':
                Miclase->SSRMS.theta6=Miclase->SSRMS.theta6+dtheta;
                Miclase->SSRMS.AplicarTHz(   Miclase->SSRMS.theta6, {0,0,6});     //b6
                Miclase->SSRMS.THList[12]=Miclase->SSRMS.THz;
                break;

                    case 'L':

                Miclase->SSRMS.theta6=Miclase->SSRMS.theta6-dtheta;
                Miclase->SSRMS.AplicarTHz(   Miclase->SSRMS.theta6,{0,0,6});      //b6
                Miclase->SSRMS.THList[12]=Miclase->SSRMS.THz;
                break;

                case 'P':
                Miclase->SSRMS.theta7=Miclase->SSRMS.theta7-dtheta;
                Miclase->SSRMS.AplicarTHz(   Miclase->SSRMS.theta7, {0,0,0});     //b7, efector final
                Miclase->SSRMS.THList[14]=Miclase->SSRMS.THz;
                break;

// Camera movements

                case 'A':
                Miclase->thetaCamera=Miclase->thetaCamera+0.05;

                break;
                  case 'D':
                Miclase->thetaCamera=Miclase->thetaCamera-0.05;

                break;

                case 'W':
                Miclase->phiCamera=Miclase->phiCamera+0.05;

                break;
                  case 'S':
                Miclase->phiCamera=Miclase->phiCamera-0.05;

                break;
                case VK_ESCAPE:
                    PostQuitMessage(0);
                break;
            }
        }
        break;

        default:
            return DefWindowProc(hwnd, uMsg, wParam, lParam);
    }

    return 0;
}

void EnableOpenGL(HWND hwnd, HDC* hDC, HGLRC* hRC)
{
    PIXELFORMATDESCRIPTOR pfd;

    int iFormat;

    /* get the device context (DC) */
    *hDC = GetDC(hwnd);

    /* set the pixel format for the DC */
    ZeroMemory(&pfd, sizeof(pfd));

    pfd.nSize = sizeof(pfd);
    pfd.nVersion = 1;
    pfd.dwFlags = PFD_DRAW_TO_WINDOW |
                  PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER;
    pfd.iPixelType = PFD_TYPE_RGBA;
    pfd.cColorBits = 24;
    pfd.cDepthBits = 16;
    pfd.iLayerType = PFD_MAIN_PLANE;

    iFormat = ChoosePixelFormat(*hDC, &pfd);

    SetPixelFormat(*hDC, iFormat, &pfd);

    /* create and enable the render context (RC) */
    *hRC = wglCreateContext(*hDC);

    wglMakeCurrent(*hDC, *hRC);
}

void DisableOpenGL (HWND hwnd, HDC hDC, HGLRC hRC)
{
    wglMakeCurrent(NULL, NULL);
    wglDeleteContext(hRC);
    ReleaseDC(hwnd, hDC);
}

