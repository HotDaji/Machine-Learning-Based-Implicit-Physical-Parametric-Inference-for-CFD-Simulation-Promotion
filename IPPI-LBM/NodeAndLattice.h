#pragma once
#ifndef NODEANDLATTICE_H
#define NODEANDLATTICE_H

#include "Point2D.h"
#include "Parameter.h"
 // for the use of 'cout'

using namespace std;
class Particlenode
{
public:
    Particlenode();
	~Particlenode();
public:
#ifdef D2Q9
    Point2D position;
    Point2D position_ref;
    Point2D vel_ref;//��ֵ���ٶ�
    Point2D vel;
    Point2D force;
    Point2D center;
 
    
#endif // D2Q9
    void UpdateParticleBoundary();//����IB���������յ��λ��
};

class Lattice 
{
public:
    Lattice();
    ~Lattice();


#ifdef D2Q9
public:
    
  
  
    REAL rho;
    Point2D u={0,0};
    Position2D pos;//
    REAL p;//press
    REAL tau;//
    REAL TauAtBeginning;
    Point2D force;//
    vector<REAL> feq;
    vector<REAL> fcol;
    vector<REAL> fnew;
    vector<double> filteredFnew;
    vector<REAL> lat_force;
    double Cs = 0.0;
    double avg_u = 0;
    double avg_v = 0;
    double avg_density = 0;
    double dis[9] = {0,0,0,0,0,0,0,0,0};//YMS插值时用到的
    double normalDis;//动态cs时用到的
    char style = 'f';//初始化-1
    double cp = 0;
    double CoherentCS = 0;
    
    Point2D absU;
    double absRHO;
    double abstau;
    double xflowtau;
    bool NNflag= false;
public:

    void InitialLatticeFeq();
    void Compute_feq();
    void Compute_Disforce();
   
    void ComputeBGKCollision();
  
    void CalMacro();

#ifdef RRLBM
    void ComputeRRCollision();
#endif
#ifdef MRT
    void ComputeMRTCollision();
    void MRTInitial();
#endif

#endif // D2Q9

};
#endif
