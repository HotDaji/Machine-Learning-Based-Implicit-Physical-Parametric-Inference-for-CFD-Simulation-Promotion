#include"NodeAndLattice.h"
#ifdef D2Q9
const int e[9][2] ={
    {0,0}, 
    {1,0}, 
    {0,1}, 
    {-1,0}, 
    {0,-1}, 
    {1,1}, 
    {-1,1},
    {-1,-1},
    {1,-1},             
    };
const REAL weight[9] = {4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36.};

#endif
#ifdef MRT
const double tm[9][9] = { 
 { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 },
 { -4.0, -1.0, -1.0, -1.0, -1.0, 2.0, 2.0 ,2.0 ,2.0 }, 
 { 4.0, -2.0, -2.0, -2.0, -2.0, 1.0, 1.0, 1.0, 1.0 }, 
 { 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0 }, 
 { 0.0, -2.0, 0.0, 2.0, 0.0, 1.0, -1.0, -1.0, 1.0 }, 
 { 0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0 }, 
 { 0.0, 0.0, -2.0, 0.0, 2.0, 1.0, 1.0, -1.0, -1.0 },
  { 0.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0 }, 
  { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0 } };
const double a1 = 1.0 / 36.0;
const double tminv[9][9] = { { 4 * a1, -4 * a1, 4 * a1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, { 4 * a1, -a1, -2 * a1, 6 * a1, -6 * a1, 0.0, 0.0, 9 * a1, 0.0 }, { 4 * a1, -a1, -2 * a1, 0.0, 0.0, 6 * a1, -6 * a1, -9 * a1, 0.0 }, { 4 * a1, -a1, -2 * a1, -6 * a1, 6 * a1, 0.0, 0.0, 9 * a1, 0.0 },{ 4 * a1, -a1, -2 * a1, 0.0, 0.0, -6 * a1, 6 * a1, -9 * a1, 0.0 }, {4 * a1, 2 * a1, a1, 6 * a1, 3 * a1, 6 * a1, 3 * a1, 0.0, 9 * a1 },{ 4 * a1, 2 * a1, a1, -6 * a1, -3 * a1, 6 * a1, 3 * a1, 0.0, -9 * a1 }, { 4 * a1, 2 * a1, a1, -6 * a1, -3 * a1, -6 * a1, -3 * a1, 0.0, 9 * a1 }, { 4 * a1, 2 * a1, a1, 6 * a1, 3 * a1, -6 * a1, -3 * a1, 0.0, -9 * a1 } };


// const double tm[9][9] = {	{  1,	 1,	 1,	 1,	 1,	 1,	 1,	 1,	 1 },
// 							        { -4,	-1,	-1,	-1,	-1,	 2,	 2,	 2,	 2 },
// 							        {  4,	-2,	-2,	-2,	-2,	 1,	 1,	 1,	 1 },
// 							        {  0,	 1,	 0,	-1,	 0,	 1,	-1,	-1,	 1 },
// 							        {  0,	-2,	 0,	 2,	 0,	 1,	-1,	-1,	 1 },
// 							        {  0,	 0,	 1,	 0,	-1,	 1,	 1,	-1,	-1 },
// 							        {  0,	 0,	-2,	 0,	 2,	 1,	 1,	-1,	-1 },
// 							        {  0,	 1,	-1,	 1,	-1,	 0,	 0,	 0,	 0 },
// 							        {  0,	 0,	 0,	 0,  0,  1,	-1,	 1,	-1 }  };							
// const REAL tminv[9][9] = {	{ 1.0/9.0,	-1.0/9.0 ,	 1.0/9.0 ,	 0.0	,	 0.0	 ,	 0.0	,	 0.0	 ,	 0.0	,	 0.0	 },
// 								        { 1.0/9.0,	-1.0/36.0,	-1.0/18.0,	 1.0/6.0,	-1.0/6.0 ,	 0.0	,	 0.0	 ,	 1.0/4.0,	 0.0	 },
// 								        { 1.0/9.0,	-1.0/36.0,	-1.0/18.0,	 0.0	,	 0.0	 ,	 1.0/6.0,	-1.0/6.0 ,	-1.0/4.0,	 0.0	 },
// 								        { 1.0/9.0,	-1.0/36.0,	-1.0/18.0,	-1.0/6.0,	 1.0/6.0 ,	 0.0	,	 0.0	 ,	 1.0/4.0,	 0.0	 },
// 								        { 1.0/9.0,	-1.0/36.0,	-1.0/18.0,	 0.0	,	 0.0	 ,	-1.0/6.0,	 1.0/6.0 ,	-1.0/4.0,	 0.0	 },
// 								        { 1.0/9.0,	 1.0/18.0,	 1.0/36.0,	 1.0/6.0,	 1.0/12.0,	 1.0/6.0,	 1.0/12.0,	 0.0	,	 1.0/4.0 },
// 								        { 1.0/9.0,	 1.0/18.0,	 1.0/36.0,	-1.0/6.0,	-1.0/12.0,	 1.0/6.0,	 1.0/12.0,	 0.0	,	-1.0/4.0 },
// 								        { 1.0/9.0,	 1.0/18.0,	 1.0/36.0,	-1.0/6.0,	-1.0/12.0,	-1.0/6.0,	-1.0/12.0,	 0.0	,	 1.0/4.0 },
// 								        { 1.0/9.0,	 1.0/18.0,	 1.0/36.0,	 1.0/6.0,	 1.0/12.0,	-1.0/6.0,	-1.0/12.0,	 0.0	,	-1.0/4.0 }  };
// #endif
#endif
Particlenode::Particlenode() {

}
Particlenode::~Particlenode(){}
Lattice::Lattice()
{

}
Lattice::~Lattice()
{

}
void Particlenode::UpdateParticleBoundary()
{
	position_ref.x += vel_ref.x;
	position_ref.y += vel_ref.y; 
}
void Lattice::InitialLatticeFeq()
{
	if (feq.empty())
	{
		for (int k = 0; k < 9; k++)
		{
			feq.push_back(0.0);
		}
	}

	if (fnew.empty())
	{
		for (int k = 0; k < 9; k++)
		{
			fnew.push_back(0.0);
		}
	}

	if (fcol.empty())
	{
		for (int k = 0; k < 9; k++)
		{
			fcol.push_back(0.0);
		}
	}


	if (filteredFnew.empty())
	{
		for (int k = 0; k < 9; k++)
		{
			filteredFnew.push_back(0.0);
		}
	}//初始化滤波器的fnew

	REAL w0 = rho * 4.0 / 9.0;
	REAL w1 = rho / 9.0;
	REAL w2 = rho / 36.0;

	REAL uu = 1 - 1.5 * (u.x * u.x + u.y * u.y);

	
	fcol[0] = fnew[0] = feq[0] = w0*( uu );
	fcol[1] = fnew[1] = feq[1] = w1*( uu + 3.0*(u.x) + 4.5*(u.x)*(u.x) );
	fcol[2] = fnew[2] = feq[2] = w1*( uu + 3.0*(u.y) + 4.5*(u.y)*(u.y) );
	fcol[3] = fnew[3] = feq[3] = w1*( uu - 3.0*(u.x) + 4.5*(u.x)*(u.x) );
	fcol[4] = fnew[4] = feq[4] = w1*( uu - 3.0*(u.y) + 4.5*(u.y)*(u.y) );
	fcol[5] = fnew[5] = feq[5] = w2*( uu + 3.0*(u.x+u.y) + 4.5*(u.x+u.y)*(u.x+u.y) );
	fcol[6] = fnew[6] = feq[6] = w2*( uu - 3.0*(u.x-u.y) + 4.5*(u.x-u.y)*(u.x-u.y) );
	fcol[7] = fnew[7] = feq[7] = w2*( uu - 3.0*(u.x+u.y) + 4.5*(u.x+u.y)*(u.x+u.y) );
	fcol[8] = fnew[8] = feq[8] = w2*( uu + 3.0*(u.x-u.y) + 4.5*(u.x-u.y)*(u.x-u.y) );

	
}

void Lattice::Compute_feq() 
{
	if (feq.empty())
	{
		for (int k = 0; k < 9; k++)
		{
			feq.push_back(0.0);
		}
	}

	REAL w0 = rho * (4.0 / 9.0);
	REAL w1 = rho / 9.0;
	REAL w2 = rho / 36.0;

	REAL uu = 1 - 1.5 * (u.x * u.x + u.y * u.y);

	feq[0] = w0*( uu );
	feq[1] = w1*( uu + 3.0*(u.x) + 4.5*(u.x)*(u.x) );
	feq[2] = w1*( uu + 3.0*(u.y) + 4.5*(u.y)*(u.y) );
	feq[3] = w1*( uu - 3.0*(u.x) + 4.5*(u.x)*(u.x) );
	feq[4] = w1*( uu - 3.0*(u.y) + 4.5*(u.y)*(u.y) );
	feq[5] = w2*( uu + 3.0*(u.x+u.y) + 4.5*(u.x+u.y)*(u.x+u.y) );
	feq[6] = w2*( uu - 3.0*(u.x-u.y) + 4.5*(u.x-u.y)*(u.x-u.y) );
	feq[7] = w2*( uu - 3.0*(u.x+u.y) + 4.5*(u.x+u.y)*(u.x+u.y) );
	feq[8] = w2*( uu + 3.0*(u.x-u.y) + 4.5*(u.x-u.y)*(u.x-u.y) );
}
void Lattice::ComputeBGKCollision() 
{

	Compute_feq();
	for (int k = 0; k < 9; k++)
	{
		REAL FCOL = (1.0 - 1.0 / tau) * fnew[k] + feq[k] / tau ; 
		fcol[k] = FCOL;

	}

}
void Lattice::Compute_Disforce() 
{

	if (lat_force.empty()) {
		for (int k = 0; k < 9; k++)
		{
			lat_force.push_back(0.0);
		}
	}
	
	// REAL w0 = 4.0 / 9.0;
	// REAL w1 = 1.0 / 9.0;
	// REAL w2 = 1.0 / 36.0;

	REAL omega = 1.0 - 0.5*( 1.0 / TauAtBeginning);
	// REAL omega = 1.0 - 0.5*( 1.0 / tau);
	for(int k = 0; k<9 ; k++)
	{
		double fx,fy,eu;
		eu = (e[k][0]*u.x + e[k][1]*u.y);
		fx = weight[k]*omega*(3.0 * (e[k][0] - u.x) + 9.0 * eu * e[k][0] )* force.x;
		fy = weight[k] * omega * (3.0 * (e[k][1] - u.y) + 9.0 * eu * e[k][1] )* force.y;
		lat_force[k] = fx + fy;
	}
	// lat_force[0] = omega * w0 * (3 * ((0 - u.x) * force.x + (-u.y) * force.y));
	// lat_force[1] = omega * w1 * (3 * ((1 - u.x) * force.x + (-u.y) * force.y) + 9 * (u.x) * force.x);
	// lat_force[2] = omega * w1 * (3 * ((-1 - u.x) * force.x + (-u.y) * force.y) + 9 * (-u.x) * (-force.x));
	// lat_force[3] = omega * w1 * (3 * ((0 - u.x) * force.x + (1 - u.y) * force.y) + 9 * u.y * force.y);
	// lat_force[4] = omega * w1 * (3 * ((0 - u.x) * force.x + (-1 - u.y) * force.y) + 9 * u.y * force.y);
	// lat_force[5] = omega * w2 * (3 * ((1 - u.x) * force.x + (1 - u.y) * force.y) + 9 * (u.x + u.y) * (force.x + force.y));
	// lat_force[6] = omega * w2 * (3 * ((-1 - u.x) * force.x + (-1 - u.y) * force.y) + 9 * (u.x + u.y) * (force.x + force.y));
	// lat_force[7] = omega * w2 * (3 * ((1 - u.x) * force.x + (-1 - u.y) * force.y) + 9 * (u.x - u.y) * (force.x - force.y));
	// lat_force[8] = omega * w2 * (3 * ((-1 - u.x) * force.x + (1 - u.y) * force.y) + 9 * (-u.x + u.y) * (-force.x + force.y));
}
void Lattice::CalMacro()
{
	double trho = 0.0;
	Point2D tu = {0.0, 0.0};
	int e[9][2] = {{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};
	for (int k=0; k<9; k++)
	{
		trho += fnew[k];
		tu.x += fnew[k]*e[k][0];
		tu.y += fnew[k]*e[k][1];
	}
	// tu.x+=0.5*1*1*force.x;
	// tu.y+=0.5*1*1*force.y;
	if (trho < 1.e-6)
	{
		tu.x = 0.0;
		tu.y = 0.0;
	}
	else
	{
		tu.x /= trho;
		tu.y /= trho;
	}
	
	rho = trho;
	u = tu;
	p = rho/3.0;


	
}
#ifdef RRLBM
void Lattice::ComputeRRCollision(){
	Compute_feq();//计算feq
	
	REAL MFTneqXX,MFTneqXY,MFTneqYX,MFTneqYY;
	REAL fone[9];
	REAL mftxx = 0;
	REAL mfttempxx = 0;

	REAL mftxy = 0;
	REAL mfttempxy = 0;

	REAL mftyx= 0;
	REAL mfttempyx = 0;

	REAL mftyy = 0;
	REAL mfttempyy = 0;
	for(int k =0 ; k< 9;k++)
	{
		// mftxx += fnew[k]*e[k][0]*e[k][0];
		// mfttempxx += (feq[k]*e[k][0]*e[k][0] - 1./2 * lat_force[k]);

		// mftxy += fnew[k] * e[k][0] * e[k][1];
		// mfttempxy += (feq[k]* e[k][0] * e[k][1] - 1./2 * lat_force[k]);//xy

		// mftyx += fnew[k] * e[k][1] * e[k][0];
		// mfttempyx += (feq[k] * e[k][1] * e[k][0]- 1./2 * lat_force[k]);//yx

		// mftyy += fnew[k] * e[k][1] * e[k][1];
		// mfttempyy+= (feq[k] * e[k][1] * e[k][1] - 1./2 * lat_force[k]);//yy

		mftxx += fnew[k]*e[k][0]*e[k][0];
		mfttempxx += feq[k]*e[k][0]*e[k][0] ;

		mftxy += fnew[k] * e[k][0] * e[k][1];
		mfttempxy += feq[k]* e[k][0] * e[k][1];//xy

		mftyx += fnew[k] * e[k][1] * e[k][0];
		mfttempyx += feq[k] * e[k][1] * e[k][0];//yx

		mftyy += fnew[k] * e[k][1] * e[k][1];
		mfttempyy+= feq[k] * e[k][1] * e[k][1];//yy
	}
	MFTneqXX = mftxx - mfttempxx;
	MFTneqXY = mftxy - mfttempxy;
	MFTneqYX = mftyx - mfttempyx;
	MFTneqYY = mftyy - mfttempyy;

	for(int k =0; k<9;k++)
	{
		double Q[4];
				Q[0] = e[k][0] * e[k][0] - 1. / 3; //xx
				Q[1] = e[k][0] * e[k][1];//xy
				Q[2] = e[k][1] * e[k][0];//yx
				Q[3] = e[k][1] * e[k][1] - 1. / 3;//yy
		fone[k] = 4.5 *weight[k]*(Q[0]*MFTneqXX + Q[1]* MFTneqXY +Q[2]*MFTneqYX + Q[3]*MFTneqYY);
		fcol[k] = feq[k] + (1.0 - 1.0/tau)*fone[k] ;

	}
}
#endif

#ifdef MRT
void Lattice::ComputeMRTCollision(){
	Compute_feq();
	double fmom[9];
	double fmeq[9];
	double stmiv[9][9];
	// double sm[9] =  { 1.0, 1.4, 1.4, 1.0, 1.2, 1.0, 1.2, tau, tau };
	//double sm[9] = {1.0, 1.4, 1.2, 1.0, 1.0, 1.0, 1.0, tau, tau};
	// double s4 = 8.0*(2.0-1.0/tau)/(8.0 - 1.0/tau);
	// double sm[9] =  { 0, 1./tau, 1./tau, 0, s4,0, s4, 1./tau, 1./tau };



	
	// double sm[9] =  { 0, 1./tau, 1.8, 0, s4,0, s4, 1./tau, 1./tau };


	double sm[9] =  {0, 1.2, 1.1, 0, 1.2,0, 1.2, 1./tau, 1./tau };

	// double sm[9] =  {1, 1.4, 1.4, 1, 1.2,1, 1.2, 1./tau, 1./tau };
	
	for(int i = 0; i <9; i++)
		for(int j =0 ; j <9; j++)
		{
			stmiv[i][j] = tminv[i][j] * sm[j];
		}
	double t1 = u.x*u.x + u.y*u.y;
	double t2 = u.x*u.x - u.y*u.y;

//------------计算平衡时间--------------------//
	fmeq[0] = rho;
	fmeq[1] = rho*(-2.0 + 3.0*t1);
	fmeq[2] = rho*(1.0 - 3.0*t1);
	fmeq[3] = rho*u.x;
	fmeq[4] = -rho*u.x;
	fmeq[5] = rho*u.y;
	fmeq[6] = -rho*u.y;
	fmeq[7] = rho*t2;
	fmeq[8] = rho*u.x*u.y;
	// fmeq[0] = rho;
	// fmeq[1] = -2.0*rho + 3.0*t1;
	// fmeq[2] = rho - 3.0*t1;
	// fmeq[3] = rho*u.x;
	// fmeq[4] = -u.x;
	// fmeq[5] = rho*u.y;
	// fmeq[6] = -u.y;
	// fmeq[7] = t2;
	// fmeq[8] = u.x*u.y;
//--------------------------------------------//
//-----------------计算时间------------------//
	for(int k =0; k <  9 ;k ++)
	{
		double suma = 0.0;
		for(int l =0 ; l<9;l++)
		{
			suma += tm[k][l]*fnew[l];
		}
		fmom[k] = suma;
	}
//-----------------计算时间------------------//
//-----------------计算时间空间中的碰撞------------------//
	for(int k =0; k <  9 ;k ++)
	{
		double sumb = 0.;
		for(int l =0 ; l<9;l++)
		{
			sumb += stmiv[k][l]*(fmom[l] - fmeq[l]);
		}
		fcol[k] = fnew[k]-sumb;
	}	
//-----------------计算时间空间中的碰撞------------------//
}
void Lattice::MRTInitial()
{
	if (feq.empty())
	{
		for (int k = 0; k < 9; k++)
		{
			feq.push_back(0.0);
		}
	}

	if (fnew.empty())
	{
		for (int k = 0; k < 9; k++)
		{
			fnew.push_back(0.0);
		}
	}

	if (fcol.empty())
	{
		for (int k = 0; k < 9; k++)
		{
			fcol.push_back(0.0);
		}
	}

	ComputeMRTCollision();

	// for(int k = 0; k< 9 ; k++)
	// {
	// 	fnew[k] = weight[k] * rho;
	// 	fcol[k] = weight[k] * rho;
	// }

}
#endif
