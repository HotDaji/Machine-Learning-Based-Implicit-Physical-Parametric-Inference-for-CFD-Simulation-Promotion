#include "mesh.h"
using namespace std;
int main(int argc, char* argv[])
{

	int RankID,CoreNum;//ProcessID and number of Processes

   	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&RankID);
   	MPI_Comm_size(MPI_COMM_WORLD,&CoreNum);

	REAL Cs = 0.173;//les系数
    REAL Pr = 1;
    
	

	// int Nx = 2401;
	// int Ny = 501;
	//不同分辨率测试
	int Nx = 2401;
	int Ny = 501;

	


	int StepH = 0;//台阶的长、高
	int StepL = 0;
	const int ParticleNumNodes = 256;
	const REAL Radius = 20;//半径
	const REAL ParticlePositionX = 400;
	const REAL ParticlePositionY = 620;
	Point2D CylinderCenter;
	CylinderCenter.x = ParticlePositionX;
	CylinderCenter.y = ParticlePositionY;
	Point2D ParticleVel;
	ParticleVel.x = 0.0;
	ParticleVel.y = 0.0;
	double Re = 5e5;
    Point2D StartVel;
	StartVel.x = 0.1;
	StartVel.y = 0.0;
	const REAL Rho = 1.0;
	const REAL tau = 0.5 + 6.0*((Ny* StartVel.x) / Re);
	
	const int NumStep = 200000;
	const int SaveStep = 20000;
	const int startcp = 30000;
	const int ManageDataTime = 10000;
	
	
	
	Mesh mesh;
	
	mesh.startAvg = 30000;
	mesh.center = CylinderCenter;
	mesh.radius =Radius;
	mesh.startTau = tau;
	mesh.startU= StartVel;
	mesh.MyrankID = RankID;
	mesh.MyCoreNum = CoreNum;
	mesh.stepH = StepH;
	mesh.stepL = StepL;
	mesh.StartCs = Cs;
	mesh.LatmeshSizeNx = Nx;
	mesh.LatmeshSizeNy = Ny;
	if(RankID == 0)
	{
		mesh.InitialFile();
	}
	
	int Nx_p;
	int PartitionRemainderFlag = 0;
	
	if(Nx % CoreNum == 0) //Nx正好可以被等距划分
	{	
		Nx_p = Nx / CoreNum;//每个进程X方向划分的长度
	}
	else
	{
		Nx_p = Nx / CoreNum;
		PartitionRemainderFlag = Nx % CoreNum;
		PartitionRemainderFlag = Nx_p + PartitionRemainderFlag;//最后一个进程在加上余出的内容
	}
	

	mesh.InitialLatMesh(Nx_p , Ny , tau, StartVel, Rho, RankID, CoreNum , PartitionRemainderFlag); //并行生成网格
	MPI_Barrier(MPI_COMM_WORLD);
	mesh.JudgeStyle(CylinderCenter,Radius,RankID,CoreNum);
	mesh.ParallelSize = Nx_p;
	MPI_Barrier(MPI_COMM_WORLD);
	if(RankID == 0)
	{
	cout << "Initial Completed" << endl;
	cout << "simulation parameters" << endl;
	cout << "=====================" << endl;
	cout << "L = " << StepL << endl;
	cout << "U = " << StartVel.x << endl;
	cout << "Re = " << Re << endl;
	cout << "Tau = " << tau << endl;
	cout << endl;
	cout << "Start Computing" << endl;
	}
	

	for (int t = 1; t <= NumStep; ++t)
	{

		
		mesh.OneLBMStep(RankID , CoreNum , Cs, Pr,t);
		MPI_Barrier(MPI_COMM_WORLD);
		
		if(t%SaveStep==0)
		{
			mesh.WriteDatFluid(t ,RankID,Nx,CoreNum);
			MPI_Barrier(MPI_COMM_WORLD);
			if(RankID == 0){
				mesh.MergeDatFluid(t,CoreNum);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		
		// if((t%100==0)&&(t>100000))
		// {
		// 		mesh.WriteDatFluid(t ,RankID,Nx,CoreNum);
		// 		MPI_Barrier(MPI_COMM_WORLD);
		// 		if(RankID == 0){
		// 			mesh.MergeDatFluid(t,CoreNum);
		// 		}
		// 		MPI_Barrier(MPI_COMM_WORLD);
		// }
		
		
		//mesh.XflowBoundary(t);
		//MPI_Barrier(MPI_COMM_WORLD);
		// if(t == 1)
		// {
		// 	mesh.ManageDataFromXflow(t);
		//  	MPI_Barrier(MPI_COMM_WORLD);
		// }
		// if(t%100==0)
		// {
		// 	mesh.XflowABS(t);
		// 	MPI_Barrier(MPI_COMM_WORLD);

		// 	mesh.WriteDatFluid(t ,RankID,Nx,CoreNum);
		// 	MPI_Barrier(MPI_COMM_WORLD);
		// 	if(RankID == 0){
		// 		mesh.MergeDatFluid(t,CoreNum);
		// 	}
		// 	MPI_Barrier(MPI_COMM_WORLD);

		// }
		// if((t==1) )
		// {
		// 	mesh.ManageDataFromXflow(t);
		// 	MPI_Barrier(MPI_COMM_WORLD);
		// 	if(RankID == 0) cout<<"ManageSuccess in Time Step !!!!!!!!!!!!!!!!!!! "<<t<<endl;


		// 	mesh.WriteDatFluid(t ,RankID,Nx,CoreNum);
		// 	MPI_Barrier(MPI_COMM_WORLD);
		// 	if(RankID == 0){
		// 		mesh.MergeDatFluid(t,CoreNum);
		// 	}
		// 	MPI_Barrier(MPI_COMM_WORLD);
			
			

		// }

		// if((t==1000) )
		// {
		// 	mesh.XflowABS(t);
		// 	MPI_Barrier(MPI_COMM_WORLD);
		// 	if(RankID == 0) cout<<"ManageSuccess in Time Step !!!!!!!!!!!!!!!!!!! "<<t<<endl;


		// 	mesh.WriteDatFluid(t ,RankID,Nx,CoreNum);
		// 	MPI_Barrier(MPI_COMM_WORLD);
		// 	if(RankID == 0){
		// 		mesh.MergeDatFluid(t,CoreNum);
		// 	}
		// 	MPI_Barrier(MPI_COMM_WORLD);
			
			

		// }

		// if((t!=0)&&(t%10 == 0))
		// {
			
		// 		mesh.WriteDatFluid(t ,RankID,Nx,CoreNum);
		// 		MPI_Barrier(MPI_COMM_WORLD);

		// 		if(RankID == 0){
		// 			mesh.MergeDatFluid(t,CoreNum);
		// 		}
		// 		MPI_Barrier(MPI_COMM_WORLD);
		// }
		// mesh.ComputeAvgVelocity(t);
		// MPI_Barrier(MPI_COMM_WORLD);


		
		// if(t%SaveStep==0)
		// {
		// 	mesh.WriteDatFluid(t ,RankID,Nx,CoreNum);
		// 		MPI_Barrier(MPI_COMM_WORLD);

		// 		if(RankID == 0){
		// 			mesh.MergeDatFluid(t,CoreNum);
		// 		}
		// 		MPI_Barrier(MPI_COMM_WORLD);
		// }
		mesh.ComputeAvgVelocity(t);
		MPI_Barrier(MPI_COMM_WORLD);

		if(RankID == 0)
		{
			cout << "Compute Step:"<<" " << t << endl;
		}

	}
	

	MPI_Finalize();
	return 0;

}