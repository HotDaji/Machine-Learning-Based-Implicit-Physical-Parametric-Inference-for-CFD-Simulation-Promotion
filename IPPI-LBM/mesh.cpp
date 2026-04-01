#include"mesh.h"

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

//--------------------------xflowmanger---------------------------------//


    struct SaveData
    {
        double ux[9] = {0,0,0,0,0,0,0,0,0}; //宏观量
        double uy[9] = {0,0,0,0,0,0,0,0,0};
        double rho[9] = {0,0,0,0,0,0,0,0,0};
    //中心点
    
        double centerPosX;
        double centerPosY;
        //double centerU,centerV,centerRho;
        
    };

struct hash_pair {
    template <class T1, class T2>
    size_t operator()(const pair<T1, T2>& p) const {
        auto hash1 = hash<T1>{}(p.first);
        auto hash2 = hash<T2>{}(p.second);
        return hash1 ^ hash2;
    }
};

// 设置格式化输出的函数，四舍五入到小数点后5位，忽略后面的数字
void RoundToThreeDecimalPlaces(PointData &point) {
    
    point.vx = std::round(point.vx * 1000.0) / 1000.0;
    point.vy = std::round(point.vy * 1000.0) / 1000.0;
    
}
//--------------------------xflowmanger---------------------------------//

Mesh::Mesh()
{
}
Mesh::~Mesh()
{

}




void Mesh::InitialLatMesh(int Nx_p, int Ny, double tau, Point2D StartVel, double Rho , int RankID , int CoreNum, int PartitionRemainderFlag)
{
   
    int NX;
    if(PartitionRemainderFlag == 0)
    {
        NX = Nx_p + 2 ;//划分长度加2，X方向创建两个缓冲区
    }
    else
    {
        if(RankID == CoreNum - 1)
        {
            NX = PartitionRemainderFlag+2;
        }
        else
        {
            NX = Nx_p + 2 ;
        }
        
    }
    
    LatMesh.resize(NX);

    for (int i = 0; i < NX; i++)
    {
        LatMesh[i].resize(Ny);
    }

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < Ny ; j++) {
            LatMesh[i][j].tau = tau;
            LatMesh[i][j].TauAtBeginning = tau;  
            LatMesh[i][j].rho = Rho;
            LatMesh[i][j].pos.x = i;
            LatMesh[i][j].pos.y = j;
            LatMesh[i][j].force.x = 0.0;
            LatMesh[i][j].force.y = 0.0;
            LatMesh[i][j].u.x = 0.0;
            LatMesh[i][j].u.y = 0.0;
            LatMesh[i][j].Cs = 0.0;
            
        }
    }
    for (int i = 1; i < NX-1; i++) {
        for (int j = 0; j < Ny; j++) {
            LatMesh[i][j].pos.x = (i-1)+Nx_p*RankID;
            LatMesh[i][j].pos.y = j;
        }
    }
    for(int j = 0; j < Ny; j++) //缓冲区也需要坐标 为速度的Interplote和spread提供邻居进程的信息
    {
        LatMesh[0][j].pos.x = -1 + Nx_p*(RankID); //
        LatMesh[NX-1][j].pos.x = Nx_p*(RankID+1);
    }
   
    // if(RankID == 7){
    // cerr<<-1 + Nx_p*(RankID)<<endl;
    // cerr<<LatMesh[1][0].pos.x<<endl;
    // cerr<<LatMesh[NX-2][0].pos.x<<endl;
    // cerr<<Nx_p*(RankID+1)<<endl;
    
    // }
    // if(RankID == 6){
    //     cerr<<"6:"<<LatMesh[NX-2][0].pos.x<<endl;
    // }
    // if(RankID ==8){
    //     cerr<<"8:"<<LatMesh[1][0].pos.x<<endl;
    // }
    
    // if(!RankID){
       
    //     for (int j = 1; j < Ny - 1; j++)
    //     {
    //     LatMesh[1][j].u = StartVel;
                
    //     }                        
        
    // }
    for(int i = 0 ; i < NX; i++)
    {
       for (int j = 0; j < Ny; j++)
        {
        
        LatMesh[i][j].u = {0.0,0.0};  
        
        }
 
    }
        
#ifdef MRT
    for (int i = 0; i < NX; i++)
    {
         for (int j = 0; j < Ny ; j++)
         {
           LatMesh[i][j].InitialLatticeFeq();
           //LatMesh[i][j].MRTInitial();
         }
    }
#else
    for (int i = 0; i < NX; i++)
    {
         for (int j = 0; j < Ny ; j++)
         {
            
            LatMesh[i][j].InitialLatticeFeq();
         }
    }
#endif
    

   
}

void Mesh::InitialIBParticle(int NumNodes, double Radius, double ParticlePositionX, double ParticlePositionY , Point2D ParticleVel)
{
#ifdef D2Q9
   
    IBParticle.resize(NumNodes);
    for (int i = 0; i < NumNodes; i++) 
    {
        IBParticle[i].position.x = ParticlePositionX + Radius * sin(2. * M_PI * (double)i / NumNodes);
        IBParticle[i].position.y = ParticlePositionY + Radius * cos(2. * M_PI * (double)i / NumNodes);
        IBParticle[i].position_ref.x = IBParticle[i].position.x;
        IBParticle[i].position_ref.y = IBParticle[i].position.y;
        IBParticle[i].vel = ParticleVel;
        IBParticle[i].vel_ref.x = 0.0;
        IBParticle[i].vel_ref.y = 0.0;
        IBParticle[i].force.x = 0;
        IBParticle[i].force.y = 0;
    }
#endif // D2Q9
}
void Mesh::Streaming(int RankID, int CoreNum)//D2Q9
{
    int Nx = LatMesh.size();
    int Ny = LatMesh[0].size();
    int ip,jp;


    
    for(int i = 1; i<Nx - 1; i++)
    {
        for(int j = 1; j <Ny - 1  ; j++)
        {
            if(LatMesh[i][j].style=='f')
            {
                for(int k =0; k < 9; k++)
                {
                    ip = i - e[k][0];// x方向一次演化所移动的距离
                    jp = j - e[k][1];// y方向一次演化所移动的距离
                    LatMesh[i][j].fnew[k] = LatMesh[ip][jp].fcol[k];
               

                } 
                LatMesh[i][j].CalMacro();
            }

        }
    }
    

  
    
   
    
    
}

void Mesh::BoundaryCondition(int RankID, int CoreNum)
{

    int Nx = LatMesh.size();
    int Ny = LatMesh[0].size();
    REAL U = startU.x;
     //--------------------固体边界YMS---------------------------------------------------------//
    //const int inv[9]  = {0, 3, 4, 1, 2, 7, 8, 5, 6};
    const int inv[9]  = {0, 3, 4, 1, 2, 7, 8, 5, 6};


    if(RankID==0)//i
    {
        for(int j = 1 ; j <Ny-1;j++)
        {
                LatMesh[1][j].rho = 1.0;
                LatMesh[1][j].u = {0.1,0.0};
                LatMesh[1][j].Compute_feq();
                LatMesh[2][j].Compute_feq();
                for(int k = 0; k < 9 ; k++)
                {
                    LatMesh[1][j].fnew[k] = LatMesh[1][j].feq[k] + LatMesh[2][j].fnew[k] - LatMesh[2][j].feq[k];
                }
        }

        //处理四角点
                // LatMesh[1][0].rho = 1.0;
                // LatMesh[1][0].u = LatMesh[2][1].u;
                // LatMesh[1][0].Compute_feq();
                // LatMesh[2][1].Compute_feq();
                // for(int k = 0; k < 9 ; k++)
                // {
                //     LatMesh[1][0].fnew[k] = LatMesh[1][0].feq[k] + LatMesh[2][1].fnew[k] - LatMesh[2][1].feq[k];
                // }

                // LatMesh[1][Ny-1].rho = 1.0;
                // LatMesh[1][Ny-1].u = LatMesh[2][Ny-2].u;
                // LatMesh[1][Ny-1].Compute_feq();
                // LatMesh[2][Ny-2].Compute_feq();
                // for(int k = 0; k < 9 ; k++)
                // {
                //     LatMesh[1][Ny-1].fnew[k] = LatMesh[1][Ny-1].feq[k] + LatMesh[2][Ny-2].fnew[k] - LatMesh[2][Ny-2].feq[k];
                // }

    }
    if(RankID == CoreNum - 1 )//o
    {
        for(int j = 1 ; j <Ny-1 ;j++)
        {
            LatMesh[Nx-2][j].rho = 1.0 ;
            LatMesh[Nx-2][j].u =  LatMesh[Nx-3][j].u;
            LatMesh[Nx-2][j].Compute_feq();
            LatMesh[Nx-3][j].Compute_feq();
            for(int k = 0; k < 9 ; k++)
            {
                LatMesh[Nx-2][j].fnew[k] = LatMesh[Nx-2][j].feq[k] + LatMesh[Nx-3][j].fnew[k] - LatMesh[Nx-3][j].feq[k];
            }
            
        }

            // LatMesh[Nx-2][0].rho = 1.0 ;
            // LatMesh[Nx-2][0].u =  LatMesh[Nx-3][1].u;
            // LatMesh[Nx-2][0].Compute_feq();
            // LatMesh[Nx-3][1].Compute_feq();
            // for(int k = 0; k < 9 ; k++)
            // {
            //     LatMesh[Nx-2][0].fnew[k] = LatMesh[Nx-2][0].feq[k] + LatMesh[Nx-3][1].fnew[k] - LatMesh[Nx-3][1].feq[k];
            // }

            // LatMesh[Nx-2][Ny-1].rho = 1.0 ;
            // LatMesh[Nx-2][Ny-1].u =  LatMesh[Nx-3][Ny-2].u;
            // LatMesh[Nx-2][Ny-1].Compute_feq();
            // LatMesh[Nx-3][Ny-2].Compute_feq();
            // for(int k = 0; k < 9 ; k++)
            // {
            //     LatMesh[Nx-2][Ny-1].fnew[k] = LatMesh[Nx-2][Ny-1].feq[k] + LatMesh[Nx-3][Ny-2].fnew[k] - LatMesh[Nx-3][Ny-2].feq[k];
            // }
    }

    for(int i = 1;i < Nx-1 ;i++ ) //上边界处理Zouhe
        {
            LatMesh[i][Ny-1].rho =1.0;
            LatMesh[i][Ny-1].u =LatMesh[i][Ny-2].u;
            LatMesh[i][Ny-1].Compute_feq();
            LatMesh[i][Ny-2].Compute_feq();

            for(int k =0; k<9 ; k++)
            {
               
                LatMesh[i][Ny-1].fnew[k] = LatMesh[i][Ny-1].feq[k] + LatMesh[i][Ny-2].fnew[k] - LatMesh[i][Ny-2].feq[k];
            }

        }

        for(int i = 1 ; i < Nx - 1 ; i++) //下边界
        {
            int j=0;
            if(LatMesh[i][j].style=='m')
            {
                LatMesh[i][j].fnew[0] = LatMesh[i][j].fcol[0];
                LatMesh[i][j].fnew[1] = LatMesh[i-1][j].fcol[1];
                LatMesh[i][j].fnew[2] = LatMesh[i][j+1].fcol[4];
                LatMesh[i][j].fnew[3] = LatMesh[i+1][j].fcol[3];
                LatMesh[i][j].fnew[4] = LatMesh[i][j+1].fcol[4];
                LatMesh[i][j].fnew[5] = LatMesh[i-1][j+1].fcol[8];
                LatMesh[i][j].fnew[6] = LatMesh[i+1][j+1].fcol[7];
                LatMesh[i][j].fnew[7] = LatMesh[i+1][j+1].fcol[7];
                LatMesh[i][j].fnew[8] = LatMesh[i-1][j+1].fcol[8];
                LatMesh[i][j].CalMacro();
            }
            if(LatMesh[i][j].style=='d')
            {
                LatMesh[i][j].rho =  LatMesh[i][j+1].rho;
                LatMesh[i][j].u = {0.0,0.0};

                LatMesh[i][j].Compute_feq();
                LatMesh[i][j+1].Compute_feq();

                for(int k =0; k<9 ; k++)
                {
                LatMesh[i][j].fnew[k] = LatMesh[i][j].feq[k] + LatMesh[i][j+1].fnew[k] - LatMesh[i][j+1].feq[k];
                }
            }
        }
    
    
  

}
void Mesh::OneLBMStep(int RankID, int CoreNum,REAL Cs,REAL Pr,int t)
{
   
    int Nx = LatMesh.size();
    int Ny = LatMesh[0].size();
#ifdef SRTLBM
//-------------lOAD nn --------------------------//
        if(MyrankID == 0)
        {
            NNmodel_Type00(t);
            NNmodel_Type01(t);
            NNmodel_Type02(t);
            
        }
        if(MyrankID == MyCoreNum-1)
        {
            NNmodel_Type04(t);
            NNmodel_Type05(t);
            NNmodel_Type06(t);
            
        }
        NNmodel_Type03(t);
        NNmodel_Type07(t);
        MPI_Barrier(MPI_COMM_WORLD); 
//-------------lOAD nn --------------------------//        
        

        
    
     for (int i = 1; i < Nx-1 ; ++i)
        {
            for (int j = 0; j < Ny ; ++j)
            {
                
                ErrorCatch(i,j);
                if( LatMesh[i][j].style!='s'){    
                    if(LatMesh[i][j].NNflag == false)
                    {
                        CaculateLESTau(i , j, Cs, Pr,RankID,CoreNum);//神经网络标志位为负的继续用LES计算
                    }      
                
                }
                LatMesh[i][j].ComputeBGKCollision();

            }
        }
        // if((t>100000))
        // {
        //     NNmodel_Type03(t);
        //     NNmodel_Type07(t);
                    
        // }
       
#endif
#ifdef RRLBM

    for (int i = 1; i < Nx-1 ; ++i)
        {
            for (int j = 0; j < Ny ; ++j)
            {
                if( LatMesh[i][j].style!='s'){
                CaculateLESTau(i , j, Cs, Pr,RankID,CoreNum);
                LatMesh[i][j].ComputeRRCollision();
                }
                
            }
        }
#endif
#ifdef MRT
    for (int i = 1; i < Nx-1 ; ++i)
        {
            for (int j = 0; j < Ny ; ++j)
            {
                ErrorCatch(i,j);
                if( LatMesh[i][j].style!='s'){
                CaculateLESTau(i , j, Cs, Pr,RankID,CoreNum);
                LatMesh[i][j].ComputeMRTCollision();
                }
                
            }
        }
#endif
   
    int cuont = (Ny) *12;
    vector<double> LeftBufSend;
    vector<double> RightBufSend;
    vector<double> LeftBufRec(cuont);
    vector<double> RightBufRec(cuont);
    //左边发送为tag 0，右边接收为0 ，右边发送为1，左边接收为1
   
    for(int j = 0 ; j < Ny ; j++) //数据压栈
    {
        for(int k =0;k<9;k++)
        {
            LeftBufSend.push_back(LatMesh[1][j].fcol[k]);
            RightBufSend.push_back(LatMesh[Nx-2][j].fcol[k]);
        }
        LeftBufSend.push_back(LatMesh[1][j].u.x);
        LeftBufSend.push_back(LatMesh[1][j].u.y);
        LeftBufSend.push_back(LatMesh[1][j].rho);
        RightBufSend.push_back(LatMesh[Nx-2][j].u.x);
        RightBufSend.push_back(LatMesh[Nx-2][j].u.y);
        RightBufSend.push_back(LatMesh[Nx-2][j].rho);


    }
    //压栈的数据包括了feq和rho，ux,uy
  
    MPI_Barrier(MPI_COMM_WORLD); 

   if(RankID == 0){
    MPI_Send(RightBufSend.data(), cuont, MPI_DOUBLE, RankID + 1, 1 , MPI_COMM_WORLD );
    MPI_Recv(LeftBufRec.data(),cuont, MPI_DOUBLE,CoreNum - 1, 1 , MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Send(LeftBufSend.data(), cuont, MPI_DOUBLE, CoreNum - 1, 0 , MPI_COMM_WORLD );
    MPI_Recv(RightBufRec.data(),cuont, MPI_DOUBLE,RankID + 1, 0 , MPI_COMM_WORLD,MPI_STATUS_IGNORE);
   }else if (RankID == CoreNum - 1)
   {
    MPI_Recv(LeftBufRec.data(),cuont, MPI_DOUBLE,RankID - 1, 1 , MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Send(RightBufSend.data(), cuont, MPI_DOUBLE, 0 , 1 , MPI_COMM_WORLD );
    MPI_Recv(RightBufRec.data(),cuont, MPI_DOUBLE,0 , 0 , MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Send(LeftBufSend.data(), cuont, MPI_DOUBLE, RankID - 1, 0 , MPI_COMM_WORLD );
   }else
   {  
    MPI_Recv(LeftBufRec.data(),cuont, MPI_DOUBLE,RankID - 1, 1 , MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Send(RightBufSend.data(), cuont, MPI_DOUBLE, RankID + 1, 1 , MPI_COMM_WORLD );
    MPI_Recv(RightBufRec.data(),cuont, MPI_DOUBLE,RankID + 1, 0 , MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Send(LeftBufSend.data(), cuont, MPI_DOUBLE, RankID - 1, 0 , MPI_COMM_WORLD );
   }

    for(int j = 0 ; j < Ny ; j++) //更新缓冲区
    {
        for(int k =0;k<9;k++)
        {
            LatMesh[0][j].fcol[k] = LeftBufRec[j*12+k];
            LatMesh[Nx-1][j].fcol[k] = RightBufRec[j*12+k];
        }
        LatMesh[0][j].u.x = LeftBufRec[j*12+9];
        LatMesh[0][j].u.y = LeftBufRec[j*12+10];
        LatMesh[0][j].rho = LeftBufRec[j*12+11];
        LatMesh[Nx-1][j].u.x = RightBufRec[j*12+9];
        LatMesh[Nx-1][j].u.y = RightBufRec[j*12+10];
        LatMesh[Nx-1][j].rho = RightBufRec[j*12+11];
    }

    MPI_Barrier(MPI_COMM_WORLD);


    Streaming(RankID , CoreNum);
    
    BoundaryCondition(RankID , CoreNum);
    

    
    MPI_Barrier(MPI_COMM_WORLD); 

   
}
void Mesh::CalParticleForce(int Radius)
{
#ifdef DirectForce
    REAL area = 2 * M_PI * Radius / IBParticle.size();
    
    for (int i = 0; i < IBParticle.size(); i++)
    {
        IBParticle[i].force.x = 0;
        IBParticle[i].force.y = 0;
    }

    for (int i = 0; i < IBParticle.size(); i++) 
    {
        IBParticle[i].force.x = 2 * (IBParticle[i].vel.x - IBParticle[i].vel_ref.x);//*area;
        IBParticle[i].force.y = 2 * (IBParticle[i].vel.y - IBParticle[i].vel_ref.y);//*area;
        // cerr<<IBParticle[i].vel.x<<" "<<IBParticle[i].vel.y<<endl;
    }
    
#endif

#ifdef Stiffness
    REAL area = 2 * M_PI * Radius / IBParticle.size();
    REAL s = 0.1;
    for (int i = 0; i < IBParticle.size(); i++)
    {
        IBParticle[i].force.x = 0;
        IBParticle[i].force.y = 0;
    }

    for (int i = 0; i < IBParticle.size(); i++) 
    {
        IBParticle[i].force.x =  -s * (IBParticle[i].position.x - IBParticle[i].position_ref.x ) * area;
        IBParticle[i].force.y =  -s * ( IBParticle[i].position.y - IBParticle[i].position_ref.y) * area;
        // cerr<<IBParticle[i].vel.x<<" "<<IBParticle[i].vel.y<<endl;
    }
#endif
    
}
void Mesh::InterpolateVel(int RankID, int CoreNum , int Nx_p)
{

#ifdef Dirac_Delta
    int ParticleNodeNum = IBParticle.size();
    int Nx = LatMesh.size();
    int Ny = LatMesh[0].size();

    for (int i = 0; i < ParticleNodeNum; i++)
    {   
        IBParticle[i].vel_ref.x = 0;
        IBParticle[i].vel_ref.y = 0; 
        if((IBParticle[i].position.x<LatMesh[1][0].pos.x )||(IBParticle[i].position.x >= LatMesh[Nx-1][0].pos.x )){ continue;}
        
        int XInt = (int)(IBParticle[i].position.x  );
        int YInt = (int)(IBParticle[i].position.y  );
        XInt= XInt - Nx_p*RankID + 1;
        // if(XInt == Nx-1){XInt = XInt -1;}
        // if(XInt == (Nx-1)){cerr<<"Boundary"<<endl;}
        for (int X = XInt-1 ; X <= XInt+2; X++) 
        {
            for(int Y = YInt-1 ; Y <= YInt+2 ; Y++)
            {
                REAL DistX = IBParticle[i].position.x - LatMesh[X][Y].pos.x ;
                REAL DistY = IBParticle[i].position.y - LatMesh[X][Y].pos.y ;
                REAL DeltaValue = DiracDeltaD2Q9(DistX, DistY);
                // cerr<<"Weight:"<<WeightX<<" "<<WeightY<<endl;
                // cerr<<"DIST:"<<DistX<<" "<<DistY<<endl;
                IBParticle[i].vel_ref.x += (LatMesh[X][Y].u.x * DeltaValue);
                IBParticle[i].vel_ref.y += (LatMesh[X][Y].u.y * DeltaValue);
                //cerr<<"Xvel:"<<LatMesh[X][Y].u.x<<" Yvel:"<<LatMesh[X][Y].u.y<<"RHO:"<<LatMesh[X][Y].rho<<endl;
                
            }
        }

        // for (int X = XInt  ; X < XInt +2; X++) 
        // {
        //     for(int Y = YInt ; Y < YInt + 2 ; Y++)
        //     {
        //         REAL DistX = IBParticle[i].position.x - LatMesh[X][Y].pos.x -0.5 ;
        //         REAL DistY = IBParticle[i].position.y - LatMesh[X][Y].pos.y +0.5;
        //         REAL WeightX = 1 - abs(DistX);
        //         REAL WeightY = 1 - abs(DistY);
        //         // cerr<<"Weight:"<<WeightX<<" "<<WeightY<<endl;
        //         // cerr<<"DIST:"<<DistX<<" "<<DistY<<endl;
        //         IBParticle[i].vel_ref.x += (LatMesh[X][Y].u.x * WeightX * WeightY);
        //         IBParticle[i].vel_ref.y += (LatMesh[X][Y].u.y * WeightX * WeightY);
        //         //cerr<<"Xvel:"<<LatMesh[X][Y].u.x<<" Yvel:"<<LatMesh[X][Y].u.y<<"RHO:"<<LatMesh[X][Y].rho<<endl;
                
        //     }
        // }
        
        // cerr<<IBParticle[i].vel_ref.x<<" "<<IBParticle[i].vel_ref.y<<" "<<endl;
    }
#else
    int ParticleNodeNum = IBParticle.size();
    int Nx = LatMesh.size();
    int Ny = LatMesh[0].size();
    
    for (int i = 0; i < ParticleNodeNum; i++)
    {   

        IBParticle[i].vel_ref.x = 0;
        IBParticle[i].vel_ref.y = 0; 
        if((IBParticle[i].position.x<LatMesh[1][0].pos.x)||(IBParticle[i].position.x >= LatMesh[Nx-1][0].pos.x)){ continue;}
        
        int XInt = (int)(IBParticle[i].position.x );
        int YInt = (int)(IBParticle[i].position.y );
        
        XInt= XInt - Nx_p*RankID+1 ;
       
        for (int X = XInt; X <= XInt + 1; X++) 
        {
            for(int Y = YInt; Y <= YInt + 1 ; Y++)
            {
                REAL DistX = IBParticle[i].position.x - LatMesh[X][Y].pos.x ;
                REAL DistY = IBParticle[i].position.y - LatMesh[X][Y].pos.y ;
                REAL WeightX = 1 - abs(DistX);
                REAL WeightY = 1 - abs(DistY);
                // cerr<<"Weight:"<<WeightX<<" "<<WeightY<<endl;
                // cerr<<"DIST:"<<DistX<<" "<<DistY<<endl;
                IBParticle[i].vel_ref.x += (LatMesh[X][Y].u.x * WeightX * WeightY);
                IBParticle[i].vel_ref.y += (LatMesh[X][Y].u.y * WeightX * WeightY);
                //cerr<<"Xvel:"<<LatMesh[X][Y].u.x<<" Yvel:"<<LatMesh[X][Y].u.y<<"RHO:"<<LatMesh[X][Y].rho<<endl;
                
            }
        }
        
        // cerr<<IBParticle[i].vel_ref.x<<" "<<IBParticle[i].vel_ref.y<<" "<<endl;
    }
    
    
#endif

#ifdef Stiffness
    for (int i = 0; i < ParticleNodeNum; i++){
        IBParticle[i].UpdateParticleBoundary();
    }
    MPI_Barrier(MPI_COMM_WORLD);

#endif
    

}
void Mesh::SpreadForce(int RankID, int CoreNum, int Nx_p)
{
    int NeedSendFlag = 0;
#ifdef Dirac_Delta
    int Nx = LatMesh.size();
    int Ny = LatMesh[0].size();

    int ParticleNodeNum = IBParticle.size();

    for (int X = 0; X < Nx; ++X)//Reset Force
    {
        for (int Y = 1; Y < Ny - 1; ++Y)
        {
            LatMesh[X][Y].force.x = 0;
            LatMesh[X][Y].force.y = 0;
        }
    }
    for (int i = 0; i < ParticleNodeNum; ++i)
    {
       
       if((IBParticle[i].position.x< LatMesh[1][0].pos.x  )||(IBParticle[i].position.x>=LatMesh[Nx-1][0].pos.x )){continue;}
        

        int XInt = (int)(IBParticle[i].position.x );
        
        int YInt = (int)(IBParticle[i].position.y );
        XInt= XInt - Nx_p*RankID + 1;

        for (int X = XInt-1  ; X <=XInt+ 2; ++X)
        {
            for (int Y = YInt-1 ; Y <=YInt + 2; ++Y)
            {
                REAL DistX = IBParticle[i].position.x - LatMesh[X][Y].pos.x ;
                REAL DistY = IBParticle[i].position.y - LatMesh[X][Y].pos.y ;
                //cerr<<DistX<< " "<< DistY<<endl;
                REAL DeltaValue = DiracDeltaD2Q9(DistX,DistY);
                
                LatMesh[X][Y].force.x += (IBParticle[i].force.x * DeltaValue);
                LatMesh[X][Y].force.y += (IBParticle[i].force.y * DeltaValue);
                // if(X == Nx-1){cerr<<"yes:"<<RankID<<endl;
                // cerr<<LatMesh[X][Y].force.x<<" "<<LatMesh[X][Y].force.y<<endl;
                // }
                
            }
           
        }


        //  for (int X = XInt ; X < XInt +2 ; ++X)
        // {
        //     for (int Y = YInt ; Y < YInt + 2; ++Y)
        //     {
        //         REAL DistX = IBParticle[i].position.x - LatMesh[X][Y].pos.x ;
        //         REAL DistY = IBParticle[i].position.y - LatMesh[X][Y].pos.y  ;
        //         REAL WeightX = 1 - abs(DistX);
        //         REAL WeightY = 1 - abs(DistY);
        //         // cerr<<"Weight:"<<WeightX<<" "<<WeightY<<endl;
        //         //cerr<<"DIST:"<<DistX<<" "<<DistY<<endl;make
        //         LatMesh[X][Y].force.x += (IBParticle[i].force.x * WeightX * WeightY);
        //         LatMesh[X][Y].force.y += (IBParticle[i].force.y * WeightX * WeightY);
        //         // if(X == Nx-1){cerr<<"yes:"<<RankID<<endl;
        //         // cerr<<LatMesh[X][Y].force.x<<" "<<LatMesh[X][Y].force.y<<endl;
        //         // }
                
        //     }
           
        // }
        //  if((XInt == Nx - 2)){
        //                 cerr<<"右缓冲区插值:"<<LatMesh[Nx-1][YInt].force.x<<" "<<LatMesh[Nx-1][YInt].force.y
        //                 <<" "<<LatMesh[Nx-1][YInt].pos.x<<" "<<LatMesh[Nx-1][YInt].pos.y
        //                  <<endl;
                        
        //             }
        //     if((XInt == 0)){
        //                         cerr<<"左缓冲区插值:"<<LatMesh[1][YInt].force.x<<" "<<LatMesh[1][YInt].force.y
        //                         <<" "<<LatMesh[1][YInt].pos.x<<" "<<LatMesh[1][YInt].pos.y
        //                         <<endl;

        //                     }
    }

    // for(int j =0 ; j <Ny;j++){//测试force

    //     if((LatMesh[Nx-1][j].force.x!=0)||(LatMesh[Nx-1][j].force.y!=0))
    //     cerr<<"yes"<<endl;
    // }
#else
    int Nx = LatMesh.size();
    int Ny = LatMesh[0].size();

    int ParticleNodeNum = IBParticle.size();

    for (int X = 0; X < Nx ; ++X)//
    {
        for (int Y = 0; Y < Ny ; ++Y)
        {
            LatMesh[X][Y].force.x = 0;
            LatMesh[X][Y].force.y = 0;
        }
    }
    
    for (int i = 0; i < ParticleNodeNum; ++i)//
    {

       if((IBParticle[i].position.x<LatMesh[1][0].pos.x )||(IBParticle[i].position.x>=LatMesh[Nx-1][0].pos.x )){continue;}
        
        int XInt = (int)(IBParticle[i].position.x );
        
        int YInt = (int)(IBParticle[i].position.y );
        XInt= XInt - Nx_p*RankID + 1 ;
        // if((XInt  == 1) || (XInt+1 == Nx-1)){NeedSendFlag = 1;}
        for (int X = XInt; X <= XInt + 1; ++X)
        {
            for (int Y = YInt; Y <= YInt + 1; ++Y)
            {
                REAL DistX = IBParticle[i].position.x - LatMesh[X][Y].pos.x ;
                REAL DistY = IBParticle[i].position.y - LatMesh[X][Y].pos.y ;
                REAL WeightX = 1 - abs(DistX);
                REAL WeightY = 1 - abs(DistY);
                // cerr<<"Weight:"<<WeightX<<" "<<WeightY<<endl;
                //cerr<<"DIST:"<<DistX<<" "<<DistY<<endl;make
                LatMesh[X][Y].force.x += (IBParticle[i].force.x * WeightX * WeightY);
                LatMesh[X][Y].force.y += (IBParticle[i].force.y * WeightX * WeightY);
            }
        }
    }
#endif
    //SPREAD环节过后开始MPI通信
   
    int cuont = Ny*2 ;// 初始化数据栈 force.x,force.y
    vector<REAL> LeftBufSend(cuont);
    vector<REAL> RightBufSend(cuont);
    vector<REAL> LeftBufRec(cuont);
    vector<REAL> RightBufRec(cuont);
    
    for(int j =0 ; j <Ny;j++){//数据压栈
        LeftBufSend[j] = LatMesh[0][j].force.x;
        RightBufSend[j] = LatMesh[Nx-1][j].force.x;
        //cerr<<LatMesh[Nx-1][j].force.y<<endl;
    }
    for(int j =0 ; j <Ny;j++){//数据压栈
        LeftBufSend[j+Ny] = LatMesh[0][j].force.y;
        RightBufSend[j+Ny] = LatMesh[Nx-1][j].force.y;
    }

    if(!RankID){
    MPI_Send(RightBufSend.data(), cuont, MPI_DOUBLE, RankID + 1, 1 , MPI_COMM_WORLD );
    MPI_Recv(LeftBufRec.data(),cuont, MPI_DOUBLE,CoreNum - 1, 1 , MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Send(LeftBufSend.data(), cuont, MPI_DOUBLE, CoreNum - 1, 0 , MPI_COMM_WORLD );
    MPI_Recv(RightBufRec.data(),cuont, MPI_DOUBLE,RankID + 1, 0 , MPI_COMM_WORLD,MPI_STATUS_IGNORE);
   }else if (RankID == CoreNum - 1)
   {
    MPI_Recv(LeftBufRec.data(),cuont, MPI_DOUBLE,RankID - 1, 1 , MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Send(RightBufSend.data(), cuont, MPI_DOUBLE, 0 , 1 , MPI_COMM_WORLD );
    MPI_Recv(RightBufRec.data(),cuont, MPI_DOUBLE,0 , 0 , MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Send(LeftBufSend.data(), cuont, MPI_DOUBLE, RankID - 1, 0 , MPI_COMM_WORLD );
   }else
   {  
    MPI_Recv(LeftBufRec.data(),cuont, MPI_DOUBLE,RankID - 1, 1 , MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Send(RightBufSend.data(), cuont, MPI_DOUBLE, RankID + 1, 1 , MPI_COMM_WORLD );
    MPI_Recv(RightBufRec.data(),cuont, MPI_DOUBLE,RankID + 1, 0 , MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Send(LeftBufSend.data(), cuont, MPI_DOUBLE, RankID - 1, 0 , MPI_COMM_WORLD );
   }
    
    
    for(int j = 0 ; j < Ny ; j++) //装载宏观量
    {
            LatMesh[1][j].force.x += LeftBufRec[j];
            LatMesh[1][j].force.y += LeftBufRec[j+Ny];
            LatMesh[Nx-2][j].force.x += RightBufRec[j];
            LatMesh[Nx-2][j].force.y += RightBufRec[j+Ny];
        
    }

 MPI_Barrier(MPI_COMM_WORLD); 

}
void Mesh::CaculateLESTau(int i, int j,REAL Cs, REAL Pr , int RankID, int CoreNum)
{
 

    // REAL dux_dx = (LatMesh[i+1][j].fnew[1] - LatMesh[i-1][j].fnew[2]) / (2.0 * 1);
    // REAL duy_dy = (LatMesh[i][j+1].fnew[3] - LatMesh[i][j-1].fnew[4]) / (2.0 * 1);
    // REAL dux_dy = (LatMesh[i][j+1].fnew[1] - LatMesh[i][j-1].fnew[1]) / (2.0 * 1);
    // REAL duy_dx = (LatMesh[i+1][j].fnew[3] - LatMesh[i-1][j].fnew[3]) / (2.0 * 1);

    // REAL Sxx = dux_dx;
    // REAL Syy = duy_dy;
    // REAL Sxy = 0.5 * (dux_dy + duy_dx);
    // REAL S = (REAL)sqrt(2.0 * (Sxx * Sxx + Syy * Syy + 2.0 * Sxy * Sxy));
    // // REAL nuSGS = (Cs * dx) * (Cs * dx) * S;
    // REAL nuSGS = (Cs * 1) * (Cs * 1) * S;

    // // double tau_t = (3.0 * (1.0 / (tau + nu_t)) + 0.5) / dt;
    // REAL tauSGS = (0.5 + 3.0 * ((LatMesh[i][j].tau - 0.5)/3.0 + nuSGS))/1.0;//由CFL数决定dt
    // LatMesh[i][j].tau = tauSGS;
#ifdef SmagorinskyByHybrid
    

    // REAL dux_dx = (LatMesh[i+1][j].fnew[1] - LatMesh[i-1][j].fnew[2]) / (2.0 * 1);
    // REAL duy_dy = (LatMesh[i][j+1].fnew[3] - LatMesh[i][j-1].fnew[4]) / (2.0 * 1);
    // REAL dux_dy = (LatMesh[i][j+1].fnew[1] - LatMesh[i][j-1].fnew[1]) / (2.0 * 1);
    // REAL duy_dx = (LatMesh[i+1][j].fnew[3] - LatMesh[i-1][j].fnew[3]) / (2.0 * 1);

    // REAL Sxx = dux_dx;
    // REAL Syy = duy_dy;
    // REAL Sxy = 0.5 * (dux_dy + duy_dx);
    // REAL S = (REAL)sqrt(2.0 * (Sxx * Sxx + Syy * Syy + 2.0 * Sxy * Sxy));
    // // REAL nuSGS = (Cs * dx) * (Cs * dx) * S;
    // REAL nuSGS = (Cs * 1) * (Cs * 1) * S;

    // // double tau_t = (3.0 * (1.0 / (tau + nu_t)) + 0.5) / dt;
    // REAL tauSGS = (0.5 + 3.0 * ((LatMesh[i][j].tau - 0.5)/3.0 + nuSGS))/1.0;//由CFL数决定dt
    // LatMesh[i][j].tau = tauSGS;    
    // int Nx = LatMesh.size();
    // int Ny =LatMesh[0].size();

    
    
//-------------------由微观量计算-------------------------------//
    
    REAL delta = 1.0;
    REAL rhoTemp = 1.0;
    rhoTemp *= LatMesh[i][j].rho;
    LatMesh[i][j].Compute_feq();

    double Tau = startTau;
    //double Tau = LatMesh[i][j].tau;
   

    REAL Sxx = 0.0;
    REAL Sxy = 0.0;
    REAL Syx = 0.0;
    REAL Syy = 0.0;

    for(int k =0; k< 9; k++)
    {
        Sxx += e[k][0] * e[k][0] * (LatMesh[i][j].fnew[k] - LatMesh[i][j].feq[k]);
        Sxy += e[k][0] * e[k][1] * (LatMesh[i][j].fnew[k] - LatMesh[i][j].feq[k]);
        Syy += e[k][1] * e[k][1] * (LatMesh[i][j].fnew[k] - LatMesh[i][j].feq[k]);
    }
    Syx = Sxy;

    Sxx *= -3.0 / (2 * rhoTemp * Tau);
    Sxy *= -3.0 / (2 * rhoTemp * Tau);
    Syx *= -3.0 / (2 * rhoTemp * Tau);
    Syy *= -3.0 / (2 * rhoTemp * Tau);
    

    REAL S2 = Sxx * Sxx + Sxy * Sxy + Syx * Syx + Syy * Syy;
    double Q = sqrt(2.0*S2);
    double g1 = startTau*startTau + 18.0*Cs*Cs*Q/(rhoTemp);
    double vt = (sqrt(g1)+startTau)/2.0;
    //REAL vt = rhoTemp * Cs * Cs * delta * delta * sqrt(2 * S2);
    LatMesh[i][j].tau = vt;

 

//---------------------由速度梯度计算-------------------------//

    
    // double cs      = 0.208;
    // double delta   = 1.0;
    // double rhoTemp = 1.0;

    // rhoTemp *= LatMesh[i][j].rho;
	
    // double Sxx = 0.5 * (GetDuxDx(i, j) + GetDuxDx(i, j));
    // double Sxy = 0.5 * (GetDuxDy(i, j) + GetDuyDx(i, j));
    // double Syx = Sxy;
    // double Syy = 0.5 * (GetDuyDy(i, j) + GetDuyDy(i, j));
	
    // double S2 = Sxx * Sxx + Sxy * Sxy + Syx * Syx + Syy * Syy;
	
    // double vt = rhoTemp * cs * cs * delta * delta * sqrt(2 * S2);
    // double Taueff= startTau + 3 * vt;	
    // LatMesh[i][j].tau = Taueff;
        
//---------------------WALE-------------------------//
    // double Cw    = 0.45;
    // double delta = 1.0;
	
    // double Sxx = 0.5 * (GetDuxDx(i, j) + GetDuxDx(i, j));
    // double Sxy = 0.5 * (GetDuxDy(i, j) + GetDuyDx(i, j));
    // double Syx = Sxy;
    // double Syy = 0.5 * (GetDuyDy(i, j) + GetDuyDy(i, j));
	
    // double SS = Sxx * Sxx + Sxy * Sxy + Syx * Syx + Syy * Syy;
	
    // double gxx = GetDuxDx(i, j);
    // double gxy = GetDuxDy(i, j);
    // double gyx = GetDuyDx(i, j);
    // double gyy = GetDuyDy(i, j);
	
    // double gxx2 = gxx * gxx + gxy * gyx;
    // double gxy2 = gxx * gxy + gxy * gyy;
    // double gyx2 = gyx * gxx + gyy * gyx;
    // double gyy2 = gyx * gxy + gyy * gyy;
	
    // double grr2 = gxx * gxx + gxy * gyx + gyx * gxy + gyy * gyy;
	
    // double Gxx = 0.5 * (gxx2 + gxx2) - grr2 / 3.0;
    // double Gxy = 0.5 * (gxy2 + gyx2);
    // double Gyx = 0.5 * (gyx2 + gxy2);
    // double Gyy = 0.5 * (gyy2 + gyy2) - grr2 / 3.0;
	
    // double GG = Gxx * Gxx + Gxy * Gxy + Gyx * Gyx + Gyy * Gyy;
	
    // double vt = Cw * Cw * delta * delta * pow(GG, 1.5) / ( pow(SS, 2.5) + pow(GG, 1.25) + 1e-12);
	
    // double Taueff= startTau + 3 * vt;	
    // LatMesh[i][j].tau = Taueff;
//-----------------动态CS-------------------------------//
    // LatMesh[i][j].Compute_feq();
	// double SC0 = 3;
	// double a = 6.0;
	// double b = 1.0;
	// double SC = SC0/(1+pow(1.71828,a-b* LatMesh[i][j].normalDis));
	
    // // double SC = 0.206;
    // double delta = 1.0;
    // double rhoTemp = 1.0;
    // rhoTemp *= LatMesh[i][j].rho;
    // double Qxx = 0.0; 
    // double Qxy = 0.0;
    // double Qyx = 0.0;
    // double Qyy = 0.0;
    // double Qzx = 0.0;
    // double Qzy = 0.0;

    
    // for(int k =0; k< 9; k++)
    // {
    //     Qxx += e[k][0] * e[k][0] * (LatMesh[i][j].fnew[k] - LatMesh[i][j].feq[k]);
    //     Qxy += e[k][0] * e[k][1] * (LatMesh[i][j].fnew[k] - LatMesh[i][j].feq[k]);
    //     Qyy += e[k][1] * e[k][1] * (LatMesh[i][j].fnew[k] - LatMesh[i][j].feq[k]);
    // }
    // Qyx = Qxy;
    // double Q = sqrt(2 * (Qxx * Qxx + Qxy * Qxy + Qyx * Qyx + Qyy * Qyy));    
    // double tauT = (pow(startTau * startTau + 18.0 * delta * SC * SC * Q / rhoTemp, 0.5) - startTau) / 2.0;
    // double tauEff = startTau + tauT;
	
    // LatMesh[i][j].tau = tauEff;
//----------------------------微观量liu------------------//
    // REAL delta = 1.0;
    // REAL rhoTemp = 1.0;
    // rhoTemp *= LatMesh[i][j].rho;
    // LatMesh[i][j].Compute_feq();

    // double Tau = startTau;
    // // double Tau = LatMesh[i][j].tau;
   

    // REAL Sxx = 0.0;
    // REAL Sxy = 0.0;
    // REAL Syx = 0.0;
    // REAL Syy = 0.0;

    // for(int k =0; k< 9; k++)
    // {
    //     Sxx += e[k][0] * e[k][0] * (LatMesh[i][j].fnew[k] - LatMesh[i][j].feq[k]);
    //     Sxy += e[k][0] * e[k][1] * (LatMesh[i][j].fnew[k] - LatMesh[i][j].feq[k]);
    //     Syy += e[k][1] * e[k][1] * (LatMesh[i][j].fnew[k] - LatMesh[i][j].feq[k]);
    // }
    // Syx = Sxy;

    // // Sxx *= -3.0 / (2 * rhoTemp * Tau);
    // // Sxy *= -3.0 / (2 * rhoTemp * Tau);
    // // Syx *= -3.0 / (2 * rhoTemp * Tau);
    // // Syy *= -3.0 / (2 * rhoTemp * Tau);
  

    // double Q = Sxx * Sxx + Sxy * Sxy + Syx * Syx + Syy * Syy;
    // double Qave = Tau*Tau + 18*(LatMesh[i][j].Cs*LatMesh[i][j].Cs)*sqrt(Q)/rhoTemp;
    // double S2 = (sqrt(Qave) - Tau)/(6);
    // double vt = S2;
    
    // LatMesh[i][j].tau = Tau + 3*vt;

//------------------CoHerent-Model------------------------//
    // double delta = 1.0;
    // REAL rhoTemp = 1.0;
    // rhoTemp *= LatMesh[i][j].rho;
    // LatMesh[i][j].Compute_feq();

    // double Tau = startTau;
    // //double Tau = LatMesh[i][j].tau;
   

    // REAL Sxx = 0.0;
    // REAL Sxy = 0.0;
    // REAL Syx = 0.0;
    // REAL Syy = 0.0;

    // for(int k =0; k< 9; k++)
    // {
    //     Sxx += e[k][0] * e[k][0] * (LatMesh[i][j].fnew[k] - LatMesh[i][j].feq[k]);
    //     Sxy += e[k][0] * e[k][1] * (LatMesh[i][j].fnew[k] - LatMesh[i][j].feq[k]);
    //     Syy += e[k][1] * e[k][1] * (LatMesh[i][j].fnew[k] - LatMesh[i][j].feq[k]);
    // }
    // Syx = Sxy;

    // // Sxx *= -3.0 / (2 * rhoTemp * Tau);
    // // Sxy *= -3.0 / (2 * rhoTemp * Tau);
    // // Syx *= -3.0 / (2 * rhoTemp * Tau);
    // // Syy *= -3.0 / (2 * rhoTemp * Tau);
    
   
    // REAL S2 = Sxx * Sxx + Sxy * Sxy + Syx * Syx + Syy * Syy;
    // double Ccsm = ComputeDynamicCs(i,j);
    // LatMesh[i][j].CoherentCS = Ccsm;
    // REAL vt =Ccsm * delta * delta * sqrt(2 * S2);
    // LatMesh[i][j].tau = startTau + 3*vt;

    // if(RankID == 0)
    //     {
    //         cerr<<Ccsm<<endl;
    //     }


//-------------------hatfilter----------------------------//
    // int ip,jp;
    // int Nx = LatMesh.size();
    // int Ny = LatMesh[0].size();
    // LatMesh[i][j].Compute_feq();
    // int hatsize = 1;
    // int hatsizefeq = 1;
    // double feq[9];

    // for(int k =0; k < 9; k++)
    // {
    //     ip = i - e[k][0];// x方向一次演化所移动的距离
    //     jp = j - e[k][1];// y方向一次演化所移动的距离
    //     if(!((ip<0)||(ip>Nx-1)||(jp<0)||(jp>Ny-1)))
    //     {
    //         hatsize++;
    //         for(int m =0;m<9;m++)
    //         {
    //             LatMesh[i][j].filteredFnew[m] += LatMesh[ip][jp].fnew[m];
    //         }
            
    //     }
    // } 

    // for(int k = 0; k <9 ; k++)
    // {
    //     LatMesh[i][j].filteredFnew[k] = LatMesh[i][j].filteredFnew[k]/hatsize;
    // }

    // for(int k = 0; k <9 ; k++)
    // {
    //     feq[k] = LatMesh[i][j].feq[k];
    // }

    // for(int k =0; k < 9; k++)
    // {
    //     ip = i - e[k][0];// x方向一次演化所移动的距离
    //     jp = j - e[k][1];// y方向一次演化所移动的距离
    //     if(!((ip<0)||(ip>Nx-1)||(jp<0)||(jp>Ny-1)))
    //     {
    //         LatMesh[ip][jp].Compute_feq();
    //         hatsizefeq++;
    //         for(int m =0;m<9;m++)
    //         {
    //             feq[m] += LatMesh[ip][jp].feq[m];
    //         }
            
    //     }
    // } 

    // for(int k = 0; k <9 ; k++)
    // {
    //     feq[k] = feq[k]/hatsizefeq;
    // }


    // REAL delta = 1.0;
    // REAL rhoTemp = 1.0;
    // rhoTemp *= LatMesh[i][j].rho;
    

    // double Tau = startTau;
    // //double Tau = LatMesh[i][j].tau;
   

    // REAL Sxx = 0.0;
    // REAL Sxy = 0.0;
    // REAL Syx = 0.0;
    // REAL Syy = 0.0;

    // for(int k =0; k< 9; k++)
    // {
    //     Sxx += e[k][0] * e[k][0] * (LatMesh[i][j].filteredFnew[k] - feq[k]);
    //     Sxy += e[k][0] * e[k][1] * (LatMesh[i][j].filteredFnew[k] - feq[k]);
    //     Syy += e[k][1] * e[k][1] * (LatMesh[i][j].filteredFnew[k] - feq[k]);
    // }
    // Syx = Sxy;

    // Sxx *= -3.0 / (2 * rhoTemp * Tau);
    // Sxy *= -3.0 / (2 * rhoTemp * Tau);
    // Syx *= -3.0 / (2 * rhoTemp * Tau);
    // Syy *= -3.0 / (2 * rhoTemp * Tau);
    

    // REAL S2 = Sxx * Sxx + Sxy * Sxy + Syx * Syx + Syy * Syy;

    // REAL vt = rhoTemp * Cs * Cs * delta * delta * sqrt(2 * S2);
    // LatMesh[i][j].tau = startTau + 3*vt;


    
    
#endif

}
void Mesh::InitialFile()
{
    
    int ignore;
    ignore = system("mkdir -p fluid");
    ignore = system("mkdir -p particle");
    ignore = system("mkdir -p cp");
    
}
void Mesh::WriteDatFluid(int t , int RankID,int NX , int CoreNum )
{
   

    int time = t;
    stringstream output_filename;
    output_filename << "fluid/fluid_t" << time <<"_"<<RankID<< ".dat";
    ofstream output_file;
    int Nx = LatMesh.size();
    int Ny = LatMesh[0].size();
    /// Open file
    output_file.open(output_filename.str().c_str());

    /// Write Tecplot header
    if(RankID == 0){
    output_file << "TITLE = \"Simulation of Fluid State at Time " << time << "\"\n";
    output_file << "VARIABLES = \"X\" \"Y\" \"Density\" \"U\" \"V\" \"avgu\" \"avgv\" \"tau\" \"avgRho\" \n";

    /// 这样输出的可以直接合并
    output_file << "ZONE T=\"Time " << time << "\", ";
    output_file << "I=" << Ny << ", J=" << NX  << ", K=1, ";
    output_file << "DATAPACKING=POINT, ZONETYPE=Ordered\n";
    }

    for (int X = 1; X < Nx - 1; ++X)
        {
            for (int Y = 0; Y < Ny ; ++Y)
            {
                output_file << LatMesh[X][Y].pos.x << " " << LatMesh[X][Y].pos.y << " "; // X and Y coordinates
                output_file << LatMesh[X][Y].rho << " ";
                output_file << LatMesh[X][Y].u.x << " ";
                output_file << LatMesh[X][Y].u.y << " ";
                output_file << LatMesh[X][Y].avg_u << " ";
                output_file << LatMesh[X][Y].avg_v << " ";
                output_file << LatMesh[X][Y].tau << " ";
                output_file << LatMesh[X][Y].avg_density <<"\n"; // Density difference
                // output_file <<LatMesh[X][Y].CoherentCS<< "\n";
                // if(LatMesh[X][Y].style=='f')
                // {
                //     output_file <<0<< "\n";
                // }
                // else if(LatMesh[X][Y].style=='s')
                // {
                //     output_file <<1<< "\n";
                // }
                // else if(LatMesh[X][Y].style=='b')
                // {
                //     output_file <<2<< "\n";
                // }
                // else{
                //     output_file <<3<< "\n";
                // }
                //output_file << (int)LatMesh[X][Y].style<< "\n"; // Velocity components
            }
        } 
    

    /// Close file
    output_file.close();

    return;

    
}
void Mesh::UpdateVelocity(int RankID, int CoreNum)
{//只更新除了缓冲区内的所有欧拉点的值
    int Nx = LatMesh.size();
    int Ny = LatMesh[0].size();
    if(RankID == 0)
    {
        for(int i = 2;i < Nx - 1; i++)
        {
            for(int j = 0; j < Ny ; j++)
            {
                LatMesh[i][j].u.x = LatMesh[i][j].u.x + 0.5 * 1 * LatMesh[i][j].force.x/LatMesh[i][j].rho;
                LatMesh[i][j].u.y = LatMesh[i][j].u.y + 0.5 * 1 * LatMesh[i][j].force.y/LatMesh[i][j].rho;
            }
        }
    }
    else if(RankID == CoreNum - 1)
    {
        for(int i = 1;i < Nx - 2; i++)
        {
            for(int j = 0; j < Ny ; j++)
            {
                LatMesh[i][j].u.x = LatMesh[i][j].u.x + 0.5 * 1 * LatMesh[i][j].force.x/LatMesh[i][j].rho;
                LatMesh[i][j].u.y = LatMesh[i][j].u.y + 0.5 * 1 * LatMesh[i][j].force.y/LatMesh[i][j].rho;
            }
        }
    }
    else
    {
        for(int i = 1;i < Nx - 1; i++)
        {
            for(int j = 0; j < Ny ; j++)
            {
                LatMesh[i][j].u.x = LatMesh[i][j].u.x + 0.5 * 1 * LatMesh[i][j].force.x/LatMesh[i][j].rho;
                LatMesh[i][j].u.y = LatMesh[i][j].u.y + 0.5 * 1 * LatMesh[i][j].force.y/LatMesh[i][j].rho;
            }
        }
    }
}

void Mesh::WriteDatParticle(int t)
{
    int time = t;
    stringstream output_filename;
    output_filename << "particle/particle_t" << time << ".dat";
    ofstream output_file;

    /// Open file
    output_file.open(output_filename.str().c_str());

    /// Write Tecplot header
    output_file << "TITLE = \"Particle State at Time " << time << "\"\n";
    output_file << "VARIABLES = \"X\" \"Y\" \"Z\"\n";

    /// Assuming particles form a series of lines or a polyline
    output_file << "ZONE T=\"Particle\", N=" << IBParticle.size() << ", E=" << IBParticle.size();
    output_file << ", DATAPACKING=POINT, ZONETYPE=FELINESEG\n";
    output_file << "Nodes=" << IBParticle.size() << ", Elements=" << IBParticle.size() << "\n";

    /// Write node positions
    for (int n = 0; n < IBParticle.size(); ++n) {
        output_file << IBParticle[n].position.x << " " << IBParticle[n].position.y << " 0\n";
    }

    /// Write connectivity of lines
    for (int n = 0; n < IBParticle.size() - 1; ++n) {
        // Tecplot uses 1-based indexing for elements
        output_file << n + 1 << " " << n + 2 << "\n";
    }
    output_file << 1 << " " << IBParticle.size() << "\n";
    // Optionally, connect the last node to the first for a closed loop
    // output_file << particle_num_nodes << " 1\n";

    /// Close file
    output_file.close();
    return;
}
void Mesh::MergeDatFluid(int t , int CoreNum)
{
    int time = t;
    stringstream output_filename;
    output_filename << "fluid/fluid_Merged_" << time << ".dat";
    stringstream BaseFilename;
    BaseFilename<< "fluid/fluid_t" << time <<"_";

    std::ofstream outFile(output_filename.str());

    for (int i = 0; i <= CoreNum; ++i)
    {
        stringstream oss;
        oss << BaseFilename.str() << i << ".dat";
        std::string inputFilename = oss.str();
        // 打开源文件
        std::ifstream inFile(inputFilename);
        outFile << inFile.rdbuf();
        
        // 关闭源文件
        inFile.close();
        
        // 删除源文件
        if (std::remove(inputFilename.c_str()) != 0) {
            std::cerr << "Error deleting file: " << inputFilename << std::endl;
        }
    }
    outFile.close();
    
    std::cout << "Files merged and original files deleted successfully." << std::endl;



}
REAL Mesh::DiracDelta(REAL x){
    //  if (std::abs(x) >= 2 ) {//2*deltaX,其中deltaX是格子距离
    //     return (REAL)0;
    // } else if((std::abs(x)>=1) && (abs(x)<=2)){
    //     REAL S= 0.5;
    //     S= S - DiracDelta(2. - std::abs(x));
    //     return (REAL)S;
    // } else{
    //     REAL S;
    //     REAL R;
    //     R = 1.0 + 4.0*abs(x) - 4*(x*x);
    //     S = (1.0/8.0)*(3.0 - 2.0*std::abs(x) + sqrt(R));
    // }
    if(abs(x)<1){
        REAL S;
        REAL R;
        R = 1.0 + 4.0*abs(x) - 4*(x*x);
        S = (1.0/8.0)*(3.0 - 2.0*std::abs(x) + sqrt(R));
        return S;
    }
    else if(abs(x)<2)
    {
        REAL S;
        REAL R;
        R = -7.0 + 12.0*abs(x) - 4*(x*x);
        S = (1.0/8.0)*(5.0 - 2.0*std::abs(x) - sqrt(R));
        return S;
    }
    else{
        return 0;
    }
}
REAL Mesh::DiracDeltaD2Q9(REAL x,REAL y ){//

    return DiracDelta(x) * DiracDelta(y);
}

void Mesh::JudgeStyle(Point2D CylinderCenter1, REAL Radius , int RankID ,int CoreNum)
{
#ifdef D2Q9
   
    int Nx = LatMesh.size();
    int Ny = LatMesh[0].size();

    for(int i = 0; i < Nx ;i++)
    {
        for(int j = 0;j < Ny ; j++)
        {
            LatMesh[i][j].style = 'f';
            if((RankID == 0)&&(i == 1))
            {
                LatMesh[i][j].style = 'i';
            }
            else if((RankID == CoreNum - 1)&&( i == Nx - 2))
            {
                LatMesh[i][j].style = 'o';
            }
            else if( j == 0)
            {
                if(LatMesh[i][j].pos.x<400)//这里需要更改如果更换分辨率
                {LatMesh[i][j].style = 'm';}
                else
                {
                    LatMesh[i][j].style = 'd';
                }
            }
            else if ( j == Ny-1)
            {
                LatMesh[i][j].style = 'u';
            }
        }
    }
    

    for(int i = 0; i<Nx;i++)
        for(int j = 0; j < Ny ; j++)
        {
            LatMesh[i][j].rho = 1.0;
            if(LatMesh[i][j].style == 'f')
            {
                LatMesh[i][j].u = {0.0,0.0};
            }
            else
            {
                LatMesh[i][Ny].u = {0.0,0.0};
            }
            
            
        }




    cout<<"Core "<<RankID<<"Judge Style Done!"<<endl;

 
#endif
    
    
   
}
REAL Mesh::GetDis2D(REAL x1, REAL y1, REAL x2, REAL y2)
{
    double distance = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
    distance = (REAL)distance;
    return distance;
}
void Mesh::OutputCd(int RankID,int CoreNum, int time,double Radius)
{
    int oppo[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

    int Nx = LatMesh.size();
    int Ny = LatMesh[0].size();
    double fx = 0.0,fy =0.0;
    for (int i = 1; i < Nx-1; i++)
    {
        for (int j = 1; j < Ny - 1 ; j++)
        {
            if(LatMesh[i][j].style =='b')
            {
                for(int k = 0;k <9;k++)
                {
                    int ip = i + e[k][0];
                    int jp = j + e[k][1];
                    if (LatMesh[ip][jp].style == 's')
                    {
                    fx += e[k][0] * (LatMesh[i][j].fnew[oppo[k]] + LatMesh[i][j].fcol[k]);//动量交换法
                    fy += e[k][1] * (LatMesh[i][j].fnew[oppo[k]] + LatMesh[i][j].fcol[k]);
                    }
                    
                }
            }
        }
    }
    

    // fx = fx / (0.5*0.1*0.1*100);
    // fy = fy / (0.5*0.1*0.1*100);
   
    
    double cd = 0;
    double cl = 0;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&fx, &cd, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    MPI_Allreduce(&fy, &cl, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   
    cd = cd /(0.5*startU.x*startU.x*2*Radius) ; 
    cl = cl /(0.5*startU.x*startU.x*2*Radius) ;
     MPI_Barrier(MPI_COMM_WORLD);
    if(!RankID)
    {
        cerr<<"cd:"<<cd<<" "<<"cl:"<<cl<<endl;  
        int t = time;
        stringstream output_filename;
        output_filename << "particle/CD_CL"<< ".dat";
        ofstream output_file;
        output_file.open(output_filename.str().c_str(),fstream::app);
        output_file<<t<<" "<<cd<<" "<<cl<<endl;
        output_file.close();
    }
    // if(fx!=0 || fy!=0){
    //     cerr<<"cd:"<<fx<<" "<<"cl:"<<fy<<endl;  
    //     int t = time;
    //     stringstream output_filename;
    //     output_filename << "particle/CD_CL"<< ".dat";
    //     ofstream output_file;
    //     output_file.open(output_filename.str().c_str(),fstream::app);
    //     output_file<<t<<" "<<fx<<" "<<fy<<endl;
    //     output_file.close();
    // }
  
}
void Mesh::UpdateBoundary()
{
    for(int i =0; i < IBParticle.size();i++)
    {
        IBParticle[i].position.x += IBParticle[i].vel_ref.x;
        IBParticle[i].position.y += IBParticle[i].vel_ref.y;
    }
}

double Mesh::ComputeDistance( double x0, double y0, double x1, double y1 )
{
    
    double LocateX,LocateY;
    LocateX = center.x;
    LocateY = center.y;

    double x3 = x1 - x0;
	double x4 = x0 - LocateX;
	double y3 = y1 - y0;
	double y4 = y0 - LocateY;

	double a = x3*x3 + y3*y3;
	double b = 2 * ( x3*x4 + y3*y4 );
	double c = x4*x4 +y4*y4 - radius*radius;

    if( b*b - 4*a*c >=0 )
	{
		double m1 = ( -1*b + sqrt( b*b - 4*a*c ) ) / (2*a);
		double m2 = ( -1*b - sqrt( b*b - 4*a*c ) ) / (2*a);
		if( m1>=0-eps && m1<=1+eps )
		{
			return m1;
		}
		else if( m2>=0-eps && m2<=1+eps )
		{
			return m2;
		}
		else
		{
			return 0;
		}
	}
	else
	{
		return 0;
	}
}

void Mesh::OutputCp(int RankID,int CoreNum, int time,double Radius)
{

    int Nx = LatMesh.size();
    int Ny = LatMesh[0].size();
    stringstream CpFileName;
    CpFileName << "cp/cp_"<<time<<"_" << RankID << ".dat";
    ofstream cpfile;
    /// Open file
    cpfile.open(CpFileName.str().c_str());
    for(int i = 1 ; i<Nx - 1; i ++)
    {
        for(int j = 0 ; j <Ny;j++)
        {
            if(LatMesh[i][j].style=='b')
            {
                cpfile<<LatMesh[i][j].pos.x<<" "<<LatMesh[i][j].pos.y<<" "<<LatMesh[i][j].cp<<endl;
            }
        }
    }
    cpfile.close();
    return ;
}
void Mesh::MergeCp(int t, int CoreNum)
{
     
    int time = t;
    stringstream output_filename;
    output_filename << "cp/cp_Merged_" << time << ".dat";
    stringstream BaseFilename;
    BaseFilename<< "cp/cp_" << time <<"_";

    std::ofstream outFile(output_filename.str());

    for (int i = 0; i < CoreNum; ++i)
    {
        stringstream oss;
        oss << BaseFilename.str() << i << ".dat";
        std::string inputFilename = oss.str();
        // 打开源文件
        std::ifstream inFile(inputFilename);

        inFile.seekg(0, std::ios::end);
        if (inFile.tellg() == 0) 
        {
        //std::cerr << "Warning: Input file " << inputFilename << " is empty. Skipping." << std::endl;
        std::remove(inputFilename.c_str());
        inFile.close();
        continue;
        }
        inFile.seekg(0, std::ios::beg);
        outFile << inFile.rdbuf();
        
        // 关闭源文件
        inFile.close();
        
        // 删除源文件
        if (std::remove(inputFilename.c_str()) != 0) {
            std::cerr << "Error deleting file: " << inputFilename << std::endl;
        }
    }
    
    outFile.close();
    
    std::cout << "Files merged and original files deleted successfully." << std::endl;

}
void Mesh::ComputeCp(int startcp,int t)
{
   
    if(t < startcp)
        return;
	int Nx = LatMesh.size();
    int Ny = LatMesh[0].size();
    for(int i = 1; i < Nx-1; i++)
        for(int j = 0; j < Ny; j++)
            if (LatMesh[i][j].style == 'b')
            {
                double temp = (LatMesh[i][j].p  - 1.0 / 3.0) / (0.5*startU.x*startU.x);
                LatMesh[i][j].cp = ( LatMesh[i][j].cp * (t - startcp) + temp) / (t - startcp + 1);
            }
}


double Mesh::GetDuxDx(int i , int j)
{
    int imax = i + 1;
	int Nx = LatMesh.size();
    int Ny = LatMesh[0].size();
	if( ( LatMesh[imax][j].style == 's') ||((MyrankID == MyCoreNum-1)&&(imax==Nx-1)))
	{
		return 0.0;
	}
	else
	{
		
        return (LatMesh[imax][j].u.x - LatMesh[i][j].u.x)/1.0;
	}	
    // if(LatMesh[i][j].u.x>0.0)
    // {
    //     if(isAvailable(i - 1,j))
    //     {
    //         if(isAvailable(i-2,j))
    //         {
    //             return(3.0*LatMesh[i][j].u.x - 4.0*LatMesh[i-1][j].u.x + LatMesh[i-2][j].u.x)/2.0;
    //         }
    //         else
    //         {
    //             return(LatMesh[i][j].u.x - LatMesh[i-1][j].u.x);
    //         }
    //     }
    //     else
    //     {
    //         if(isAvailable(i+1, j))
    //         {
    //             return (LatMesh[i+1][j].u.x - LatMesh[i][j].u.x);
    //         }
    //         else
    //         {
    //             return 0.0;
    //         }
    //     }
    // }
    // else
    // {
    //     if(isAvailable(i + 1,j))
    //     {
    //         if(isAvailable(i+2,j))
    //         {
    //             return -(3.0 * LatMesh[i][j].u.x - 4.0*LatMesh[i+1][j].u.x + LatMesh[i+2][j].u.x)/2.0;
    //         }
    //         else
    //         {
    //             return ( LatMesh[i+1][j].u.x - LatMesh[i][j].u.x);
    //         }
    //     }
    //     else
    //     {
    //         if(isAvailable(i-1 , j))
    //         {
    //             return (LatMesh[i][j].u.x - LatMesh[i-1][j].u.x);
    //         }
    //         else
    //         {
    //             return 0.0;
    //         }
    //     }
    // }

}
double Mesh::GetDuxDy(int i , int j)
{
    int jmax = j + 1;
	int Nx = LatMesh.size();
    int Ny = LatMesh[0].size();

  
	
	if( ((jmax >= Ny-1) || LatMesh[i][jmax].style == 's') )
	{
		return 0.0;
	}
	else
	{
		
        return (LatMesh[i][jmax].u.x - LatMesh[i][j].u.x)/1.0;
	}	

    //  if(LatMesh[i][j].u.y>0.0)
    // {
    //     if(isAvailable(i ,j - 1))
    //     {
    //         if(isAvailable(i,j - 2))
    //         {
    //             return(3.0*LatMesh[i][j].u.x - 4.0*LatMesh[i][j-1].u.x + LatMesh[i][j-2].u.x)/2.0;
    //         }
    //         else
    //         {
    //             return(LatMesh[i][j].u.x - LatMesh[i][j-1].u.x);
    //         }
    //     }
    //     else
    //     {
    //         if(isAvailable(i, j+1))
    //         {
    //             return (LatMesh[i][j+1].u.x - LatMesh[i][j].u.x);
    //         }
    //         else
    //         {
    //             return 0.0;
    //         }
    //     }
    // }
    // else
    // {
    //     if(isAvailable(i ,j+1))
    //     {
    //         if(isAvailable(i,j+2))
    //         {
    //             return -(3.0 * LatMesh[i][j].u.x - 4.0*LatMesh[i][j+1].u.x + LatMesh[i][j+2].u.x)/2.0;
    //         }
    //         else
    //         {
    //             return ( LatMesh[i][j+1].u.x - LatMesh[i][j].u.x);
    //         }
    //     }
    //     else
    //     {
    //         if(isAvailable(i , j-1))
    //         {
    //             return (LatMesh[i][j].u.x - LatMesh[i][j-1].u.x);
    //         }
    //         else
    //         {
    //             return 0.0;
    //         }
    //     }
    // }
	
}
double Mesh::GetDuyDx(int i , int j)
{
    

    int imax = i + 1;
	int Nx = LatMesh.size();
    int Ny = LatMesh[0].size();
	if( ( LatMesh[imax][j].style == 's') ||((MyrankID == MyCoreNum-1)&&(imax==Nx-1)))
	{
		return 0.0;
	}
	else
	{
		
        return (LatMesh[imax][j].u.y - LatMesh[i][j].u.y)/1.0;
	}
    //  if(LatMesh[i][j].u.x>0.0)
    // {
    //     if(isAvailable(i - 1,j))
    //     {
    //         if(isAvailable(i-2,j))
    //         {
    //             return(3.0*LatMesh[i][j].u.y - 4.0*LatMesh[i-1][j].u.y + LatMesh[i-2][j].u.y)/2.0;
    //         }
    //         else
    //         {
    //             return(LatMesh[i][j].u.y - LatMesh[i-1][j].u.y);
    //         }
    //     }
    //     else
    //     {
    //         if(isAvailable(i+1, j))
    //         {
    //             return (LatMesh[i+1][j].u.y - LatMesh[i][j].u.y);
    //         }
    //         else
    //         {
    //             return 0.0;
    //         }
    //     }
    // }
    // else
    // {
    //     if(isAvailable(i + 1,j))
    //     {
    //         if(isAvailable(i+2,j))
    //         {
    //             return -(3.0 * LatMesh[i][j].u.y - 4.0*LatMesh[i+1][j].u.y + LatMesh[i+2][j].u.y)/2.0;
    //         }
    //         else
    //         {
    //             return ( LatMesh[i+1][j].u.y - LatMesh[i][j].u.y);
    //         }
    //     }
    //     else
    //     {
    //         if(isAvailable(i-1 , j))
    //         {
    //             return (LatMesh[i][j].u.y - LatMesh[i-1][j].u.y);
    //         }
    //         else
    //         {
    //             return 0.0;
    //         }
    //     }
    // }



}
double Mesh::GetDuyDy(int i , int j)
{
    int jmax = j + 1;
	int Nx = LatMesh.size();
    int Ny = LatMesh[0].size();

  
	
	if( ((jmax >= Ny - 1) || LatMesh[i][jmax].style == 's') )
	{
		return 0.0;
	}
	else
	{
		
        return (LatMesh[i][jmax].u.y - LatMesh[i][j].u.y)/1.0;
	}

    //   if(LatMesh[i][j].u.y>0.0)
    // {
    //     if(isAvailable(i ,j - 1))
    //     {
    //         if(isAvailable(i,j - 2))
    //         {
    //             return(3.0*LatMesh[i][j].u.y - 4.0*LatMesh[i][j-1].u.y + LatMesh[i][j-2].u.y)/2.0;
    //         }
    //         else
    //         {
    //             return(LatMesh[i][j].u.y - LatMesh[i][j-1].u.y);
    //         }
    //     }
    //     else
    //     {
    //         if(isAvailable(i, j+1))
    //         {
    //             return (LatMesh[i][j+1].u.y - LatMesh[i][j].u.y);
    //         }
    //         else
    //         {
    //             return 0.0;
    //         }
    //     }
    // }
    // else
    // {
    //     if(isAvailable(i ,j+1))
    //     {
    //         if(isAvailable(i,j+2))
    //         {
    //             return -(3.0 * LatMesh[i][j].u.y - 4.0*LatMesh[i][j+1].u.y + LatMesh[i][j+2].u.y)/2.0;
    //         }
    //         else
    //         {
    //             return ( LatMesh[i][j+1].u.y - LatMesh[i][j].u.y);
    //         }
    //     }
    //     else
    //     {
    //         if(isAvailable(i , j-1))
    //         {
    //             return (LatMesh[i][j].u.y - LatMesh[i][j-1].u.y);
    //         }
    //         else
    //         {
    //             return 0.0;
    //         }
    //     }
    // }	

}
void Mesh::SetNormalDis()
{
    int Nx = LatMesh.size();
    int Ny = LatMesh[0].size();

    for(int i = 0; i <Nx ; i++)
        for(int j = 0 ; j <Ny;j++)
        {
            if(LatMesh[i][j].style =='s')
                continue;
            double x = i;
            double y = j;
            double lx = center.x;
            double ly = center.y;
            double len = sqrt( (x-lx)*(x-lx) + (y-ly)*(y-ly) ) - radius;
			len = len/1.0;
            LatMesh[i][j].normalDis = len;
        }
}
void Mesh::ComputeAvgVelocity(int t){

    
    
    if(t<startAvg)
    {
        return;
    }
        

    int Nx = LatMesh.size();
    int Ny = LatMesh[0].size();
    for(int i = 0; i <Nx ; i++)
        for(int j = 0 ; j <Ny;j++)
        {
            LatMesh[i][j].avg_u = ( LatMesh[i][j].avg_u*(t- startAvg)+  LatMesh[i][j].u.x)/(t - startAvg + 1);
            LatMesh[i][j].avg_v = ( LatMesh[i][j].avg_v*(t- startAvg) +  LatMesh[i][j].u.y)/(t - startAvg + 1);
            LatMesh[i][j].avg_density = (LatMesh[i][j].avg_density*(t- startAvg) +  LatMesh[i][j].rho)/(t - startAvg + 1);
        }

}
void Mesh::GetSt(int t){

    int Nx = LatMesh.size();
    int Ny = LatMesh[0].size();
    for(int i = 1; i <Nx-1; i++)
        for(int j = 0 ; j <Ny;j++)
        {
            if(LatMesh[i][j].pos.x == (9.33*stepH))//选择监测点位置（5.33H）
            {
                double yps = LatMesh[i][j].pos.y*0.005;
                
                stringstream output_filename;
                output_filename << "particle/5.33H_"<<t<<".dat";
                ofstream output_file;
                output_file.open(output_filename.str().c_str(),fstream::app);
                
                output_file<<yps<<" "<<LatMesh[i][j].u.x<<" "<<LatMesh[i][j].u.y<<" "<<LatMesh[i][j].avg_u<<" "<<LatMesh[i][j].avg_v<<endl;
                output_file.close();
            }
            
            if(LatMesh[i][j].pos.x == (12*stepH))//选择监测点位置（8H）
            {
                
                double yps2 = LatMesh[i][j].pos.y*0.005;
                stringstream output_filename;
                output_filename << "particle/8.0H_"<<t<<".dat";
                ofstream output_file;
                output_file.open(output_filename.str().c_str(),fstream::app);
                output_file<<yps2<<" "<<LatMesh[i][j].u.x<<" "<<LatMesh[i][j].u.y<<" "<<LatMesh[i][j].avg_u<<" "<<LatMesh[i][j].avg_v<<endl;
                output_file.close();
            }

        }
}
void Mesh::WriteAvgu(int t , int RankID,int NX , int CoreNum )
{
     int time = t;
    stringstream output_filename;
    output_filename << "fluid/Avgu_" << time <<"_"<<RankID<< ".dat";
    ofstream output_file;
    int Nx = LatMesh.size();
    int Ny = LatMesh[0].size();
    /// Open file
    output_file.open(output_filename.str().c_str());

    /// Write Tecplot header
    if(RankID == 0){
    output_file << "TITLE = \"Simulation of Fluid State at Time " << time << "\"\n";
    output_file << "VARIABLES = \"X\" \"Y\" \"AvgdDenstiy\" \"U\" \"V\" \"avgy\" \"avgv\" \"Style\" \n";

    /// 这样输出的可以直接合并
    output_file << "ZONE T=\"Time " << time << "\", ";
    output_file << "I=" << Ny-2 << ", J=" << NX  << ", K=1, ";
    output_file << "DATAPACKING=POINT, ZONETYPE=Ordered\n";
    }

    for (int X = 1; X < Nx - 1; ++X)
        {
            for (int Y = 1; Y < Ny - 1; ++Y)
            {
                output_file << LatMesh[X][Y].pos.x << " " << LatMesh[X][Y].pos.y << " "; // X and Y coordinates
                output_file << LatMesh[X][Y].rho << " ";
                output_file << LatMesh[X][Y].u.x << " ";
                output_file << LatMesh[X][Y].u.y << " ";
                // output_file << LatMesh[X][Y].avg_u << " ";
                // output_file << LatMesh[X][Y].avg_v << " ";
                output_file << LatMesh[X][Y].avg_u << " ";
                output_file << LatMesh[X][Y].avg_v << " "; // Density difference
                if(LatMesh[X][Y].style=='f')
                {
                    output_file <<0<< "\n";
                }
                else if(LatMesh[X][Y].style=='s')
                {
                    output_file <<1<< "\n";
                }
                else if(LatMesh[X][Y].style=='b')
                {
                    output_file <<2<< "\n";
                }
                else{
                    output_file <<3<< "\n";
                }
                //output_file << (int)LatMesh[X][Y].style<< "\n"; // Velocity components
            }
        } 
    

    /// Close file
    output_file.close();

    return;

}
void Mesh::MergeAvgu(int t, int CoreNum)
{
    int time = t;
    stringstream output_filename;
    output_filename << "fluid/fluid_avg_" << time << ".dat";
    stringstream BaseFilename;
    BaseFilename<< "fluid/Avgu_" << time <<"_";

    std::ofstream outFile(output_filename.str());

    for (int i = 0; i <= CoreNum; ++i)
    {
        stringstream oss;
        oss << BaseFilename.str() << i << ".dat";
        std::string inputFilename = oss.str();
        // 打开源文件
        std::ifstream inFile(inputFilename);
        outFile << inFile.rdbuf();
        
        // 关闭源文件
        inFile.close();
        
        // 删除源文件
        if (std::remove(inputFilename.c_str()) != 0) {
            std::cerr << "Error deleting file: " << inputFilename << std::endl;
        }
    }
    outFile.close();
    
    std::cout << "Files merged and original files deleted successfully." << std::endl;

}
void Mesh::ErrorCatch(int i , int j)
{
    if((LatMesh[i][j].rho)>10||(LatMesh[i][j].rho<0))
    {
        cerr<<"Simulation Error!"<<endl;
        errorflag = 1;
        MPI_Abort(MPI_COMM_WORLD, 1);
       
    }

}

double Mesh::ComputeDynamicCs(int i, int j)
{

    double Sxx = 0.5 * (GetDuxDx(i, j) + GetDuxDx(i, j));
    double Sxy = 0.5 * (GetDuxDy(i, j) + GetDuyDx(i, j));
    double Syx = Sxy;
    double Syy = 0.5 * (GetDuyDy(i, j) + GetDuyDy(i, j));
	
    double SS = Sxx * Sxx + Sxy * Sxy + Syx * Syx + Syy * Syy;
	
    double Wxy = 0.5*(GetDuxDy(i,j) - GetDuyDx(i,j));
	
	
    double WW = 2*Wxy*Wxy;

     double CoherentCs = 0.0 ;
    double Fcs;
    double Q = 0.5*(WW - SS);
    double E = 0.5*(WW + SS);
    if(E == 0){
        return 0.016;//原始的系数
    }
    else
    {
        Fcs = Q/E;
    }
    Fcs = abs(Fcs);
    // CoherentCs = (1/22.0)*pow(Fcs,1.5)*(1 - Fcs);
    CoherentCs = (1/11.0)*pow(Fcs,1.5)*(1 - Fcs);
    // if(CoherentCs< 0){
    //     CoherentCs = 0;
    // }
    return CoherentCs;

    // double CoherentCs = 0.0 ;
    // double Fcs;
    // double Q = -0.5*GetDuxDy(i,j)*GetDuyDx(i,j);
    // double E = 0.5*GetDuyDx(i,j)*GetDuyDx(i,j);
    // if(E == 0){
    //     Fcs = 0;//原始的系数
    // }
    // else
    // {
    //     Fcs = Q/E;
    // }
    // Fcs = abs(Fcs);
    // CoherentCs = (1/22.0)*pow(Fcs,1.5)*(1 - Fcs);
    // if(CoherentCs< 0){
    //     CoherentCs = 0;
    // }
    // return CoherentCs;

}
void Mesh::FilterLBM(int i, int j)
{

}
bool Mesh::isAvailable(int i, int j)
{
    int Nx =LatMesh.size();
    int Ny =LatMesh[0].size();


    if(i<0||i>(Nx-1)||j<0||j>(Ny-1) )
    {
        return false;
    }
    if(LatMesh[i][j].style == 's')
    {
        return false;
    }

    return true;
}


void Mesh::ManageDataFromXflow(int t)
{
    
    int step =t+10014;

    int Nx = LatmeshSizeNx;
    int Ny = LatmeshSizeNy;
    string XFlowDataName;
    XFlowDataName = "Xflowdata_pingban_every_1_step/pb_" + std::to_string(step) + ".csv";
    ifstream file(XFlowDataName);
    vector<PointData> data;

    if (!file.is_open()) {
        cerr << "错误: 无法打开文件: " << XFlowDataName << endl;
    }

    string line;

    // 跳过第一行（标题行）
    if (getline(file, line)) {
        // 不需要处理标题行
    }

    // 按行读取数据
    while (getline(file, line)) {
        stringstream ss(line); // 使用字符串流以便于逐个提取数据
        PointData point;       // 创建点数据对象
        
        //id,x,y,z,VelocityModule,StaticPressure,TotalPressure,Velocity_X,Velocity_Y,Velocity_Z,Velocity_Magnitude
        // 尝试读取每个数据项
        char comma; // 用于清除分隔符
        double unused; // 用于存放要跳过的第一列数据
       
        //id,x,y,z,VelocityModule,StaticPressure,TotalPressure,Velocity_X,Velocity_Y,Velocity_Z,Velocity_Magnitude
        if (                     !(ss>>unused>>comma
                                    >> point.vx >> comma 
                                    >> point.vy >> comma 
                                    >> point.vz >> comma 
                                    >> point.VelocityModule>>comma
                                    >> point.StaticPress >> comma
                                    >> unused>>comma 
                                    >> point.TotalPressure >> comma
                                    >> unused>>comma 
                                    >> point.viscosity>>comma  
                                    >> point.ux >> comma 
                                    >> point.uy >> comma 
                                    >> point.uz >> comma 
                                    )) {
            //cerr << "错误: 数据格式无效在行: " << line << endl;
            continue; // 跳过此行
        }
//txt "vtkOriginalPointIds","VelocityModule","StaticPressure","TotalPressure","Velocity:0","Velocity:1","Velocity:2","Points:0","Points:1","Points:2"

        //  if (                     !(ss>>unused>>comma
        //                             >> point.VelocityModule >> comma 
        //                             >> point.StaticPress >> comma 
        //                             >> point.TotalPressure >> comma 
        //                             >> point.ux>>comma
        //                             >> point.uy >> comma 
        //                             >> point.uz >> comma 
        //                             >> point.vx >> comma 
        //                             >> point.vy >> comma 
        //                             >> point.vz >> comma 
        //                             )) {
        //     //cerr << "错误: 数据格式无效在行: " << line << endl;
        //     continue; // 跳过此行
        // }

        // 将有效数据添加到向量中
        RoundToThreeDecimalPlaces(point);
        //cout<<point.vx<<" "<<point.vy<<" "<<point.vz<<endl;
        data.push_back(point);

    }
    file.close();
    if(MyrankID == 0){
        cerr<<data[0].vx<<" "<<data[0].vy<<endl;
    cerr<<data[1].vx<<" "<<data[1].vy<<endl;
    }
    
        //cout<<MyrankID<<" "<<"start Managing!!!!!!!!!!"<<"datasize:"<<data.size()<<endl;
    
    double DX = 0.001;
    double dv = 340.112;
    for ( auto& point : data) {
        
        point.vx = point.vx + 0.4;
        point.vx = point.vx / DX;
        point.vy = point.vy / DX;
        point.vz = point.vz / DX;
        point.VelocityModule = point.VelocityModule/dv/sqrt(3.0);
        point.ux = point.ux/dv/sqrt(3.0);
        point.uy = point.uy/dv/sqrt(3.0);
        point.uz = point.uz/dv/sqrt(3.0);
        point.TotalPressure = (point.StaticPress+101325)/101325;
          
    }
   if(MyrankID == 0)
   {
        cerr<<data[0].vx<<" "<<data[0].vy<<endl;
        cerr<<data[1].vx<<" "<<data[1].vy<<endl;
    }
    unordered_map<pair<double, double>, PointData, hash_pair> combinedData;
     for (const auto& point : data) 
    {
        auto key = make_pair(point.vx, point.vy);

        // 如果已经存在相同的 vx 和 vy, 只进行一次合并，而不是简单的相加
        if (combinedData.find(key) == combinedData.end()) {
            combinedData[key] = point; // 第一次遇到，直接存储
        } else {
            // 当再次遇到相同的 vx 和 vy 时，合并数据并剔除当前项
            combinedData[key].VelocityModule += point.VelocityModule;
            combinedData[key].StaticPress += point.StaticPress;
           
            combinedData[key].TotalPressure += point.TotalPressure;
          
            combinedData[key].ux += point.ux;
            combinedData[key].uy += point.uy;
            combinedData[key].uz += point.uz;
            combinedData[key].vz += point.vz;
            combinedData[key].viscosity += point.viscosity;


            combinedData[key].VelocityModule = combinedData[key].VelocityModule/2.0;
            combinedData[key].StaticPress =  combinedData[key].StaticPress / 2.0;
            
            combinedData[key].TotalPressure = combinedData[key].TotalPressure/2.0;
         
            combinedData[key].ux = combinedData[key].ux/2.0;
            combinedData[key].uy = combinedData[key].uy/2.0;
            combinedData[key].uz = combinedData[key].uz/2.0;
            combinedData[key].vz = combinedData[key].vz/2.0;
            combinedData[key].viscosity = combinedData[key].viscosity/2.0;
            

            // 在这里，我们不需要将当前项保留，因为已合并
        }
    }
    //cerr<<combinedData.size()<<endl;
    // 列出最终结果
    vector<PointData> result;
    
    for (const auto& entry : combinedData) {
        result.push_back(entry.second); // 添加每个合并后的对象
    }
    // if(MyrankID==0)
    // {
    //     stringstream output_filename;
    //     output_filename << "Result.dat";
    //     ofstream output_file;
    //     output_file.open(output_filename.str().c_str());
    //     for(const auto& point : result)
    //     {
    //         output_file<<point.vx<<" "<<point.vy<<endl;
    //     }
    //     output_file.close();
    // }
   

//-----------------------------------剪裁一下result，让它仅保留该进程有的区域--------------//
    int MyThreadSizeStart = MyrankID * ParallelSize;
    int MyThreadSizeEnd = (MyrankID + 1) * ParallelSize+2;

    if(MyrankID == MyCoreNum-1) MyThreadSizeEnd = LatmeshSizeNx;
    if(!MyrankID) MyThreadSizeStart = MyThreadSizeStart - 2;
    
   
    result.erase(std::remove_if(result.begin(), result.end(),
                            [&](const PointData& pd) {
                                // 删除在区间外的数据，即保留在区间内的数据
                                return pd.vx < (MyThreadSizeStart-0.5) || pd.vx > (MyThreadSizeEnd+0.5);
                            }),
             result.end());
    
   

    //cerr<<result.size()<<endl;
   
// -------------------- 最后遍历一下这个进程中的所有lattice----------------//
    
    int NX = LatMesh.size()-1;
    int NY = LatMesh[0].size() - 1;
    
    int temp = 0;
    for(int i = 0; i <=NX  ; i++)
    { 
        for(int j = 0; j <=NY ; j++)
        {
            for(const auto& point : result)//120W
            {   
                double a1,a2;
                a1 = abs(point.vx - LatMesh[i][j].pos.x);
                a2 = abs(point.vy - LatMesh[i][j].pos.y);
                if((a1<DX)&&(a2<DX))
                {

                    LatMesh[i][j].absU.x =  (LatMesh[i][j].u.x - point.ux);
                    LatMesh[i][j].absU.y =  (LatMesh[i][j].u.y - point.uy);
                    LatMesh[i][j].absRHO =  (LatMesh[i][j].rho - point.TotalPressure);

                    LatMesh[i][j].u.x = point.ux;
                    LatMesh[i][j].u.y = point.uy;
                    LatMesh[i][j].rho = point.TotalPressure;
                    LatMesh[i][j].tau = 3*point.viscosity+0.5;
                   
                    temp+=1;
                        
                    
                }
                else{continue;}
            }
            
        }
    }

    //cout<<temp<<endl;
    for(int i = 0; i <=NX  ; i++)
    { 
        for(int j = 0; j <=NY ; j++)
        {
            LatMesh[i][j].InitialLatticeFeq();
        }
    }
}
vector<PointData> Mesh::ReadXflowdata(int t)
{
    int step = t+10014;//这里要根据文件起始做调整
    int Nx = LatmeshSizeNx;
    int Ny = LatmeshSizeNy;
    string XFlowDataName;
    XFlowDataName = "Xflowdata_pingban_every_1_step/pb_" + std::to_string(step) + ".csv";
    ifstream file(XFlowDataName);
    vector<PointData> data;
    if (!file.is_open()) {
        cerr << "错误: 无法打开文件: " << XFlowDataName << endl;
    }
    string line;
    // 跳过第一行（标题行）
    if (getline(file, line)) {
        // 不需要处理标题行
    }
    // 按行读取数据
    while (getline(file, line)) {
        stringstream ss(line); // 使用字符串流以便于逐个提取数据
        PointData point;       // 创建点数据对象
        //id,x,y,z,VelocityModule,StaticPressure,TotalPressure,Velocity_X,Velocity_Y,Velocity_Z,Velocity_Magnitude
        // 尝试读取每个数据项
        char comma; // 用于清除分隔符
        double unused; // 用于存放要跳过的第一列数据
       
        //id,x,y,z,VelocityModule,StaticPressure,TotalPressure,Velocity_X,Velocity_Y,Velocity_Z,Velocity_Magnitude
       if (                     !(ss>>unused>>comma
                                    >> point.vx >> comma 
                                    >> point.vy >> comma 
                                    >> point.vz >> comma 
                                    >> point.VelocityModule>>comma
                                    >> point.StaticPress >> comma
                                    >> unused>>comma 
                                    >> point.TotalPressure >> comma
                                    >> unused>>comma 
                                    >> point.viscosity>>comma  
                                    >> point.ux >> comma 
                                    >> point.uy >> comma 
                                    >> point.uz >> comma 
                                    )) {
            //cerr << "错误: 数据格式无效在行: " << line << endl;
            continue; // 跳过此行
        }

        RoundToThreeDecimalPlaces(point);
        
        data.push_back(point);

    }
    file.close();

     
    double DX = 0.001;
    double dv = 340.112;
    for ( auto& point : data) {
        
        point.vx = point.vx + 0.4;
        point.vx = point.vx / DX;
        point.vy = point.vy / DX;
        point.vz = point.vz / DX;
        point.VelocityModule = point.VelocityModule/dv/sqrt(3.0);
        point.ux = point.ux/dv/sqrt(3.0);
        point.uy = point.uy/dv/sqrt(3.0);
        point.uz = point.uz/dv/sqrt(3.0);
        point.TotalPressure = (point.StaticPress+101325)/101325;
          
    }
    
    
    unordered_map<pair<double, double>, PointData, hash_pair> combinedData;
     for (const auto& point : data) 
    {
        
       
        auto key = make_pair(point.vx, point.vy);

        // 如果已经存在相同的 vx 和 vy, 只进行一次合并，而不是简单的相加
        if (combinedData.find(key) == combinedData.end()) {
            combinedData[key] = point; // 第一次遇到，直接存储
        } else {
            // 当再次遇到相同的 vx 和 vy 时，合并数据并剔除当前项
            combinedData[key].VelocityModule += point.VelocityModule;
            combinedData[key].StaticPress += point.StaticPress;
           
            combinedData[key].TotalPressure += point.TotalPressure;
          
            combinedData[key].ux += point.ux;
            combinedData[key].uy += point.uy;
            combinedData[key].uz += point.uz;
            combinedData[key].vz += point.vz;
            combinedData[key].viscosity += point.viscosity;


            combinedData[key].VelocityModule = combinedData[key].VelocityModule/2.0;
            combinedData[key].StaticPress =  combinedData[key].StaticPress / 2.0;
            
            combinedData[key].TotalPressure = combinedData[key].TotalPressure/2.0;
         
            combinedData[key].ux = combinedData[key].ux/2.0;
            combinedData[key].uy = combinedData[key].uy/2.0;
            combinedData[key].uz = combinedData[key].uz/2.0;
            combinedData[key].vz = combinedData[key].vz/2.0;
            combinedData[key].viscosity = combinedData[key].viscosity/2.0;
    
            // 在这里，我们不需要将当前项保留，因为已合并
        }
    }
    //cerr<<combinedData.size()<<endl;
    // 列出最终结果
    vector<PointData> result;
    
    for (const auto& entry : combinedData) {
        result.push_back(entry.second); // 添加每个合并后的对象
    }
   

//-----------------------------------剪裁一下result，让它仅保留该进程有的区域--------------//
    int MyThreadSizeStart = MyrankID * ParallelSize;
    int MyThreadSizeEnd = (MyrankID + 1) * ParallelSize+2;

    if(MyrankID == MyCoreNum-1) MyThreadSizeEnd = LatmeshSizeNx;
    if(!MyrankID) MyThreadSizeStart = MyThreadSizeStart - 2;
    
   
    result.erase(std::remove_if(result.begin(), result.end(),
                            [&](const PointData& pd) {
                                // 删除在区间外的数据，即保留在区间内的数据
                                return pd.vx < (MyThreadSizeStart-0.5) || pd.vx > (MyThreadSizeEnd+0.5);
                            }),
             result.end());

    return result;         
}
void Mesh::XflowBoundary(int t)
{
   vector<PointData> result =ReadXflowdata(t);
   double DX = 0.001;
// -------------------- 最后遍历一下这个进程中的所有lattice----------------//
    
    int NX = LatMesh.size()-1;
    int NY = LatMesh[0].size() - 1;
 
    for(int i = 0; i <=NX  ; i++)
    { 
        for(int j = 0; j <=NY ; j++)
        {
            if(LatMesh[i][j].style == 'f') continue;
            for(const auto& point : result)//120W
            {   
                double a1,a2;
                a1 = abs(point.vx - LatMesh[i][j].pos.x);
                a2 = abs(point.vy - LatMesh[i][j].pos.y);
                if((a1<DX)&&(a2<DX))
                {
                    LatMesh[i][j].u.x = point.ux;
                    LatMesh[i][j].u.y = point.uy;
                    LatMesh[i][j].rho = point.TotalPressure;
                }
                
            }
            //LatMesh[i][j].InitialLatticeFeq();
            
        }
    }

    //cout<<temp<<endl;
   
}


void Mesh::XflowABS(int t)
{
    vector<PointData> result = ReadXflowdata(t);
    double DX = 0.001;
// -------------------- 最后遍历一下这个进程中的所有lattice----------------//
    
    int NX = LatMesh.size()-1;
    int NY = LatMesh[0].size() - 1;
 
    for(int i = 0; i <=NX  ; i++)
    { 
        for(int j = 0; j <=NY ; j++)
        {
            
            for(const auto& point : result)//120W
            {   
                double a1,a2;
                a1 = abs(point.vx - LatMesh[i][j].pos.x);
                a2 = abs(point.vy - LatMesh[i][j].pos.y);
                if((a1<DX)&&(a2<DX))
                {

                  LatMesh[i][j].absU.x =  (LatMesh[i][j].u.x - point.ux);
                  LatMesh[i][j].absU.y =  (LatMesh[i][j].u.y - point.uy);
                  LatMesh[i][j].absRHO =  (LatMesh[i][j].rho - point.TotalPressure);
                  LatMesh[i][j].xflowtau=point.viscosity*3.0+0.5;
                  LatMesh[i][j].abstau = LatMesh[i][j].tau-LatMesh[i][j].xflowtau;
                    
                }
                
            }
            LatMesh[i][j].InitialLatticeFeq();
            
        }
    }

    //cout<<temp<<endl;
   
}

void Mesh::FindAndWriteTau(int i, int j, int t)
{
    double deltaTau;
    double myfcol[9];
    double DX = 0.001;
    double XflowRho,XflowU,XflowV;
    int inei[9],jnei[9];
    for(int k =0; k < 9; k++)//目标格子和它的可计算邻居的节点编号
    {
        inei[k] = i + e[k][0];
        jnei[k] = j + e[k][1];
    } 

    
    vector<PointData> Xflowdata = ReadXflowdata(t+1);//读取这一时刻的xflowdata
     for(const auto& point : Xflowdata)//120W
            {   
                double a1,a2;
                a1 = abs(point.vx - LatMesh[i][j].pos.x);
                a2 = abs(point.vy - LatMesh[i][j].pos.y);
                if((a1<DX)&&(a2<DX))
                {
                  XflowRho = point.TotalPressure;
                  XflowU = point.ux;
                  XflowV = point.uy;
                    
                }
                
            }
    double bestDeltaTau = 1e10;
    double bestABS = 1e10;
    for(int m = -10000;m<10001;m++)
    {
        deltaTau = (double)m*1e-5;

        
        myfcol[0] = (1.0-1.0/(LatMesh[inei[0]][jnei[0]].tau+deltaTau))*LatMesh[inei[0]][jnei[0]].fnew[0] + LatMesh[inei[0]][jnei[0]].feq[0]/(LatMesh[inei[0]][jnei[0]].tau+deltaTau);
        myfcol[1] = (1.0-1.0/(LatMesh[inei[3]][jnei[3]].tau+deltaTau))*LatMesh[inei[3]][jnei[3]].fnew[1] + LatMesh[inei[3]][jnei[3]].feq[1]/(LatMesh[inei[3]][jnei[3]].tau+deltaTau);
        myfcol[2] = (1.0-1.0/(LatMesh[inei[4]][jnei[4]].tau+deltaTau))*LatMesh[inei[4]][jnei[4]].fnew[2] + LatMesh[inei[4]][jnei[4]].feq[2]/(LatMesh[inei[4]][jnei[4]].tau+deltaTau);
        myfcol[3] = (1.0-1.0/(LatMesh[inei[1]][jnei[1]].tau+deltaTau))*LatMesh[inei[1]][jnei[1]].fnew[3] + LatMesh[inei[1]][jnei[1]].feq[3]/(LatMesh[inei[1]][jnei[1]].tau+deltaTau);
        myfcol[4] = (1.0-1.0/(LatMesh[inei[2]][jnei[2]].tau+deltaTau))*LatMesh[inei[2]][jnei[2]].fnew[4] + LatMesh[inei[2]][jnei[2]].feq[4]/(LatMesh[inei[2]][jnei[2]].tau+deltaTau);
        myfcol[5] = (1.0-1.0/(LatMesh[inei[7]][jnei[7]].tau+deltaTau))*LatMesh[inei[7]][jnei[7]].fnew[5] + LatMesh[inei[7]][jnei[7]].feq[5]/(LatMesh[inei[7]][jnei[7]].tau+deltaTau);
        myfcol[6] = (1.0-1.0/(LatMesh[inei[8]][jnei[8]].tau+deltaTau))*LatMesh[inei[8]][jnei[8]].fnew[6] + LatMesh[inei[8]][jnei[8]].feq[6]/(LatMesh[inei[8]][jnei[8]].tau+deltaTau);
        myfcol[7] = (1.0-1.0/(LatMesh[inei[5]][jnei[5]].tau+deltaTau))*LatMesh[inei[5]][jnei[5]].fnew[7] + LatMesh[inei[5]][jnei[5]].feq[7]/(LatMesh[inei[5]][jnei[5]].tau+deltaTau);
        myfcol[8] = (1.0-1.0/(LatMesh[inei[6]][jnei[6]].tau+deltaTau))*LatMesh[inei[6]][jnei[6]].fnew[8] + LatMesh[inei[6]][jnei[6]].feq[8]/(LatMesh[inei[6]][jnei[6]].tau+deltaTau);

        double trho=0;
        double tu=0;
        double tv=0;
        double TargetAbs=0;
        for (int k=0; k<9; k++)//计算出宏观量
	    {
		trho += myfcol[k];
		tu += myfcol[k]*e[k][0];
		tv += myfcol[k]*e[k][1];
	    }
        TargetAbs = abs(trho - XflowRho)+abs(tu - XflowU)+abs(tv - XflowV);
        if(TargetAbs<bestABS)
        {
            bestABS = TargetAbs;
            bestDeltaTau = deltaTau;
        }

    }

    for(int k=0;k<9;k++)
    {
        LatMesh[inei[k]][jnei[k]].tau += bestDeltaTau;
    }
    // for(int k=1;k<9;k++)
    // {
    //     LatMesh[inei[k]][jnei[k]].tau += bestDeltaTau;
    // }

    if(!MyrankID)
    {
        cerr<<"Timestep:"<<t<<" "<<"BestTau:"<<bestDeltaTau<<endl;
    
        stringstream output_filename;
        output_filename << "cp/BestTau"<< ".dat";
        ofstream output_file;
        output_file.open(output_filename.str().c_str(),fstream::app);
        output_file<<t<<" "<<bestDeltaTau<<" ";
        for(int k = 0; k<9;k++)
        {
            output_file<<LatMesh[inei[k]][jnei[k]].tau<<" ";
        }
        output_file<<endl;
        output_file.close();
    }
}

void Mesh::ChangeTauFromXflow(int t)
{
   
    double DX = 0.001;
    double XflowViscosity;
    int Nx = LatMesh.size();
    int Ny = LatMesh[0].size();
    int LatNum = Nx*Ny;
    vector<PointData> Xflowdata = ReadXflowdata(t);
    double TotalErrorV=0;
    double TempErrorV;
    double TotalErrorR=0;
    double TempErrorR;
    for(int i = 0; i < Nx; i++)
    {
        for(int j = 0 ; j <Ny;j++)
        {
            for(const auto& point : Xflowdata)//120W
            {   
                    double a1,a2;
                    a1 = abs(point.vx - LatMesh[i][j].pos.x);
                    a2 = abs(point.vy - LatMesh[i][j].pos.y);
                    if((a1<DX)&&(a2<DX))
                    {
                    XflowViscosity = point.viscosity;
                    TempErrorV = (LatMesh[i][j].u.x - point.ux)*(LatMesh[i][j].u.x - point.ux)+(LatMesh[i][j].u.y - point.uy)*(LatMesh[i][j].u.y - point.uy);
                    TempErrorR = (LatMesh[i][j].rho - point.TotalPressure)*(LatMesh[i][j].rho - point.TotalPressure);
                    }
                    
            }   
            TotalErrorV += TempErrorV;
            TotalErrorR += TempErrorR;
            LatMesh[i][j].tau = 0.5+3.0*XflowViscosity;
        }
    }
   double VelocityMSE = TotalErrorV/LatNum;
   double RhoMSE =  TotalErrorR/LatNum;
    
        
            stringstream output_filename;
            output_filename << "cp/InletMSE"<< ".dat";
            ofstream output_file;
            output_file.open(output_filename.str().c_str(),fstream::app);
            output_file<<t<<" "<<VelocityMSE<<" "<<RhoMSE;
           
            output_file<<endl;
            output_file.close();
        
  
}
int Mesh::NNmodel_Type03(int t)
{
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(MyrankID, &cpuset);

    pthread_t current_thread = pthread_self();
    pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);

    //at::set_num_threads(1); // 正确设置每个进程的线程数

    const std::string model_path = "/public/home/chaosuan/SHU_24820214/NNmodel/AutoT03.pt";
    torch::Device device(torch::kCPU);
    torch::jit::script::Module module = torch::jit::load(model_path, device);
    torch::NoGradGuard no_grad;
    
    int Nx = LatMesh.size();
    //int j = LatMesh[0].size()-1;应该应用在流体点上
    int j = LatMesh[0].size()-2;
    for(int i = 1;i<Nx-1;i++)
    {
        double input[15] = {
        LatMesh[i+1][j].u.x,
        LatMesh[i+1][j].u.y,
        LatMesh[i+1][j].rho,
        LatMesh[i+1][j-1].u.x,
        LatMesh[i+1][j-1].u.y,
        LatMesh[i+1][j-1].rho,
        LatMesh[i][j-1].u.x,
        LatMesh[i][j-1].u.y,
        LatMesh[i][j-1].rho,
        LatMesh[i-1][j-1].u.x,
        LatMesh[i-1][j-1].u.y,
        LatMesh[i-1][j-1].rho,
        LatMesh[i-1][j].u.x,
        LatMesh[i-1][j].u.y,
        LatMesh[i-1][j].rho,
                    };
        std::vector<float> input_float32(15);
        for (int idx = 0; idx < 15; ++idx) {
            input_float32[idx] = static_cast<float>(input[idx]);
        }
        // 将数组转换为 Tensor
        torch::Tensor tensor_Input = torch::from_blob(input_float32.data(), {15}, torch::kFloat32);
        // 调整形状为 {1, 15}
        tensor_Input = tensor_Input.reshape({1, 15});
        tensor_Input = tensor_Input.to(device);
        torch::jit::IValue output = module.forward({tensor_Input});
        // 将输出转为 Tensor
        torch::Tensor output_tensor = output.toTensor();

        // 提取单个 double 值（Tensor 形状为 {1, 1}）
        double NN_Tau = static_cast<double>(output_tensor.item<float>());

        if(NN_Tau>=0.501)//强制
        {
            LatMesh[i][j].tau = NN_Tau;
            LatMesh[i][j].NNflag = true;
        }
        else
        {LatMesh[i][j].NNflag = false;}
        
    }
    return 0;
    
}
int Mesh::NNmodel_Type07(int t)
{
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(MyrankID, &cpuset);

    pthread_t current_thread = pthread_self();
    pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);

    //at::set_num_threads(1); // 正确设置每个进程的线程数

    const std::string model_path = "/public/home/chaosuan/SHU_24820214/NNmodel/AutoT07.pt";
    torch::Device device(torch::kCPU);
    torch::jit::script::Module module = torch::jit::load(model_path, device);
    torch::NoGradGuard no_grad;
    
    int Nx = LatMesh.size();
    //int j = 0;应该应用在流体点上
    int j = 1;
    for(int i = 1;i<Nx-1;i++)
    {
        double input[15] = {
        LatMesh[i-1][j].u.x,
        LatMesh[i-1][j].u.y,
        LatMesh[i-1][j].rho,
        LatMesh[i-1][j+1].u.x,
        LatMesh[i-1][j+1].u.y,
        LatMesh[i-1][j+1].rho,
        LatMesh[i][j+1].u.x,
        LatMesh[i][j+1].u.y,
        LatMesh[i][j+1].rho,
        LatMesh[i+1][j+1].u.x,
        LatMesh[i+1][j+1].u.y,
        LatMesh[i+1][j+1].rho,
        LatMesh[i+1][j].u.x,
        LatMesh[i+1][j].u.y,
        LatMesh[i+1][j].rho,
                    };
        std::vector<float> input_float32(15);
        for (int idx = 0; idx < 15; ++idx) {
            input_float32[idx] = static_cast<float>(input[idx]);
        }
        // 将数组转换为 Tensor
        torch::Tensor tensor_Input = torch::from_blob(input_float32.data(), {15}, torch::kFloat32);
        // 调整形状为 {1, 15}
        tensor_Input = tensor_Input.reshape({1, 15});
        tensor_Input = tensor_Input.to(device);
        torch::jit::IValue output = module.forward({tensor_Input});
        // 将输出转为 Tensor
        torch::Tensor output_tensor = output.toTensor();

        // 提取单个 double 值（Tensor 形状为 {1, 1}）
        double NN_Tau = static_cast<double>(output_tensor.item<float>());

        if(NN_Tau>=0.501)//强制
        {
            LatMesh[i][j].tau = NN_Tau;
            LatMesh[i][j].NNflag = true;
        }
        else
        {LatMesh[i][j].NNflag = false;}
        
    }
    return 0;
}

int Mesh::NNmodel_Type00(int t)
{
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(MyrankID, &cpuset);

    pthread_t current_thread = pthread_self();
    pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);

    //at::set_num_threads(1); // 正确设置每个进程的线程数

    const std::string model_path = "/public/home/chaosuan/SHU_24820214/NNmodel/AutoT00.pt";
    torch::Device device(torch::kCPU);
    torch::jit::script::Module module = torch::jit::load(model_path, device);
    torch::NoGradGuard no_grad;

        double input[9] = {
        LatMesh[1][1].u.x,
        LatMesh[1][1].u.y,
        LatMesh[1][1].rho,
        LatMesh[2][1].u.x,
        LatMesh[2][1].u.y,
        LatMesh[2][1].rho,
        LatMesh[2][0].u.x,
        LatMesh[2][0].u.y,
        LatMesh[2][0].rho};
       
                    
        std::vector<float> input_float32(9);
        for (int idx = 0; idx < 9; ++idx) {
            input_float32[idx] = static_cast<float>(input[idx]);
        }
        // 将数组转换为 Tensor
        torch::Tensor tensor_Input = torch::from_blob(input_float32.data(), {9}, torch::kFloat32);
        // 调整形状为 {1, 15}
        tensor_Input = tensor_Input.reshape({1, 9});
        tensor_Input = tensor_Input.to(device);
        torch::jit::IValue output = module.forward({tensor_Input});
        // 将输出转为 Tensor
        torch::Tensor output_tensor = output.toTensor();

        // 提取单个 double 值（Tensor 形状为 {1, 1}）
        double NN_Tau = static_cast<double>(output_tensor.item<float>());
        
        int i = 1;
        int j = 0;
        if(NN_Tau>=0.501)//强制
        {
            LatMesh[i][j].tau = NN_Tau;
            LatMesh[i][j].NNflag = true;
        }
        else
        {LatMesh[i][j].NNflag = false;}
        
    
    return 0;
}

int Mesh::NNmodel_Type02(int t)
{
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(MyrankID, &cpuset);

    pthread_t current_thread = pthread_self();
    pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);

    //at::set_num_threads(1); // 正确设置每个进程的线程数

    const std::string model_path = "/public/home/chaosuan/SHU_24820214/NNmodel/AutoT02.pt";
    torch::Device device(torch::kCPU);
    torch::jit::script::Module module = torch::jit::load(model_path, device);
    torch::NoGradGuard no_grad;
        int i = 1;
        int j = LatMesh[0].size()-1;//
        double input[9] = {
        LatMesh[1][j-1].u.x,
        LatMesh[1][j-1].u.y,
        LatMesh[1][j-1].rho,
        LatMesh[2][j-1].u.x,
        LatMesh[2][j-1].u.y,
        LatMesh[2][j-1].rho,
        LatMesh[2][j].u.x,
        LatMesh[2][j].u.y,
        LatMesh[2][j].rho};
       
                    
        std::vector<float> input_float32(9);
        for (int idx = 0; idx < 9; ++idx) {
            input_float32[idx] = static_cast<float>(input[idx]);
        }
        // 将数组转换为 Tensor
        torch::Tensor tensor_Input = torch::from_blob(input_float32.data(), {9}, torch::kFloat32);
        // 调整形状为 {1, 15}
        tensor_Input = tensor_Input.reshape({1, 9});
        tensor_Input = tensor_Input.to(device);
        torch::jit::IValue output = module.forward({tensor_Input});
        // 将输出转为 Tensor
        torch::Tensor output_tensor = output.toTensor();

        // 提取单个 double 值（Tensor 形状为 {1, 1}）
        double NN_Tau = static_cast<double>(output_tensor.item<float>());

       if(NN_Tau>=0.501)//强制
        {
            LatMesh[i][j].tau = NN_Tau;
            LatMesh[i][j].NNflag = true;
        }
        else
        {LatMesh[i][j].NNflag = false;}
        
    
    return 0;
}

int Mesh::NNmodel_Type04(int t)
{
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(MyrankID, &cpuset);

    pthread_t current_thread = pthread_self();
    pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);

    //at::set_num_threads(1); // 正确设置每个进程的线程数

    const std::string model_path = "/public/home/chaosuan/SHU_24820214/NNmodel/AutoT04.pt";
    torch::Device device(torch::kCPU);
    torch::jit::script::Module module = torch::jit::load(model_path, device);
    torch::NoGradGuard no_grad;
        int j = LatMesh[0].size()-1;//501-1
        int i = LatMesh.size()-2;
        
        double input[9] = {
        LatMesh[i-1][j].u.x,
        LatMesh[i-1][j].u.y,
        LatMesh[i-1][j].rho,
        LatMesh[i-1][j-1].u.x,
        LatMesh[i-1][j-1].u.y,
        LatMesh[i-1][j-1].rho,
        LatMesh[i][j-1].u.x,
        LatMesh[i][j-1].u.y,
        LatMesh[i][j-1].rho};
       
                    
        std::vector<float> input_float32(9);
        for (int idx = 0; idx < 9; ++idx) {
            input_float32[idx] = static_cast<float>(input[idx]);
        }
        // 将数组转换为 Tensor
        torch::Tensor tensor_Input = torch::from_blob(input_float32.data(), {9}, torch::kFloat32);
        // 调整形状为 {1, 15}
        tensor_Input = tensor_Input.reshape({1, 9});
        tensor_Input = tensor_Input.to(device);
        torch::jit::IValue output = module.forward({tensor_Input});
        // 将输出转为 Tensor
        torch::Tensor output_tensor = output.toTensor();

        // 提取单个 double 值（Tensor 形状为 {1, 1}）
        double NN_Tau = static_cast<double>(output_tensor.item<float>());

        if(NN_Tau>=0.501)//强制
        {
            LatMesh[i][j].tau = NN_Tau;
            LatMesh[i][j].NNflag = true;
        }
        else
        {LatMesh[i][j].NNflag = false;}
    
    return 0;
}

int Mesh::NNmodel_Type06(int t)
{
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(MyrankID, &cpuset);

    pthread_t current_thread = pthread_self();
    pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);

    //at::set_num_threads(1); // 正确设置每个进程的线程数

    const std::string model_path = "/public/home/chaosuan/SHU_24820214/NNmodel/AutoT06.pt";
    torch::Device device(torch::kCPU);
    torch::jit::script::Module module = torch::jit::load(model_path, device);
    torch::NoGradGuard no_grad;
        int j = 0;//501
        int i = LatMesh.size()-2;//因为缓冲区多减1
        double input[9] = {
        LatMesh[i][1].u.x,
        LatMesh[i][1].u.y,
        LatMesh[i][1].rho,
        LatMesh[i-1][1].u.x,
        LatMesh[i-1][1].u.y,
        LatMesh[i-1][1].rho,
        LatMesh[i-1][0].u.x,
        LatMesh[i-1][0].u.y,
        LatMesh[i-1][0].rho};
       
                    
        std::vector<float> input_float32(9);
        for (int idx = 0; idx < 9; ++idx) {
            input_float32[idx] = static_cast<float>(input[idx]);
        }
        // 将数组转换为 Tensor
        torch::Tensor tensor_Input = torch::from_blob(input_float32.data(), {9}, torch::kFloat32);
        // 调整形状为 {1, 15}
        tensor_Input = tensor_Input.reshape({1, 9});
        tensor_Input = tensor_Input.to(device);
        torch::jit::IValue output = module.forward({tensor_Input});
        // 将输出转为 Tensor
        torch::Tensor output_tensor = output.toTensor();

        // 提取单个 double 值（Tensor 形状为 {1, 1}）
        double NN_Tau = static_cast<double>(output_tensor.item<float>());

       if(NN_Tau>=0.501)//强制
        {
            LatMesh[i][j].tau = NN_Tau;
            LatMesh[i][j].NNflag = true;
        }
        else
        {LatMesh[i][j].NNflag = false;}
        
    
    return 0;
}

int Mesh::NNmodel_Type01(int t)
{
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(MyrankID, &cpuset);

    pthread_t current_thread = pthread_self();
    pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);

    //at::set_num_threads(1); // 正确设置每个进程的线程数

    const std::string model_path = "/public/home/chaosuan/SHU_24820214/NNmodel/AutoT01.pt";
    torch::Device device(torch::kCPU);
    torch::jit::script::Module module = torch::jit::load(model_path, device);
    torch::NoGradGuard no_grad;
    
    
    int i = 1;//0是缓冲区
    int Ny = LatMesh[0].size()-1;//501-1 =500
    for(int j = 1;j<Ny;j++)
    {
        double input[15] = {
        LatMesh[i][j+1].u.x,
        LatMesh[i][j+1].u.y,
        LatMesh[i][j+1].rho,
        LatMesh[i+1][j+1].u.x,
        LatMesh[i+1][j+1].u.y,
        LatMesh[i+1][j+1].rho,
        LatMesh[i+1][j].u.x,
        LatMesh[i+1][j].u.y,
        LatMesh[i+1][j].rho,
        LatMesh[i+1][j-1].u.x,
        LatMesh[i+1][j-1].u.y,
        LatMesh[i+1][j-1].rho,
        LatMesh[i][j-1].u.x,
        LatMesh[i][j-1].u.y,
        LatMesh[i][j-1].rho,
                    };
        std::vector<float> input_float32(15);
        for (int idx = 0; idx < 15; ++idx) {
            input_float32[idx] = static_cast<float>(input[idx]);
        }
        // 将数组转换为 Tensor
        torch::Tensor tensor_Input = torch::from_blob(input_float32.data(), {15}, torch::kFloat32);
        // 调整形状为 {1, 15}
        tensor_Input = tensor_Input.reshape({1, 15});
        tensor_Input = tensor_Input.to(device);
        torch::jit::IValue output = module.forward({tensor_Input});
        // 将输出转为 Tensor
        torch::Tensor output_tensor = output.toTensor();

        // 提取单个 double 值（Tensor 形状为 {1, 1}）
        double NN_Tau = static_cast<double>(output_tensor.item<float>());

        if(NN_Tau>=0.501)//强制
        {
            LatMesh[i][j].tau = NN_Tau;
            LatMesh[i][j].NNflag = true;
        }
        else
        {LatMesh[i][j].NNflag = false;}
        
    }
    return 0;
    
}

int Mesh::NNmodel_Type05(int t)
{
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(MyrankID, &cpuset);

    pthread_t current_thread = pthread_self();
    pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);

    //at::set_num_threads(1); // 正确设置每个进程的线程数

    const std::string model_path = "/public/home/chaosuan/SHU_24820214/NNmodel/AutoT05.pt";
    torch::Device device(torch::kCPU);
    torch::jit::script::Module module = torch::jit::load(model_path, device);
    torch::NoGradGuard no_grad;
    
    
    int i = LatMesh.size()-2;//LMH-1是缓冲区
    int Ny = LatMesh[0].size()-1;//501-1 =500
    for(int j = 1;j<Ny;j++)
    {
        double input[15] = {
        LatMesh[i][j+1].u.x,
        LatMesh[i][j+1].u.y,
        LatMesh[i][j+1].rho,
        LatMesh[i-1][j+1].u.x,
        LatMesh[i-1][j+1].u.y,
        LatMesh[i-1][j+1].rho,
        LatMesh[i-1][j].u.x,
        LatMesh[i-1][j].u.y,
        LatMesh[i-1][j].rho,
        LatMesh[i-1][j-1].u.x,
        LatMesh[i-1][j-1].u.y,
        LatMesh[i-1][j-1].rho,
        LatMesh[i][j-1].u.x,
        LatMesh[i][j-1].u.y,
        LatMesh[i][j-1].rho,
                    };
        std::vector<float> input_float32(15);
        for (int idx = 0; idx < 15; ++idx) {
            input_float32[idx] = static_cast<float>(input[idx]);
        }
        // 将数组转换为 Tensor
        torch::Tensor tensor_Input = torch::from_blob(input_float32.data(), {15}, torch::kFloat32);
        // 调整形状为 {1, 15}
        tensor_Input = tensor_Input.reshape({1, 15});
        tensor_Input = tensor_Input.to(device);
        torch::jit::IValue output = module.forward({tensor_Input});
        // 将输出转为 Tensor
        torch::Tensor output_tensor = output.toTensor();

        // 提取单个 double 值（Tensor 形状为 {1, 1}）
        double NN_Tau = static_cast<double>(output_tensor.item<float>());

        if(NN_Tau>=0.501)//强制
        {
            LatMesh[i][j].tau = NN_Tau;
            LatMesh[i][j].NNflag = true;
        }
        else
        {LatMesh[i][j].NNflag = false;}
        
    }
    return 0;
    
}