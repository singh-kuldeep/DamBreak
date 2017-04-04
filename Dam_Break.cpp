#include "iostream"
#include <vector>
#include <fstream>
#include <math.h>
using namespace std ;
// Calculating Euler Flux

void getEulerFFlux(vector<float> & U, vector<float> & FEuler)
{	
	float h = U[0];
	float u = U[1]/h;
	float v = U[2]/h;
	
	FEuler[0] = h*u;
	FEuler[1] = h*u*u + 0.5*9.8*h*h;
	FEuler[2] = h*u*v;
}

void getEulerGFlux(vector<float> & U, vector<float> & GEuler)
{	
	float h = U[0];
	float u = U[1]/h;
	float v = U[2]/h;
	
	GEuler[0] = h*v;
	GEuler[1] = h*u*v;
	GEuler[2] = h*v*v + 0.5*9.8*h*h;
}

// Use the Roe approximate Riemann solver to calculate Fluxes.
void getRoeFlux(vector<float> & UL, vector<float> & UR, vector<float> & FRoe)
{	
	float g = 9.8; // gravitational constant
	// Primitive and other variables.
    // Left state
    float hL = UL[0];
	float uL = UL[1]/hL;
	float vL = UL[2]/hL;
	float aL = sqrt(0.5*g*hL);
	
	// Right state
    float hR = UR[0];
	float uR = UR[1]/hR;
	float vR = UR[2]/hR;
	float aR = sqrt(0.5*g*hR);
	
	//Roe Averages
	float RT = sqrt(hR/hL);
	float uInt = (uL+RT*uR)/(1+RT);
	float vInt = (vL+RT*vR)/(1+RT);
	float aInt = sqrt(0.5*g*(hL+hR));

	// Differences in primitive variables
	std::vector<float> du(3);	
	du[0] = hR - hL;
	du[1] = uR*hR - uL*hL;
	du[2] = vR*hR - vL*hL;

	// Wave strangth (Characterstic variables)
	float dalpha[3];
	dalpha[0] = (du[0]*(uInt+aInt)-du[1])/(2.0*aInt); 
	dalpha[1] = du[2] - vInt*du[0];
	dalpha[2] = (-du[0]*(uInt - aInt)+du[1])/(2.0*aInt);

	// Absolute values of the wave speeds (Eigenvalues)
	float ws[3];
	ws[0] = fabs(uInt - aInt);
	ws[1] = fabs(uInt);
	ws[2] = fabs(uInt + aInt);

	// Modified wave speeds for nonlinear fields (the so-called entropy fix,
	// which is often implemented to remove non-physical expansion shocks).
    // There are various ways to implement the entropy fix. This is one of them
	
	#if 1
	float Da = max(0.0, 4.0*((uR-aR)-(uL-aL)));
	if (ws[0] < 0.5*Da)
	{
		ws[0] = ws[0]*ws[0]/Da + 0.25*Da ;
	}
	Da = max(0.0, 4.0*((uR+aR)-(uL+aL)));
	if (ws[2]<0.5*Da)
	{
		ws[2] = ws[2]*ws[2]/Da + 0.25*Da;
	}
	#endif

	// Right engenvectors
	float R[3][3];

	R[0][0] = 1.0;
	R[1][0] = uInt - aInt;
	R[2][0] = vInt;

	R[0][1] = 0;
	R[1][1] = 0;
	R[2][1] = 1.0;

	R[0][2] = 1.0;
	R[1][2] = uInt + aInt;
	R[2][2] = vInt;

	// Add the matrix dissipation term to complete the Roe Flux
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			// FRoe[i] = FRoe[i] - 0.5*ws[j]*dU[j]*R[i][j];
			FRoe[i] = 0.5*ws[j]*du[j]*R[i][j];
		}
	}

}

// // CFL condition to maintain stability
// float getTimeStep(float CFL, float dx, int Nx, vector< vector<float> > & U)
// {
// 	float maxSpeed = -1.0;
// 	vector<float> W(3);
	
// 	float u;
// 	float a;
// 	float lembdaMax;

// 	for (int i = 1; i < Nx-1; ++i)
// 	{
// 		U2W(U[i],W);
// 		u = W[1];
// 		a = sqrt(W[4]*W[2]/W[0]);
// 		lembdaMax = fabs(u)+a;
// 		maxSpeed = max(maxSpeed,lembdaMax);	
// 	}
// 	return CFL*dx/maxSpeed; // CFL condition
// }

// Update solution

void netfluxatXinterface(std::vector<float> UE,std::vector<float> UW, std::vector<float> & Fatface)
{
	float hE = UE[0];
	float uE = UE[1];
	float vE = UE[2];

	float hW = UW[0];
	float uW = UW[1];
	float vW = UW[2];

	std::vector<float> FRoe(3);

	getRoeFlux(UE,UW,FRoe);

	// First the Euler part
	Fatface[0] = (vE + vW - FRoe[0])/2 ; //h*v;
	Fatface[1] = (hE*uE*vE + hW*uW*vW - FRoe[1])/2;//h*u*v;
	Fatface[2] = (hE*vE*vE + 0.5*9.8*hE*hE + hW*vW*vW + 0.5*9.8*hW*hW - FRoe[2])/2; //h*v*v + 0.5*9.8*h*h;

}

void netfluxatYinterface(std::vector<float> US,std::vector<float> UN, std::vector<float> & Gatface)
{
	float hS = US[0];
	float uS = US[1];
	float vS = US[2];

	float hN = UN[0];
	float uN = UN[1];
	float vN = UN[2];

	std::vector<float> FRoe(3);

	getRoeFlux(US,UN,FRoe);

	// First the Euler part
	Gatface[0] = (vS + vN - FRoe[0])/2 ; //h*v;
	Gatface[1] = (hS*uS*vS + hN*uN*vN - FRoe[1])/2;//h*u*v;
	Gatface[2] = (hS*vS*vS + 0.5*9.8*hS*hS + hN*vN*vN + 0.5*9.8*hN*hN - FRoe[2])/2; //h*v*v + 0.5*9.8*h*h;
}

void update(vector<vector<vector<float> > > & U, vector<vector<vector<float> > > & UNew, float dt, float dx, float dy, int Nx, int Ny)
{
	std::vector<float> Fx(3);
	std::vector<float> Gy(3);
	
	for (int i = 1; i < Nx-1; ++i)
	{
		for (int j = 1; j < Ny-1; ++j)
		{			
			netfluxatXinterface(U[i-1][j],U[i][j],Fx);
			netfluxatXinterface(U[i][j-1],U[i][j],Gy);
			for (int k = 0; k < 3; ++k)
			{
				UNew[i-1][j][k] = (dt/dx)*Fx[k];
				UNew[i][j][k] = -((dt/dx)*Fx[k]+ (dt/dy)*Gy[k]);
				UNew[i][j-1][k] = (dt/dy)*Gy[k] ;
			}
		}
	}

	// Applying Boundary Condition
	for (int i = 0; i < Ny; ++i)
	{
		for (int c = 0; c < 3; ++c)
		{
			UNew[i][0][c] = UNew[i][1][c];	
			UNew[i][Nx-1][c] = UNew[i][Nx-2][c];	
		}
	}

	for (int i = 0; i < Nx; ++i)
	{
		for (int j = 0; j < Ny; ++j)
		{
			for (int k = 0; k < 3; ++k)
			{
				U[i][j][k] = UNew[i][j][k]; 
			}
		}
	}
}

void TestCase(vector<vector<vector<float> > > & U, vector<vector<vector<float> > > & UNew,int Nx, int Ny)
{
	std::vector<float> WL0(3);

	std::vector<float> WR0(3);

	// Case 1 : Simple moving 1D Shock, IVP 
	// Left (4) initial values
	WL0[0] = 1. ; // rho
	WL0[1] = 0. ; // u
	WL0[2] = 0. ; // v

	// Right(1) initial values
	WR0[0] = 0.5; 
	WR0[1] = 0. ; 
	WR0[2] = 0. ; 

	for (int i = 0; i < Nx; ++i)
	{
		for(int j = 0; j < Ny; ++j)
		{
			if (i<floor(Nx/2))
			{
				U[i][j][0] = WR0[0];
				U[i][j][1] = WR0[1];
				U[i][j][2] = WR0[2];

				UNew[i][j][0] = WR0[0];
				UNew[i][j][1] = WR0[1];
				UNew[i][j][2] = WR0[2];			
			}
			else
			{
				U[i][j][0] = WL0[0];
				U[i][j][1] = WL0[1];
				U[i][j][2] = WL0[2];

				UNew[i][j][0] = WL0[0];
				UNew[i][j][1] = WL0[1];
				UNew[i][j][2] = WL0[2];
			}
		}
	}
}

int main()
{
	//INPUTS = 1.4 ;
	int Nx = 100;
	int Ny = 10;
	float dx = 1.0/Nx; // change the dx in the matalb file too
	float dy = dx; // change the dx in the matalb file too

	float CFL = 0.0002;
	float time = 0.0; 
	int nSteps = 0; 
	float dt = CFL*dx;
	float tEnd = dt*10;//0.04;

	typedef vector<float> Dim1;
	typedef vector<Dim1> Dim2;
	typedef vector<Dim2> Dim3;

	Dim3 U(Nx,Dim2(Ny,Dim1(3))); // To store the conserved variables(h, hu, hv) 
	Dim3 UNew(Nx,Dim2(Ny,Dim1(3))); // To store the conserved variables(h, hu, hv) 
	
	TestCase(U,UNew,Nx,Ny); // setting the initial values
	
	
	// Solver
	while(time<tEnd)
	{
		for (int i = 0; i < Nx; ++i)
		{
			//check whether there are NaN 
			if (isnan(U[i][1][0]) == 1)
			{
				cout << "Can not Compute" << endl;
				return 0;
			}
		}
		/////////////////
		// dt = getTimeStep(CFL,dx,dy,Nx,U);
		time = time+dt;
		/////////////////
		update(U,UNew,dt,dx,dy,Nx,Ny); // This is very important step
		cout << "time   " << time << endl; 
	}
	
	//OUTPUTS IN CSV FORMET
	ofstream WriteU ;
	WriteU.open("conserved.csv");
	for (int i = 0; i < Nx; ++i)
	{
		WriteU << dx*i << "," << U[i][Ny/2][0] << "," << U[i][Ny/2][1] <<","<< U[i][Ny/2][2] << endl ;
	}
	
	ofstream WriteW ;
	WriteW.open("premetive.csv");
	for (int i = 0; i < Nx; ++i)
	{
		WriteW << dx*i << "," << U[i][Ny/2][0] << "," << U[i][Ny/2][1]/U[i][Ny/2][0] <<","<< U[i][Ny/2][2]/U[i][Ny/2][0] << endl;
	}
	return 0;
}
