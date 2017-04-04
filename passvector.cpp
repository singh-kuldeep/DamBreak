#include "iostream"
#include <vector>
#include <fstream>
#include <math.h>
using namespace std;

void vecwrite(std::vector<float> & v)
{
	cout << v[0] << v[1] << v [2] << endl;
}
int main()
{
	typedef vector<float> Dim1;
	typedef vector<Dim1> Dim2;
	typedef vector<Dim2> Dim3;

	Dim3 U(3,Dim2(3,Dim1(3)));

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			for (int k = 0; k < 3; ++k)
			{
				U[i][j][k] = i+j+k;
			}
		}
	}

	vecwrite(U[1][2]);

	return 0;
}