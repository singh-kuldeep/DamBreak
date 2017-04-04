#include "iostream"
#include <vector>
#include <fstream>
#include <math.h>
using namespace std ;
void assign(std::vector<float> & v, std::vector<float> w)
{
	int sz = v.size();
	for (int i = 0; i < sz; ++i)
	{
		v[i] = 10.;
		// std::cout << v[i] << std::endl;
	}
}

int main()
{


	std::vector< std::vector<float>  > v(4, std::vector<float> (3));

	std::vector<float> x;
	x.resize(3);
	x[0] = 1 ;
	x[1] = 2 ;
	x[2] = 3 ;

assign(v[2],x);

for (int i = 0; i < 4; ++i)
{
	for (int j = 0; j < 3; ++j)
	{
		// v[i][j] = i*j ; 

		cout << v[i][j]  << endl;
	}
}

return 0;
}