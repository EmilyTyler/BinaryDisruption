#include <iostream>
#include <array>
#include "random_direction.h"
using namespace std;

int main() {
	array<double, 3> arr = randomDirection();
	for (int i=0; i<3; i++)
	{
		cout << arr[i] << endl;
	}

}

