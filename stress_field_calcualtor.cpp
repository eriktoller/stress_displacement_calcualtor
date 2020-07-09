// This program is under progress and will be upadted

// Initaly some dummy code is included to create the file

#include <iostream>
#include <cmath>

double square(double x)
{
	double ret;
	ret = 1;
	for (int i = 0; i < 2; i++)
	{
		ret = ret * x;
	}
	return ret;
}

int main()
{
	double number;
	std::cout << "Input a number: ";
	std::cin >> number;
	double sq;
	sq = square(number);
	std::cout << "This is your give number squared: " << sq << std::endl;
	return 0;
}