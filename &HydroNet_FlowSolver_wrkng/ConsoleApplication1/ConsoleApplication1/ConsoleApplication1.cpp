// ConsoleApplication1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <vector>
int sum(std::vector<bool> vector) {
	int x=0;
	for (int c = 0; c < vector.size(); c++)
		x = x + vector.at(c);
	return x;
}

int main()
{
	std::vector<bool> s = { true,true,true,false,false };
	int x;

	x = sum(s);

    return 0;
}

