// ConsoleApplication2.cpp : Defines the entry point for the console application.
//

#include "ConsoleApplication2.h"
#include "stdafx.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <fstream>




template < class T>
std::vector<int> find_equal(std::vector<T> A, std::vector<T> B) {			//returns a vector whith the indexes that are equal in both vectos
	std::vector<int> C;
	for (int i = 0; i < A.size(); i++) {
		if (B.size() > 1 && A.at(i) == B.at(i)) 
			C.push_back(i);
		
		else if(A.at(i) == B.at(0))
			C.push_back(i);

	}
	return C;
}

template < class T>
std::vector<int> find_not_equal(std::vector<T> A, std::vector<T> B){		//returns a vector whith the indexes that are not equal in both vectos
	std::vector<int> C;
		for (int i = 0;i < A.size() ; i++) {
			if (B.size() > 1 && A.at(i) != B.at(i))
				C.push_back(i);
			else if (A.at(i)!= B.at(0))
				C.push_back(i);
		}
	return C;
}

template < class T>
std::vector<int> find_less(std::vector<T> A, std::vector<T> B){			//returns a vector whith the indexes that are lesser in B than A 
	std::vector<int> C;
		for (int i = 0; i < A.size() ; i++) {
			if (B.size() > 1 && A.at(i) < B.at(i))
				C.push_back(i);
			else if (A.at(i) < B.at(0))
				C.push_back(i);
		}
	return C;
}

template < class T>
std::vector<int> find_grater(std::vector<T> A, std::vector<T> B){			//returns a vector whith the indexes that are grater in B than A 
	std::vector<int> C;
		for (int i = 0; i < A.size() ; i++) {
			if (B.size() > 1 && A.at(i) > B.at(i))
				C.push_back(i);
			else if (A.at(i) > B.at(0))
				C.push_back(i);
		}
	return C;
}

template < class T>
std::vector<int> find_less_equal(std::vector<T> A, std::vector<T> B){		//returns a vector whith the indexes that are equal or lesser in B than A 
	std::vector<int> C;
		for (int i = 0; i < A.size() ; i++) {
			if (B.size() > 1 && A.at(i) <= B.at(i))
				C.push_back(i);
			else if (A.at(i) <= B.at(0))
				C.push_back(i);
		}
	return C;
}
template < class T>
std::vector<int> find_grater_equal(std::vector<T> A, std::vector<T> B){	//returns a vector whith the indexes that are equal or grater in B than A
		std::vector<int> C;
			for (int i = 0; i < A.size() ; i++) {
				if (B.size() > 1 && A.at(i) >= B.at(i))
					C.push_back(i);
				else if (A.at(i) >= B.at(0))
					C.push_back(i);
			}
		return C;
}



int main()
{
	std::vector<int> a;
	//std::vector<int>  ;
	int b,c;
	
	a= { 3,52,7 };
	
	b = 52;
	
	
	//c = find_equal(a, { b });
	c = find_equal(a, { b }).at(0);
	/*c = find_less(a, b);
	c = find_grater(a, b);
	c = find_less_equal(a, b);
	c = find_grater_equal(a, b);*/

	

	
    return 0;
}

