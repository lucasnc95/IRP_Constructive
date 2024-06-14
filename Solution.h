#ifndef SOLUTION_H
#define SOLUTION_H
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <climits>
#include "IRP.h"
#include "Location.h"



class Solution {



bool Check_Viability(const std::vector<std::vector<int>> routes, const std::vector<std::vector<int>> inventory );

void Print_Solution(const std::vector<std::vector<int>> routes, const std::vector<std::vector<int>> inventory );







};

#endif