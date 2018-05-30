//#ifndef ALLSOME_UTILITIES_H
//#define ALLSOME_UTILITIES_H
#pragma once

//#define DEBUG 1 

#ifdef DEBUG
#define dout cout
#else
#define dout 0 && cout
#endif

#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cstring>
#include <sys/stat.h>
#include <cassert>
#include <cstdlib>

using namespace std;

/*split function*/
//std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);

/*This split function only support char as delim, string as delim please boost split function*/
std::vector<std::string> split(const std::string &, char);

inline bool FileExists(const std::string& name) {
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}

template <typename T>
inline
string ToString(const T & number) {
	string String = static_cast<ostringstream*>(&(ostringstream() << number))->str();
	return String;
}

inline void dsptime()
{
	time_t nowtime;
	//nowtime = time(NULL); //get int time number
	time(&nowtime); // get current time
	struct tm * ptm = localtime(&nowtime);  //convert time to local time
	cout << ptm->tm_mon + 1 << "/" << ptm->tm_mday << "/" << ptm->tm_year + 1900 << ",";
	cout << ptm->tm_hour << ":" << ptm->tm_min << ":" << ptm->tm_sec << " ";
}

//#endif
