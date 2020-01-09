#ifndef ERRORS_H
#define ERRORS_H

enum class Errors { OK, UNKNOWN, FRmapsReturnNULL };

struct LastError { 
	static Errors lasterror; 
};

#endif