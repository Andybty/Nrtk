#ifndef DPI_TYPES_BASIC_H_
#define DPI_TYPES_BASIC_H_

/*后期扩展*/
#include <list>
#include <thread>
#include <memory>
#include <mutex>
#include <iostream>
#include <string>

typedef int int32              ;    
typedef unsigned int uint32    ;
typedef short int16            ;
typedef unsigned short uint16  ;
typedef unsigned char uint8    ;
typedef signed char int8       ;


#	if defined(_WIN32_WCE)
#     if _WIN32_WCE <= 0x400
		typedef __int64 int64;
		typedef unsigned __int64 uint64;
#	endif
#else
typedef long long int64           ;
typedef unsigned long long uint64 ;
#endif

typedef struct tagTimeStamp {
	time_t tiTime             ;        /* 用 time_t类型表达的时间信息（s) */
	double dSec               ;        /* 小数秒*/
} TimeStamp;
typedef TimeStamp gtime_t;
#endif