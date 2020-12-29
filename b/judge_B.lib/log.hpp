#ifndef HITACHI_HOKUDAI_2020_SRC_UTILS_LOG_HPP
#define HITACHI_HOKUDAI_2020_SRC_UTILS_LOG_HPP

#include <cstdio>

#ifndef DEBUG_LEVEL
#define DEBUG_LEVEL 0
#endif

#if     DEBUG_LEVEL >= 1
#define ERROR(...) do{fprintf(stderr, "[%s:% 4d] ERROR  : ", __FILE__,  __LINE__); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n");}while(0);
#else
#define ERROR(fmt, ...)
#endif

#if     DEBUG_LEVEL >= 2
#define WARNING(...) do{fprintf(stderr, "[%s:% 4d] WARNING : ", __FILE__,  __LINE__); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n");}while(0);
#else
#define WARNING(fmt, ...)
#endif

#if     DEBUG_LEVEL >= 3
#define INFO(...) do{fprintf(stderr, "[%s:% 4d] INFO   : ", __FILE__,  __LINE__); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n");}while(0);
#else
#define INFO(fmt, ...)
#endif

#if     DEBUG_LEVEL >= 4
#define DEBUG(...) do{fprintf(stderr, "[%s:% 4d] DEBUG  : ", __FILE__,  __LINE__); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n");}while(0);
#else
#define DEBUG(fmt, ...)
#endif

#endif
