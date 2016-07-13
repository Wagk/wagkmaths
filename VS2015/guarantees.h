#ifndef guarantees_h__
#define guarantees_h__


#include <cassert>

#ifdef _DEBUG
#define EXPECT(__expr__) assert((__expr__) == true)
#define ENSURE(__expr__) assert((__expr__) == true)
#else
#define EXPECT(__expr__) 
#define ENSURE(__expr__) 
#endif

#define ABORT() assert(false)


#endif

