/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#if !defined(ASSERT_HPP)
  #define ASSERT_HPP

  #if defined(_DEBUG)
    void my_assert(const char* Condition, const char* FileName, int LineNumber);
    #define ASSERT(E) if (!(E))  my_assert(#E, __FILE__, __LINE__);
    #define VERIFY(E) if (!(E))  my_assert(#E, __FILE__, __LINE__);
  #else
    #define ASSERT(E)
    #define VERIFY(E) E;
  #endif

#endif

#ifndef __STRING1
#define __STRING1(__WHAT) #__WHAT
#endif
#ifndef __STRING
#define __STRING(__WHAT) __STRING1(__WHAT)
#endif


#ifndef __MEMOK
#define __MEMOK(_WHAT) if ( ! _WHAT ){ char msg[80]=" Malloc returned NULL for '"; (void) strcat(msg, __STRING(_WHAT)) ; (void) strcat(msg,"'");\
                                         NAMD_die ( (const char*) msg );\
                                     }
#endif
